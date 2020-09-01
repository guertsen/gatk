package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import htsjdk.samtools.util.Histogram;
import picard.sam.util.Pair;

import java.io.File;
import java.io.IOException;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * ??? TODO
 */

@CommandLineProgramProperties(
        summary = GenomicConcordance.USAGE_SUMMARY,
        oneLineSummary = GenomicConcordance.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
public class GenomicConcordance extends Concordance {
    static final String USAGE_ONE_LINE_SUMMARY = "???TODO";
    static final String USAGE_SUMMARY = "???TODO";

    public static final String CONFIDENCE_HISTOGRAM_LONG_NAME = "confidence-histogram";
    public static final String CONFIDENCE_HISTOGRAM_SHORT_NAME = "ch";

    public static final String BLOCK_LENGTH_HISTOGRAM_LONG_NAME = "block-length-histogram";
    public static final String BLOCK_LENGTH_HISTOGRAM_SHORT_NAME = "blh";

    public static final String BLOCK_LENGTH_AND_CONFIDENCE_HISTOGRAM_LONG_NAME = "block-length-and-confidence-histogram";
    public static final String BLOCK_LENGTH_AND_CONFIDENCE_HISTOGRAM_SHORT_NAME = "blch";

    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME = "confidence-concordance-histogram";
    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME = "cch";

    @Argument(doc = "A table of reference confidence ???TODO",
            fullName = CONFIDENCE_HISTOGRAM_LONG_NAME,
            shortName = CONFIDENCE_HISTOGRAM_SHORT_NAME)
    protected File confidenceHistogramFile;

    @Argument(doc = "A table of reference block lengths ???TODO",
            fullName = BLOCK_LENGTH_HISTOGRAM_LONG_NAME,
            shortName = BLOCK_LENGTH_HISTOGRAM_SHORT_NAME)
    protected File blockLengthHistogramFile;

    @Argument(doc = "A table of combined reference block lengths and confidence ???TODO",
            fullName = BLOCK_LENGTH_AND_CONFIDENCE_HISTOGRAM_LONG_NAME,
            shortName = BLOCK_LENGTH_AND_CONFIDENCE_HISTOGRAM_SHORT_NAME)
    protected File blockLengthAndConfidenceHistogramFile;

    @Argument(doc = "A table of combined reference block lengths and confidence ???TODO",
            fullName = CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME,
            shortName = CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME)
    protected File confidenceConcordanceHistogramFile;

    @VisibleForTesting
    final SortedMap<String, GenomicConcordanceHistogramEntry> blockLengthHistogram = new TreeMap<>();
    @VisibleForTesting
    final SortedMap<String, GenomicConcordanceHistogramEntry> confidenceHistogram = new TreeMap<>();
    @VisibleForTesting
    final SortedMap<String, GenomicConcordanceHistogramEntry> blockLengthAndConfidenceHistogram = new TreeMap<>();
    @VisibleForTesting
    final Histogram<Pair<Integer, Integer>> confidenceConcordanceHistogram = new Histogram<>();

    private VariantContext currentTruthVariantContext = null;
    private VariantContext currentEvalVariantContext = null;
    private String currentContig = null;

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        // Explicitly allow symbolic variants
        return vc -> !vc.isFiltered() && !vc.isStructuralIndel();
    }

    private boolean isNonRef(VariantContext variantContext) {
        return variantContext.isSymbolic() && variantContext.getAlternateAllele(0).isNonRefAllele();
    }

    private void evaluateEndOfContig() {
        if (currentTruthVariantContext.getEnd() != currentEvalVariantContext.getEnd()) {
            logger.warn(String.format("For the contig '%s', the end of the last truth variant (%d) does not equal the end of the last eval variant (%d).", currentContig, currentTruthVariantContext.getEnd(), currentEvalVariantContext.getEnd()));
        }
        int blockStart = Math.max(currentTruthVariantContext.getStart(), currentEvalVariantContext.getStart());
        int blockEnd = Math.min(currentTruthVariantContext.getEnd(), currentEvalVariantContext.getEnd());
        int jointBlockLength = blockEnd - blockStart + 1;
        confidenceConcordanceHistogram.increment(new Pair<>(currentTruthVariantContext.getGenotype(0).getGQ(), currentEvalVariantContext.getGenotype(0).getGQ()), jointBlockLength);

        currentTruthVariantContext = null;
        currentEvalVariantContext = null;
        currentContig = null;
    }

    private void evaluateNewContig(TruthVersusEval truthVersusEval) {
        // If not beginning of file
        if (currentContig != null) {
            evaluateEndOfContig();
        }

        // Both VCFs must have a (at least symolic) variant here
        if (!(truthVersusEval.hasTruth() && truthVersusEval.hasEval())) {
            throw new IllegalStateException(String.format("Beginning of new contig (%s) and %s VCF has no variant.", truthVersusEval.getTruthIfPresentElseEval().getContig(), truthVersusEval.hasTruth() ? "eval" : "truth"));
        }
        currentContig = truthVersusEval.getTruth().getContig();
    }

    @Override
    protected void apply(TruthVersusEval truthVersusEval, ReadsContext readsContext, ReferenceContext refContext) {
        // TODO is this right?
        // If a truth variant is present, first evaluate the base class filter before presenting it to the base class
        if (!truthVersusEval.hasTruth() || super.makeTruthVariantFilter().evaluate(truthVersusEval.getTruth())) {
            super.apply(truthVersusEval, readsContext, refContext);
        }

        // New contig or beginning of file
        if (!truthVersusEval.getTruthIfPresentElseEval().getContig().equals(currentContig)) {
            evaluateNewContig(truthVersusEval);
        }

        // Evaluate only when currently seeing two NON_REF blocks
        if (currentTruthVariantContext != null && currentEvalVariantContext != null) {
            int blockStart = Math.max(currentTruthVariantContext.getStart(), currentEvalVariantContext.getStart());
            int blockEnd = Math.min(currentTruthVariantContext.getEnd(), currentEvalVariantContext.getEnd());
            int jointBlockLength = blockEnd - blockStart + 1;
            // It is possible that jointBlockLength is negative if there is a gap in one file and the start of a new block in the other file.
            // Since there is no overlap though, we can just skip that case.
            if (jointBlockLength > 0) {
                confidenceConcordanceHistogram.increment(new Pair<>(currentTruthVariantContext.getGenotype(0).getGQ(), currentEvalVariantContext.getGenotype(0).getGQ()), blockEnd - blockStart + 1);
            }
//            if (currentPosition - 1 != Math.min(currentTruthVariantContext.getEnd(), currentEvalVariantContext.getEnd())) {
//                //logger.warn(String.format("The current position (%s:%s) does not coincide with the end of a previous NON_REF block.", currentContig, currentPosition));
//            }

            int currentPosition = truthVersusEval.getTruthIfPresentElseEval().getStart();
            if (truthVersusEval.hasTruth() || currentPosition >= currentTruthVariantContext.getEnd()) {
                currentTruthVariantContext = null;
            }
            if (truthVersusEval.hasEval() || currentPosition >= currentEvalVariantContext.getEnd()) {
                currentEvalVariantContext = null;
            }
        }

        // Truth
        if (truthVersusEval.hasTruth() && isNonRef(truthVersusEval.getTruth())) {
            currentTruthVariantContext = truthVersusEval.getTruth();

            // TODO get length on reference or just end-start? Shouldn't matter for NON_REF?
            // The end is inclusive, thus the plus one when calculating the length
            String blockLength = String.valueOf(truthVersusEval.getTruth().getEnd() - truthVersusEval.getTruth().getStart() + 1);
            blockLengthHistogram.putIfAbsent(blockLength, new GenomicConcordanceHistogramEntry(blockLength));
            blockLengthHistogram.get(blockLength).incrementTruthValue();

            // TODO can a non_ref block ever have a number of genotypes != 1?
            if(truthVersusEval.getTruth().getGenotypes().size() != 1) {
                throw new IllegalStateException();//"The NON_REF block at ".append(truthVersusEval.getTruth().toStringDecodeGenotypes()) + " has more than one genotype, which is not supported.");
            }
            Genotype genotype = truthVersusEval.getTruth().getGenotype(0);
            String gq = String.valueOf(genotype.getGQ());
            confidenceHistogram.putIfAbsent(gq, new GenomicConcordanceHistogramEntry(gq));
            confidenceHistogram.get(gq).incrementTruthValue();

            blockLengthAndConfidenceHistogram.putIfAbsent(blockLength + "," + gq, new GenomicConcordanceHistogramEntry(blockLength + "," + gq));
            blockLengthAndConfidenceHistogram.get(blockLength + "," + gq).incrementTruthValue();
        }

        // Eval
        if (truthVersusEval.hasEval() && isNonRef(truthVersusEval.getEval())) {
            currentEvalVariantContext = truthVersusEval.getEval();

            // The end is inclusive, thus the plus one when calculating the length
            String blockLength = String.valueOf(truthVersusEval.getEval().getEnd() - truthVersusEval.getEval().getStart() + 1);
            blockLengthHistogram.putIfAbsent(blockLength, new GenomicConcordanceHistogramEntry(blockLength));
            blockLengthHistogram.get(blockLength).incrementEvalValue();

            // TODO can a non_ref block ever have a number of genotypes != 1?
            if(truthVersusEval.getEval().getGenotypes().size() != 1) {
                throw new IllegalStateException();//"The NON_REF block at ".append(truthVersusEval.getEval().toStringDecodeGenotypes()) + " has more than one genotype, which is not supported.");
            }
            Genotype genotype = truthVersusEval.getEval().getGenotype(0);
            String gq = String.valueOf(genotype.getGQ());
            confidenceHistogram.putIfAbsent(gq, new GenomicConcordanceHistogramEntry(gq));
            confidenceHistogram.get(gq).incrementEvalValue();

            blockLengthAndConfidenceHistogram.putIfAbsent(blockLength + "," + gq, new GenomicConcordanceHistogramEntry(blockLength + "," + gq));
            blockLengthAndConfidenceHistogram.get(blockLength + "," + gq).incrementEvalValue();
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

        evaluateEndOfContig();

        try(GenomicConcordanceHistogramEntry.Writer blockLengthHistogramWriter = GenomicConcordanceHistogramEntry.getWriter(blockLengthHistogramFile)) {
            blockLengthHistogramWriter.writeAllRecords(blockLengthHistogram.values());
        } catch (IOException e) {
            throw new UserException("Encountered an IO exception writing the block length histogram table", e);
        }

        try(GenomicConcordanceHistogramEntry.Writer confidenceHistogramWriter = GenomicConcordanceHistogramEntry.getWriter(confidenceHistogramFile)) {
            confidenceHistogramWriter.writeAllRecords(confidenceHistogram.values());
        } catch (IOException e) {
            throw new UserException("Encountered an IO exception writing the confidence histogram table", e);
        }

        try(GenomicConcordanceHistogramEntry.Writer confidenceHistogramWriter = GenomicConcordanceHistogramEntry.getWriter(blockLengthAndConfidenceHistogramFile)) {
            confidenceHistogramWriter.writeAllRecords(blockLengthAndConfidenceHistogram.values());
        } catch (IOException e) {
            throw new UserException("Encountered an IO exception writing the combined block length and confidence histogram table", e);
        }

        MetricsFile<?, Pair<Integer, Integer>> metricsFile = getMetricsFile();
        metricsFile.addHistogram(confidenceConcordanceHistogram);
        metricsFile.write(confidenceConcordanceHistogramFile);

        return "SUCCESS";
    }
}
