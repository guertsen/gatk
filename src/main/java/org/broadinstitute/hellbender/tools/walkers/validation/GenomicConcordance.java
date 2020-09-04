package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.sam.util.Pair;

import java.io.File;

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

    public static final String TRUTH_BLOCK_HISTOGRAM_LONG_NAME = "truth-block-histogram";
    public static final String TRUTH_BLOCK_HISTOGRAM_SHORT_NAME = "tbh";

    public static final String EVAL_BLOCK_HISTOGRAM_LONG_NAME = "eval-block-histogram";
    public static final String EVAL_BLOCK_HISTOGRAM_SHORT_NAME = "ebh";

    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME = "confidence-concordance-histogram";
    public static final String CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME = "cch";

    @Argument(doc = "A histogram of block lengths and their associated confidence scores for the truth sample",
            fullName = TRUTH_BLOCK_HISTOGRAM_LONG_NAME,
            shortName = TRUTH_BLOCK_HISTOGRAM_SHORT_NAME)
    protected File truthBlockHistogramFile;

    @Argument(doc = "A histogram of block lengths and their associated confidence scores for the eval sample",
            fullName = EVAL_BLOCK_HISTOGRAM_LONG_NAME,
            shortName = EVAL_BLOCK_HISTOGRAM_SHORT_NAME)
    protected File evalBlockHistogramFile;

    @Argument(doc = "A table of combined reference block lengths and confidence ???TODO",
            fullName = CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME,
            shortName = CONFIDENCE_CONCORDANCE_HISTOGRAM_SHORT_NAME)
    protected File confidenceConcordanceHistogramFile;

    // TODO this should be a Histogram<Pair<Integer, Integer>>, however, the MetricsFile class cannot read
    // arbitrary types, therefore, it must be converted to a String, which is probably much slower
    @VisibleForTesting
    final Histogram<String> truthBlockHistogram = new Histogram<>();
    @VisibleForTesting
    final Histogram<String> evalBlockHistogram = new Histogram<>();
    @VisibleForTesting
    final Histogram<String> confidenceConcordanceHistogram = new Histogram<>();

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
        if (currentTruthVariantContext != null && currentEvalVariantContext != null) {
            int blockStart = Math.max(currentTruthVariantContext.getStart(), currentEvalVariantContext.getStart());
            int blockEnd = Math.min(currentTruthVariantContext.getEnd(), currentEvalVariantContext.getEnd());
            int jointBlockLength = blockEnd - blockStart + 1;
            if (jointBlockLength > 0) {
                confidenceConcordanceHistogram.increment(new Pair<>(currentTruthVariantContext.getGenotype(0).getGQ(), currentEvalVariantContext.getGenotype(0).getGQ()).toString(), jointBlockLength);
            }
        }

        currentTruthVariantContext = null;
        currentEvalVariantContext = null;
        currentContig = null;
    }

    private void evaluateNewContig(TruthVersusEval truthVersusEval) {
        // If not beginning of file
        if (currentContig != null) {
            evaluateEndOfContig();
        }

        currentContig = truthVersusEval.getTruthIfPresentElseEval().getContig();
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
                confidenceConcordanceHistogram.increment(new Pair<>(currentTruthVariantContext.getGenotype(0).getGQ(), currentEvalVariantContext.getGenotype(0).getGQ()).toString(), blockEnd - blockStart + 1);
            }

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
            int blockLength = truthVersusEval.getTruth().getEnd() - truthVersusEval.getTruth().getStart() + 1;

            // TODO can a non_ref block ever have a number of genotypes != 1?
            if(truthVersusEval.getTruth().getGenotypes().size() != 1) {
                throw new IllegalStateException(String.format("The NON_REF block at %s has more than one genotype, which is not supported.", truthVersusEval.getTruth().toStringDecodeGenotypes()));
            }
            Genotype genotype = truthVersusEval.getTruth().getGenotype(0);
            int gq = genotype.getGQ();
            truthBlockHistogram.increment(new Pair<>(blockLength, gq).toString());
        }

        // Eval
        if (truthVersusEval.hasEval() && isNonRef(truthVersusEval.getEval())) {
            currentEvalVariantContext = truthVersusEval.getEval();

            // The end is inclusive, thus the plus one when calculating the length
            int blockLength = truthVersusEval.getEval().getEnd() - truthVersusEval.getEval().getStart() + 1;

            // TODO can a non_ref block ever have a number of genotypes != 1?
            if(truthVersusEval.getEval().getGenotypes().size() != 1) {
                throw new IllegalStateException(String.format("The NON_REF block at %s has more than one genotype, which is not supported.", truthVersusEval.getEval().toStringDecodeGenotypes()));
            }
            Genotype genotype = truthVersusEval.getEval().getGenotype(0);
            int gq = genotype.getGQ();
            evalBlockHistogram.increment(new Pair<>(blockLength, gq).toString());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

        evaluateEndOfContig();

        // Truth histogram
        MetricsFile<?, String> metricsFile = getMetricsFile();
        metricsFile.addHistogram(truthBlockHistogram);
        metricsFile.write(truthBlockHistogramFile);

        // Eval histogram
        metricsFile = getMetricsFile();
        metricsFile.addHistogram(evalBlockHistogram);
        metricsFile.write(evalBlockHistogramFile);

        // Confidence concordance
        metricsFile = getMetricsFile();
        metricsFile.addHistogram(confidenceConcordanceHistogram);
        metricsFile.write(confidenceConcordanceHistogramFile);

        return "SUCCESS";
    }
}