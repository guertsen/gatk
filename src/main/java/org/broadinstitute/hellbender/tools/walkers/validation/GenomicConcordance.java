package org.broadinstitute.hellbender.tools.walkers.validation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.lang.mutable.MutableLong;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

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

    @Argument(doc = "A table of reference confidence ???TODO",
            fullName = CONFIDENCE_HISTOGRAM_LONG_NAME,
            shortName = CONFIDENCE_HISTOGRAM_SHORT_NAME)
    protected File confidenceHistogramFile;

    @Argument(doc = "A table of reference block lengths ???TODO",
            fullName = BLOCK_LENGTH_HISTOGRAM_LONG_NAME,
            shortName = BLOCK_LENGTH_HISTOGRAM_SHORT_NAME)
    protected File blockLengthHistogramFile;

    @VisibleForTesting
    final SortedMap<Long, GenomicConcordanceHistogramEntry> blockLengthHistogram = new TreeMap<>();
    @VisibleForTesting
    final SortedMap<Long, GenomicConcordanceHistogramEntry> confidenceHistogram = new TreeMap<>();

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        // Explicitly allow symbolic variants
        return vc -> !vc.isFiltered() && !vc.isStructuralIndel();
    }

    @Override
    protected void apply(TruthVersusEval truthVersusEval, ReadsContext readsContext, ReferenceContext refContext) {
        super.apply(truthVersusEval, readsContext, refContext);

        if (truthVersusEval.hasTruth() && truthVersusEval.getTruth().isSymbolic()) {
            // TODO get length on reference or just end-start?
            // The end is inclusive, thus the plus one when calculating the length
            long blockLength = truthVersusEval.getTruth().getEnd() - truthVersusEval.getTruth().getStart() + 1;
            blockLengthHistogram.putIfAbsent(blockLength, new GenomicConcordanceHistogramEntry(blockLength));
            blockLengthHistogram.get(blockLength).incrementTruthValue();

            truthVersusEval.getTruth().getGenotypes().forEach(g -> {
                long gq = g.getGQ();
                confidenceHistogram.putIfAbsent(gq, new GenomicConcordanceHistogramEntry(gq));
                confidenceHistogram.get(gq).incrementTruthValue();
            });
        }

        if (truthVersusEval.hasEval() && truthVersusEval.getEval().isSymbolic()) {
            // The end is inclusive, thus the plus one when calculating the length
            long blockLength = truthVersusEval.getEval().getEnd() - truthVersusEval.getEval().getStart() + 1;
            blockLengthHistogram.putIfAbsent(blockLength, new GenomicConcordanceHistogramEntry(blockLength));
            blockLengthHistogram.get(blockLength).incrementEvalValue();

            truthVersusEval.getEval().getGenotypes().forEach(g -> {
                long gq = g.getGQ();
                confidenceHistogram.putIfAbsent(gq, new GenomicConcordanceHistogramEntry(gq));
                confidenceHistogram.get(gq).incrementEvalValue();
            });
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

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

        return "SUCCESS";
    }
}
