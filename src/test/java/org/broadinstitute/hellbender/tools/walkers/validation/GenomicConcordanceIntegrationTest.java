package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by Takuto Sato on 1/31/17.
 */
public class GenomicConcordanceIntegrationTest extends CommandLineProgramTest{
    final double epsilon = 1e-3;

    private static final String CONCORDANCE_TEST_DIR = toolsTestDir + "concordance/";

    @Test
    public void testMinimal() throws Exception {
        final File truthVcf = new File(CONCORDANCE_TEST_DIR, "testGVCF_truth.g.vcf");
        final File evalVcf = new File(CONCORDANCE_TEST_DIR, "testGVCF_eval.g.vcf");
        final Path summary = createTempPath("summary", ".txt");
        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
                "--" + GenomicConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + GenomicConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + GenomicConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        // TODO test
    }

    @Test
    public void testIdentical() throws Exception {
        final File truthVcf = new File(CONCORDANCE_TEST_DIR, "expected.testGVCFMode.gatk4.g.vcf");
        final File evalVcf = new File(CONCORDANCE_TEST_DIR, "expected.testGVCFMode.gatk4.g.vcf");
        final Path summary = createTempPath("summary", ".txt");
        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
                "--" + GenomicConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + GenomicConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + GenomicConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        // TODO test
    }

    @Test
    public void customTestTODODelete() throws Exception {
        final File truthVcf = new File("/home/broad/tmp_data", "CHMI_CHMI3_WGS1.g.vcf.gz");
        final File evalVcf = new File("/home/broad/tmp_data", "CHMI_CHMI3_WGS1gvcf.hard-filtered.gvcf.gz");
        final File summary = new File("/home/broad/Desktop", "summary.txt");
        final File truthBlockHistogramFile = new File("/home/broad/Desktop", "truth_block_histogram.tsv");
        final File evalBlockHistogramFile = new File("/home/broad/Desktop", "eval_block_histogram.tsv");
        final File confidenceConcordanceHistogramFile = new File("/home/broad/Desktop", "confidence_concordance_histogram.tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
                "--" + GenomicConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + GenomicConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + GenomicConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        // TODO test
    }
}
