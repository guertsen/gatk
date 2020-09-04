package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.testutils.BaseTest;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.Pair;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by Michael Gatzen on 8/26/20.
 */
public class GenomicConcordanceUnitTest extends CommandLineProgramTest {
    static class TestGenomicConcordanceVariant {
        private final String altAllele;
        private final int start;
        private final int stop;
        private final int confidence;

        public TestGenomicConcordanceVariant(final String altAllele, final int start, final int stop, final int confidence) {
            this.altAllele = altAllele;
            this.start = start;
            this.stop = stop;
            this.confidence = confidence;
        }

        public String getAltAllele() { return altAllele; }
        public int getStart() { return start; }
        public int getStop() { return stop; }
        public int getConfidence() { return confidence; }
    }

    private Pair<File, File> writeTestGVCFs(List<TestGenomicConcordanceVariant> truthVariants, List<TestGenomicConcordanceVariant> evalVariants) throws Exception {
        File truthFile = createTempFile("truth", ".gvcf");
        FileWriter writer = new FileWriter(truthFile);
        writer.write("##fileformat=VCFv4.2\n");
        writer.write("##contig=<ID=test_contig,length=1000>\n");
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n");
        for (TestGenomicConcordanceVariant variant : truthVariants) {
            writer.write(String.format("test_contig\t%s\t.\tA\t%s\t%s\t.\t%s\tGT:GQ\t%s:%s\n",
                    variant.start,
                    variant.getAltAllele(),
                    variant.getAltAllele().equals("<NON_REF>") ? "." : "1",
                    variant.getAltAllele().equals("<NON_REF>") ? String.format("END=%s", variant.getStop()) : ".",
                    variant.getAltAllele().equals("<NON_REF>") ? "0/0" : "0/1",
                    variant.getConfidence()));
        }
        writer.close();

        File evalFile = createTempFile("eval", ".gvcf");
        writer = new FileWriter(evalFile);
        writer.write("##fileformat=VCFv4.2\n");
        writer.write("##contig=<ID=test_contig,length=1000>\n");
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n");
        for (TestGenomicConcordanceVariant variant : evalVariants) {
            writer.write(String.format("test_contig\t%s\t.\tA\t%s\t%s\t.\t%s\tGT:GQ\t%s:%s\n",
                    variant.start,
                    variant.getAltAllele(),
                    variant.getAltAllele().equals("<NON_REF>") ? "." : "1",
                    variant.getAltAllele().equals("<NON_REF>") ? String.format("END=%s", variant.getStop()) : ".",
                    variant.getAltAllele().equals("<NON_REF>") ? "0/0" : "0/1",
                    variant.getConfidence()));
        }
        writer.close();

        return new Pair<>(truthFile, evalFile);
    }

    @DataProvider
    public Object[][] provideEvaluations() {
        return new Object[][] {
                // No non_ref blocks
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("C", 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {

                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {

                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Two non_ref blocks
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with gap then variant
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("C", 5, 5, 99),
                                new TestGenomicConcordanceVariant("<NON_REF>", 6, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "5,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 5}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with gap then non_ref
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 6, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "5,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 5}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with single block
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 3, 6, 98),
                                new TestGenomicConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "4,98", 1},
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 10}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // No overlap at all
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant("<NON_REF>", 3, 6, 98)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "4,98", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
        };
    }

    @Test(dataProvider = "provideEvaluations")
    public void test(final List<TestGenomicConcordanceVariant> truthVariants, final List<TestGenomicConcordanceVariant> evalVariants, final Map<String, Integer> expectedTruthBlockHistogram, final Map<String, Integer> expectedEvalBlockHistogram, final Map<String, Integer> expectedConfidenceConcordance) throws Exception {
        Pair<File, File> inputFiles = writeTestGVCFs(truthVariants, evalVariants);

        final Path summary = createTempPath("summary", ".txt");
        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, inputFiles.getLeft().toString(),
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, inputFiles.getRight().toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
                "--" + GenomicConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + GenomicConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + GenomicConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        if (truthBlockMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
            Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            truthBlockHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedTruthBlockHistogram.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedTruthBlockHistogram.get(bin.getId()));
                expectedTruthBlockHistogram.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
            Assert.assertEquals(expectedTruthBlockHistogram.size(), 0);
        } else {
            Assert.assertEquals(expectedTruthBlockHistogram.size(), 0);
        }

        MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));
        if (evalBlockMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);
            Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            evalBlockHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedEvalBlockHistogram.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedEvalBlockHistogram.get(bin.getId()));
                expectedEvalBlockHistogram.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
            Assert.assertEquals(expectedEvalBlockHistogram.size(), 0);
        } else {
            Assert.assertEquals(expectedEvalBlockHistogram.size(), 0);
        }

        MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        if (confidenceConcordanceMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
            Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            confidenceConcordanceHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedConfidenceConcordance.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedConfidenceConcordance.get(bin.getId()));
                expectedConfidenceConcordance.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
            Assert.assertEquals(expectedConfidenceConcordance.size(), 0);
        } else {
            Assert.assertEquals(expectedConfidenceConcordance.size(), 0);
        }

        return;
    }
}
