package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

/**
 * Created by Michael Gatzen on 8/26/20.
 */
public class GenomicConcordanceUnitTest {
    static class TestGenomicConcordanceVariant {
        private final List<String> alleles;
        private final int start;
        private final int stop;
        private final Optional<Integer> confidence;

        public TestGenomicConcordanceVariant(final List<String> alleles, final int start, final int stop, final Integer confidence) {
            this.alleles = alleles;
            this.start = start;
            this.stop = stop;
            this.confidence = confidence == null ? Optional.empty() : Optional.of(confidence);
        }

        public List<String> getAlleles() { return alleles; }
        public int getStart() { return start; }
        public int getStop() { return stop; }
        public boolean hasConfidence() { return confidence.isPresent(); }
        public int getConfidence() {
            if (confidence.isPresent()) {
                return confidence.get();
            } else {
                throw new RuntimeException("No confidence score present.");
            }
        }
    }

    @DataProvider
    public Object[][] provideEvaluations() {
        return new Object[][] {
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "C"), 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "C"), 1, 1, null)
                        ),
                        // Expected block length histogram
                        new ArrayList<GenomicConcordanceHistogramEntry>(),
                        // Expected confidence score histogram
                        new ArrayList<GenomicConcordanceHistogramEntry>(),
                },
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "<NON_REF>"), 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "C"), 1, 1, null)
                        ),
                        // Expected block length histogram
                        Arrays.asList(
                                new GenomicConcordanceHistogramEntry(1, 1, 0)
                        ),
                        // Expected confidence score histogram
                        Arrays.asList(
                                new GenomicConcordanceHistogramEntry(99, 1, 0)
                        ),
                },
                {
                        // Truth variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "<NON_REF>"), 1, 10, 99),
                                null
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "<NON_REF>"), 1, 5, 99),
                                new TestGenomicConcordanceVariant(Arrays.asList("A", "<NON_REF>"), 6, 10, 98)
                        ),
                        // Expected block length histogram
                        Arrays.asList(
                                new GenomicConcordanceHistogramEntry(10, 1, 0),
                                new GenomicConcordanceHistogramEntry(5, 0, 2)
                        ),
                        // Expected confidence score histogram
                        Arrays.asList(
                                new GenomicConcordanceHistogramEntry(99, 1, 1),
                                new GenomicConcordanceHistogramEntry(98, 0, 1)
                        ),
                },
        };
    }

    @Test(dataProvider = "provideEvaluations")
    public void testApply(final List<TestGenomicConcordanceVariant> truthVariants, final List<TestGenomicConcordanceVariant> evalVariants, final List<GenomicConcordanceHistogramEntry> expectedBlockLengthHistogram, final List<GenomicConcordanceHistogramEntry> expectedConfidenceHistogram) {
        List<AbstractConcordanceWalker.TruthVersusEval> evaluations = new ArrayList<>();

        for(int i = 0; i < truthVariants.size(); i++) {
            // Construct truth variant context object
            TestGenomicConcordanceVariant truthVariant = truthVariants.get(i);
            VariantContext truthVC = null;
            if(truthVariant != null) {
                List<Allele> truthAlleles = new ArrayList<>();
                // Add ref allele. There is no function to create a reference allele
                // from a string, hence the byte conversion. TODO create one?
                truthAlleles.add(Allele.create(truthVariant.getAlleles().get(0).getBytes(), true));
                ;

                // Use standard for loop to start at index 1
                for (int j = 1; j < truthVariant.getAlleles().size(); j++) {
                    truthAlleles.add(Allele.create(truthVariant.getAlleles().get(j)));
                }

                VariantContextBuilder builder = new VariantContextBuilder("", "1", truthVariant.start, truthVariant.stop, truthAlleles);
                GenotypeBuilder genotypeBuilder = new GenotypeBuilder("", truthAlleles);
                if (truthVariant.hasConfidence()) {
                    genotypeBuilder.GQ(truthVariant.getConfidence());
                }
                builder.genotypes(genotypeBuilder.make());
                truthVC = builder.make();
            }

            // Construct eval variant context object
            TestGenomicConcordanceVariant evalVariant = evalVariants.get(i);
            VariantContext evalVC = null;
            if(evalVariant != null) {
                List<Allele> evalAlleles = new ArrayList<>();
                evalAlleles.add(Allele.create(evalVariant.getAlleles().get(0).getBytes(), true));

                for (int j = 1; j < evalVariant.getAlleles().size(); j++) {
                    evalAlleles.add(Allele.create(evalVariant.getAlleles().get(j)));
                }

                VariantContextBuilder builder = new VariantContextBuilder("", "1", evalVariant.start, evalVariant.stop, evalAlleles);
                GenotypeBuilder genotypeBuilder = new GenotypeBuilder("", evalAlleles);
                if (evalVariant.hasConfidence()) {
                    genotypeBuilder.GQ(evalVariant.getConfidence());
                }
                builder.genotypes(genotypeBuilder.make());
                evalVC = builder.make();
            }

            // Select the right type of concordance, based on which variants have been supplied
            if(truthVariant == null) {
                evaluations.add(AbstractConcordanceWalker.TruthVersusEval.falsePositive(evalVC));
            } else if (evalVariant == null) {
                evaluations.add(AbstractConcordanceWalker.TruthVersusEval.falseNegative(truthVC));
            } else {
                evaluations.add(AbstractConcordanceWalker.TruthVersusEval.truePositive(truthVC, evalVC));
            }
        }

        // Actually run the code
        GenomicConcordance genomicConcordance = new GenomicConcordance();
        genomicConcordance.initializeConcordanceStateCounts();
        evaluations.forEach(e -> genomicConcordance.apply(e, null, null));

        // Check block lengths
        for(GenomicConcordanceHistogramEntry entry : expectedBlockLengthHistogram) {
            Assert.assertEquals(genomicConcordance.blockLengthHistogram.get(entry.bin).getTruthValue(), entry.getTruthValue());
            Assert.assertEquals(genomicConcordance.blockLengthHistogram.get(entry.bin).getEvalValue(), entry.getEvalValue());
            // Remove this entry, so we can later check if it's empty and all entries have been visited
            genomicConcordance.blockLengthHistogram.remove(entry.bin);
        }
        Assert.assertEquals(genomicConcordance.blockLengthHistogram.size(), 0);

        // Check confidence scores
        for(GenomicConcordanceHistogramEntry entry : expectedConfidenceHistogram) {
            Assert.assertEquals(genomicConcordance.confidenceHistogram.get(entry.bin).getTruthValue(), entry.getTruthValue());
            Assert.assertEquals(genomicConcordance.confidenceHistogram.get(entry.bin).getEvalValue(), entry.getEvalValue());
            // Remove this entry, so we can later check if it's empty and all entries have been visited
            genomicConcordance.confidenceHistogram.remove(entry.bin);
        }
        Assert.assertEquals(genomicConcordance.confidenceHistogram.size(), 0);

    }
}
