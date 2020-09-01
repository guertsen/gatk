package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.mutable.MutableLong;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

public class GenomicConcordanceHistogramEntry {
    private static final String BIN_COLUMN_NAME = "bin";
    private static final String TRUTH_COLUMN_NAME = "truth";
    private static final String EVAL_COLUMN_NAME = "eval";

    private static final String[] HISTOGRAM_COLUMN_HEADER =
            {BIN_COLUMN_NAME, TRUTH_COLUMN_NAME, EVAL_COLUMN_NAME};

    final String bin;
    MutableLong truthValue;
    MutableLong evalValue;

    public GenomicConcordanceHistogramEntry(final String bin) {
        this.bin = bin;
        this.truthValue = new MutableLong(0);
        this.evalValue = new MutableLong(0);
    }

    public GenomicConcordanceHistogramEntry(final String bin, final long truthValue, final long evalValue) {
        this.bin = bin;
        this.truthValue = new MutableLong(truthValue);
        this.evalValue = new MutableLong(evalValue);
    }

    public String getBin() {return bin; }
    public long getTruthValue() { return truthValue.longValue(); }
    public long getEvalValue() { return evalValue.longValue(); }

    public void incrementTruthValue() { truthValue.increment(); }
    public void incrementEvalValue() { evalValue.increment(); }

    public static class Writer extends TableWriter<GenomicConcordanceHistogramEntry> {
        private Writer(final Path output) throws IOException {
            super(output, new TableColumnCollection(HISTOGRAM_COLUMN_HEADER));
        }

        @Override
        protected void composeLine(GenomicConcordanceHistogramEntry record, DataLine dataLine) {
            dataLine.set(BIN_COLUMN_NAME, record.getBin())
                    .set(TRUTH_COLUMN_NAME, record.getTruthValue())
                    .set(EVAL_COLUMN_NAME, record.getEvalValue());
        }
    }

    public static GenomicConcordanceHistogramEntry.Writer getWriter(final File outputTable) {
        try {
            return new GenomicConcordanceHistogramEntry.Writer(IOUtils.fileToPath(outputTable));
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable), e);
        }
    }

    public static class Reader extends TableReader<GenomicConcordanceHistogramEntry> {
        public Reader(final Path summary) throws IOException {
            super(summary);
        }

        @Override
        protected GenomicConcordanceHistogramEntry createRecord(final DataLine dataLine) {
            final String bin = dataLine.get(BIN_COLUMN_NAME);
            final long truthValue = Long.parseLong(dataLine.get(TRUTH_COLUMN_NAME));
            final long evalValue = Long.parseLong(dataLine.get(EVAL_COLUMN_NAME));

            return new GenomicConcordanceHistogramEntry(bin, truthValue, evalValue);
        }
    }
}
