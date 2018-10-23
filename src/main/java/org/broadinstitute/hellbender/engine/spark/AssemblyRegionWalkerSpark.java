package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link AssemblyRegionWalker}. Subclasses should implement {@link #processAssemblyRegions(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class AssemblyRegionWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public final AssemblyRegionReadShardArgumentCollection shardingArgs = new AssemblyRegionReadShardArgumentCollection();

    @ArgumentCollection
    public final AssemblyRegionArgumentCollection assemblyRegionArgs = getAssemblyRegionArgumentCollection();

    /**
     * @return a subclass of {@link AssemblyRegionArgumentCollection} with the default values filled in.
     */
    protected abstract AssemblyRegionArgumentCollection getAssemblyRegionArgumentCollection();

    /**
     * subclasses can override this to control if reads with deletions should be included in isActive pileups
     */
    protected abstract boolean includeReadsWithDeletionsInIsActivePileups();

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    private String referenceFileName;

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    /**
     * @return The evaluator to be used to determine whether each locus is active or not. Must be implemented by tool authors.
     *         The results of this per-locus evaluator are used to determine the bounds of each active and inactive region.
     */
    public abstract AssemblyRegionEvaluator assemblyRegionEvaluator();

    private List<ShardBoundary> intervalShards;

    /**
     * Note that this sets {@code intervalShards} as a side effect, in order to add padding to the intervals.
     */
    @Override
    protected List<SimpleInterval> editIntervals(List<SimpleInterval> rawIntervals) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = rawIntervals == null ? IntervalUtils.getAllIntervalsForReference(sequenceDictionary) : rawIntervals;
        intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardingArgs.readShardSize, shardingArgs.readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        List<SimpleInterval> paddedIntervalsForReads =
                intervals.stream().map(interval -> interval.expandWithinContig(shardingArgs.readShardPadding, sequenceDictionary)).collect(Collectors.toList());
        return paddedIntervalsForReads;
    }

    /**
     * Loads assembly regions and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the assembly regions.
     *
     * @return all assembly regions as a {@link JavaRDD}, bounded by intervals if specified.
     */
    protected JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegions(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        JavaRDD<Shard<GATKRead>> shardedReads = SparkSharder.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, shardingArgs.readShardSize, shuffle);
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getAssemblyRegionsFunction(referenceFileName, bFeatureManager, getHeaderForReads(),
                assemblyRegionEvaluator(), assemblyRegionArgs, includeReadsWithDeletionsInIsActivePileups()));
    }

    private static FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext> getAssemblyRegionsFunction(
            final String referenceFileName,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMFileHeader header,
            final AssemblyRegionEvaluator evaluator,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean includeReadsWithDeletionsInIsActivePileups) {
        return (FlatMapFunction<Shard<GATKRead>, AssemblyRegionWalkerContext>) shardedRead -> {
            ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            final FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(
                    new ShardToMultiIntervalShardAdapter<>(shardedRead),
                    header, reference, features, evaluator,
                    assemblyRegionArgs.minAssemblyRegionSize, assemblyRegionArgs.maxAssemblyRegionSize,
                    assemblyRegionArgs.assemblyRegionPadding, assemblyRegionArgs.activeProbThreshold,
                    assemblyRegionArgs.maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups);
            final Iterable<AssemblyRegion> assemblyRegions = () -> assemblyRegionIter;
            return Utils.stream(assemblyRegions).map(assemblyRegion ->
                    new AssemblyRegionWalkerContext(assemblyRegion,
                        new ReferenceContext(reference, assemblyRegion.getExtendedSpan()),
                        new FeatureContext(features, assemblyRegion.getExtendedSpan()))).iterator();
        };
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferenceFileName());
        processAssemblyRegions(getAssemblyRegions(ctx), ctx);
    }

    /**
     * Process the assembly regions and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link AssemblyRegionWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processAssemblyRegions(JavaRDD<AssemblyRegionWalkerContext> rdd, JavaSparkContext ctx);

}
