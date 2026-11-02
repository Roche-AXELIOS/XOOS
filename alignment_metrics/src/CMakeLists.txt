add_default_library(alignment_metrics_core
        core/alignment.cpp
        core/hp-error-for-read.cpp
        core/interval.cpp
        core/region-labels-set.cpp
        core/pileup.cpp
        core/region-lookup.cpp
        core/supplementary-alignment.cpp
        core/substitution-lookup.cpp
        core/super-region.cpp
        io/bed-reader.cpp
        io/output-path-validation.cpp
        metadata/alignment-metadata.cpp
        metadata/dataset-metadata.cpp
        metrics/accuracy-metrics/accuracy-metrics.cpp
        metrics/accuracy-metrics/base-level-accuracy-summary.cpp
        metrics/accuracy-metrics/errors-by-cluster-size.cpp
        metrics/accuracy-metrics/errors-by-read-type.cpp
        metrics/accuracy-metrics/errors-by-substitution-type.cpp
        metrics/accuracy-metrics/hp-stat.cpp
        metrics/accuracy-metrics/qscore-stats.cpp
        metrics/coverage-metrics/coverage-distribution-summary.cpp
        metrics/coverage-metrics/coverage-histograms.cpp
        metrics/coverage-metrics/coverage-metrics.cpp
        metrics/coverage-metrics/coverage-uniformity-metrics.cpp
        metrics/coverage-metrics/gc-bias-metrics.cpp
        metrics/read-metrics/read-counts-by-cluster-size.cpp
        metrics/read-metrics/read-length-histograms.cpp
        metrics/read-metrics/read-metrics.cpp
        metrics/read-metrics/read-metrics-summary.cpp
        metrics/read-metrics/read-metrics-te.cpp
        metrics/metrics.cpp
        metrics/stratified-metrics.cpp
        util/alignment-trimming-util.cpp
        util/homopolymer-util.cpp
        util/region-util.cpp
        alignment-metrics-cli.cpp
        alignment-metrics-options.cpp
        alignment-metrics.cpp
        metrics-worker.cpp
        run-alignment-metrics.cpp
        util/labeled-region-util.cpp
)
target_include_directories(alignment_metrics_core
        PRIVATE
        .
)
target_include_directories(alignment_metrics_core
        SYSTEM
        INTERFACE
        .
)
target_link_xoos_libraries(alignment_metrics_core
        PUBLIC
        histogram
        io
        types
        yc_decode
        PRIVATE
        cli
        enum
        error
        log
        read_name
        util
)

find_package(xoos-taskflow REQUIRED)
find_package(xoos-vincentlaucsb-csv-parser REQUIRED)
find_package(Microsoft.GSL REQUIRED)
find_package(unordered_dense REQUIRED)

target_link_libraries(alignment_metrics_core
        PUBLIC
        xoos-taskflow::xoos-taskflow
        PRIVATE
        xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser
        Microsoft.GSL::GSL
        unordered_dense::unordered_dense
)

# merge_and_label_bed: a standalone tool that builds labeled BED files consumed by
# `alignment_metrics --stratify-by-label`. Reuses the region-label merge logic in alignment_metrics_core.
add_default_library(merge_and_label_bed_core
        merge-and-label-bed/merge-and-label-bed.cpp
        merge-and-label-bed/merge-and-label-bed-cli.cpp
)
target_include_directories(merge_and_label_bed_core PRIVATE .)
target_include_directories(merge_and_label_bed_core SYSTEM INTERFACE .)
target_link_libraries(merge_and_label_bed_core PUBLIC alignment_metrics_core)
target_link_xoos_libraries(merge_and_label_bed_core
        PRIVATE
        cli
        error
        io
        log
        util
)
