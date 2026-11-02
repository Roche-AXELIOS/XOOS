find_package(Eigen3 REQUIRED)
find_package(xoos-libbigwig REQUIRED)
find_package(CLI11 REQUIRED)
find_package(xoos-taskflow REQUIRED)
find_package(ZLIB REQUIRED)
find_package(xoos-cgranges REQUIRED)
find_package(xoos-hclust-cpp REQUIRED)
find_package(xoos-htslib REQUIRED)
find_package(xoos-vincentlaucsb-csv-parser REQUIRED)

add_default_library(copy_number_caller_core
        intervals/augment-intervals.cpp
        coverage/calculate-coverage-taskflow-graph.cpp
        coverage/calculate-coverage.cpp
        coverage/coverage.cpp
        gc-correct/gc-correct-algorithm.cpp
        gc-correct/gc-correct.cpp
        likelihood/flag-copy-number-calls.cpp
        likelihood/likelihood-flags.cpp
        likelihood/likelihood.cpp
        likelihood/baf-tools.cpp
        merge-segments/merge-segments.cpp
        merge-segments/merge-segments-main.cpp
        one-sample-logr/one-sample-logr.cpp
        one-sample-logr/self-normalize.cpp
        predict-sex/predict-sex.cpp
        purity-ploidy-search/purity-ploidy-search.cpp
        cli/add-copy-number-caller-options.cpp
        cli/define-copy-number-caller-options.cpp
        cli/copy-number-caller-cli.cpp
        coverage/bam-to-corrected-coverage.cpp
        intervals/check-chr-y-par.cpp
        coverage/coverage-check.cpp
        io/write-empty-outputs.cpp
        io/write-logr.cpp
        io/augmented-intervals-io.cpp
        germline-normal-wgs/germline-normal-wgs.cpp
        somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.cpp
        tumor-only/tumor-only-wgs.cpp
        io/write-metrics.cpp
        io/column-names.cpp
        seg-to-vcf/seg-to-vcf.cpp
        segmentation/cbs-util.cpp
        segmentation/cbs.cpp
        segmentation/cbs-wrapper.cpp
        segmentation/cbs-for-somatic.cpp
        segmentation/genomic-segments.cpp
        segmentation/interval-trees.cpp
        segmentation/max-t.cpp
        segmentation/parent-specific-binary-segmentation-taskflow-graph.cpp
        segmentation/pscbs.cpp
        segmentation/segment-by-segments-taskflow-graph.cpp
        segmentation/segment-pruning.cpp
        segmentation/segmentation-mode.cpp
        io/write-segments.cpp
        io/write-segment-likelihoods.cpp
        segmentation/read-segments.cpp
        segmentation/undo-segments.cpp
        two-sample-logr/calculate-logrs.cpp
        two-sample-logr/tumor-normal.cpp
        two-sample-logr/two-sample-logr.cpp
        intervals/augmented-intervals.cpp
        io/bigwig.cpp
        io/fai.cpp
        observations/filter-observations.cpp
        intervals/interval-tools.cpp
        util/mapq-util.cpp
        tumor-only/tumor-only-bafs.cpp
        tumor-only/snp-blocklist.cpp
        io/write-baf-values.cpp
        io/write-bigwig.cpp
        observations/observations.cpp
        util/parse-util.cpp
        io/write-igv-xml.cpp
)

target_include_directories(copy_number_caller_core SYSTEM INTERFACE .)
target_include_directories(copy_number_caller_core PRIVATE .)

target_link_libraries(copy_number_caller_core
        PUBLIC
        Eigen3::Eigen
        xoos-libbigwig::xoos-libbigwig
        xoos-cgranges::xoos-cgranges
        xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser
        PRIVATE
        CLI11::CLI11
        xoos-taskflow::xoos-taskflow
        xoos-hclust-cpp::xoos-hclust-cpp
)

target_link_xoos_libraries(copy_number_caller_core
        PRIVATE
        log
        io
        util
        enum
        error
        genome
        gc_correct
        types
        stats
)

target_link_xoos_libraries(copy_number_caller_core
        PUBLIC
        cli
        sex_predict
)


target_compile_definitions(copy_number_caller_core
        PRIVATE
        PROGRAM_NAME="copy_number_caller"
        VERSION="${VERSION}"
)
