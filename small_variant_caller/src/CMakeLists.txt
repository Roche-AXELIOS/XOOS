add_library(svc_core
        STATIC
        compute-bam-features/alignment-reader.cpp
        compute-bam-features/compute-bam-features-cli.cpp
        compute-bam-features/compute-bam-features.cpp
        compute-bam-features/compute-bam-region-features.cpp
        compute-bam-features/progress-meter.cpp
        compute-bam-features/read-name-utils.cpp
        compute-bam-features/read-processing.cpp
        compute-bam-features/ref-region.cpp
        compute-bam-features/region.cpp
        compute-vcf-features/compute-vcf-features-cli.cpp
        compute-vcf-features/compute-vcf-features.cpp
        compute-vcf-features/vcf-header-util.cpp
        core/bam-feature-collection.cpp
        core/column-names.cpp
        core/command-line-info.cpp
        core/config.cpp
        core/feature-normalization.cpp
        core/filtering.cpp
        core/genotype.cpp
        core/model-metadata.cpp
        core/score-calculator.cpp
        core/sequencing-protocol.cpp
        core/variant-feature-extraction.cpp
        core/variant-id.cpp
        core/variant-info-serializer.cpp
        core/variant-info.cpp
        core/variant-type.cpp
        core/workflow.cpp
        filter-variants/filter-variants-cli.cpp
        filter-variants/filter-variants.cpp
        filter-variants/filter-region.cpp
        filter-variants/germline-tagging.cpp
        filter-variants/reconcile-germline.cpp
        filter-variants/tumor-only-te-processing.cpp
        filter-variants/tumor-normal-processing.cpp
        train-model/model-trainer.cpp
        train-model/train-model-cli.cpp
        train-model/train-model.cpp
        train-model/training-data-set.cpp
        train-model/truth-label-set.cpp
        util/cli-util.cpp
        util/file-util.cpp
        util/lightgbm-util.cpp
        util/locked-tsv-writer.cpp
        util/log-util.cpp
        util/parallel-compute-utils.cpp
        util/region-util.cpp
        util/seq-util.cpp
        util/vcf-util.cpp
        vcf-to-bed/vcf-to-bed-cli.cpp
        vcf-to-bed/vcf-to-bed.cpp
)
target_include_directories(svc_core
        PRIVATE
        .
)
target_include_directories(svc_core
        SYSTEM
        INTERFACE
        .
)

target_link_xoos_libraries(svc_core
        PRIVATE
        PUBLIC
        log
        util
        io
        cli
        compress
        enum
        error
        sex_predict
        yc_decode
)

find_package(xoos-lightgbm REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Taskflow REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(vincentlaucsb-csv-parser REQUIRED)
find_package(zstd REQUIRED)

target_link_libraries(svc_core
        PUBLIC
        xoos-lightgbm::xoos-lightgbm
        ZLIB::ZLIB
        Taskflow::Taskflow
        nlohmann_json::nlohmann_json
        vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
        zstd::libzstd_static
)
