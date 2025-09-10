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
        core/column-names.cpp
        core/command-line-info.cpp
        core/config.cpp
        core/filtering.cpp
        core/genotype.cpp
        core/model-metadata.cpp
        core/score-calculator.cpp
        core/variant-feature-extraction.cpp
        core/variant-id.cpp
        core/variant-info-serializer.cpp
        core/variant-info.cpp
        core/variant-type.cpp
        filter-variants/filter-variants-cli.cpp
        filter-variants/filter-variants.cpp
        filter-variants/filter-region.cpp
        filter-variants/reconcile-germline.cpp
        train-model/model-trainer.cpp
        train-model/train-model-cli.cpp
        train-model/train-model.cpp
        util/cli-util.cpp
        util/compress/compress.cpp
        util/compress/gzip.cpp
        util/compress/ios-util.cpp
        util/compress/zstandard.cpp
        util/file-util.cpp
        util/lightgbm-util.cpp
        util/locked-tsv-writer.cpp
        util/log-util.cpp
        util/num-utils.cpp
        util/parallel-compute-utils.cpp
        util/region-util.cpp
        util/seq-util.cpp
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
        enum
        error
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
