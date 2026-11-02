include(CheckCXXCompilerFlag)

# demux_shared_core - shared between demux_core and demux_strand_index_core
add_default_library(demux_shared_core
        bloom-filter/interleaved-bloom-filter.cpp
        sequence/strand/strand-kmer-table.cpp
        sequence/kmer/kmerizer.cpp
)

target_include_directories(demux_shared_core PRIVATE .)
target_include_directories(demux_shared_core SYSTEM INTERFACE .)

target_link_xoos_libraries(demux_shared_core
        PUBLIC
        types
        util
)

find_package(gtl REQUIRED)
target_link_libraries(demux_shared_core PUBLIC gtl::gtl)

find_package(xoos-taskflow REQUIRED)
target_link_libraries(demux_shared_core PUBLIC xoos-taskflow::xoos-taskflow)

find_package(xxHash REQUIRED)
target_link_libraries(demux_shared_core PUBLIC xxHash::xxhash)

# demux_core
add_default_library(demux_core
        adapter-design/adapter-design-bundle.cpp
        adapter-design/adapter-design.cpp
        adapter-design/generate-edited-sequence.cpp
        adapter-dsl/model/adapter-model.cpp
        adapter-dsl/model/adapter-entry-json.cpp
        adapter-dsl/model/compiled-json.cpp
        adapter-dsl/dsl/dsl-compiler.cpp
        adapter-dsl/parser/read-structure-parser.cpp
        adapter-dsl/parser/builtin-parsers.cpp
        adapter-dsl/authored/authored-json.cpp
        adapter-dsl/resolve/resolved-adapter.cpp
        adapter-dsl/assemble/adapter-design-assembler.cpp
        adapters/duplex/demux-and-trim-duplex.cpp
        adapters/duplex/duplex-match.cpp
        adapters/duplex/hairpin.cpp
        adapters/duplex/hairpin-finder.cpp
        adapters/duplex/lut-bundle-duplex.cpp
        adapters/duplex/trim-duplex.cpp
        adapters/duplex/umi/demux-and-trim-duplex-umi.cpp
        adapters/duplex/umi/lut-bundle-duplex-umi.cpp
        adapters/duplex/umi/trim-duplex-umi.cpp
        adapters/duplex/stem/demux-and-trim-duplex-stem.cpp
        adapters/duplex/stem/lut-bundle-duplex-stem.cpp
        adapters/duplex/stem/trim-duplex-stem.cpp
        adapters/simplex/demux-and-trim-simplex.cpp
        adapters/simplex/dual/demux-and-trim-simplex-dual-paired.cpp
        adapters/simplex/paired/demux-and-trim-simplex-paired.cpp
        adapters/simplex-10x/lut-bundle-simplex-10x.cpp
        adapters/simplex-10x/demux-and-trim-simplex-10x.cpp
        adapters/simplex-10x/sid-fast-match.cpp
        adapters/ysu/demux-and-trim-ysu.cpp
        adapters/ysu/lut-bundle-ysu.cpp
        adapters/ysu/trim3p-ysu.cpp
        adapters/ysu/trim5p-ysu.cpp
        adapters/ys/demux-and-trim-ys.cpp
        adapters/ys/lut-bundle-ys.cpp
        adapters/ys/trim3p-ys.cpp
        adapters/ys/trim5p-ys.cpp
        cli/demux.cpp
        core/demux-and-trim-pipeline.cpp
        io/fastq/kseqpp-sequence-reader.cpp
        io/rdb/rdb-2bit-utils.cpp
        io/rdb/rdb-sequence-reader.cpp
        io/sample-sheet/sample-sheet.cpp
        io/sequence-reader.cpp
        metrics/duplex-metrics.cpp
        metrics/metrics.cpp
        metrics/simplex-metrics.cpp
        sequence/loci-range.cpp
        sequence/matcher/bitap.cpp
        sequence/matcher/dfa-classifier.cpp
        sequence/matcher/match-info.cpp
        sequence/matcher/seq-lut.cpp
        sequence/matcher/seq-matcher.cpp
        sequence/strand/strand-classifier.cpp
        simd/simd-functions.cpp
        task/alignment.cpp
        task/demux.cpp
        task/flow-context.cpp
        task/flow-manager.cpp
        task/formatter.cpp
        task/reader.cpp
        task/sink-worker.cpp
        task/sink.cpp
        task/task.cpp
        utility/alignment-util.cpp
        utility/fmt-number.cpp
        utility/string-util.cpp
)

target_include_directories(demux_core PRIVATE .)
target_include_directories(demux_core SYSTEM INTERFACE .)

target_link_libraries(demux_core PUBLIC demux_shared_core)

target_link_xoos_libraries(demux_core
        PUBLIC
        error
        PRIVATE
        cli
        concurrent
        enum
        histogram
        io
        log
        vfs
)

find_package(xoos-highway REQUIRED)
target_link_libraries(demux_core PUBLIC xoos-highway::xoos-highway)

find_package(xoos-kseqpp REQUIRED)
target_link_libraries(demux_core PUBLIC xoos-kseqpp::xoos-kseqpp)

find_package(nlohmann_json REQUIRED)
target_link_libraries(demux_core PRIVATE nlohmann_json::nlohmann_json)

find_package(xoos-vincentlaucsb-csv-parser REQUIRED)
target_link_libraries(demux_core PUBLIC xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser)

find_package(ZLIB REQUIRED)
target_link_libraries(demux_core PUBLIC ZLIB::ZLIB)

find_package(zstd REQUIRED)
target_link_libraries(demux_core PUBLIC zstd::libzstd_static)


find_package(Boost REQUIRED)
target_link_libraries(demux_core PRIVATE Boost::headers)

# required for Highway SIMD functions
add_compile_definitions(SIMD_FUNCTIONS_FILE="${CMAKE_CURRENT_SOURCE_DIR}/simd/simd-functions.cpp")

target_compile_options(demux_core PRIVATE
        $<$<CONFIG:Release>:-O3 -flto>
        $<$<CONFIG:ReleaseWithDebInfo>:-g -O3 -flto>
)

target_link_options(demux_core PRIVATE
    $<$<CONFIG:Release>:-flto>
    $<$<CONFIG:ReleaseWithDebInfo>:-flto>
)

check_cxx_compiler_flag("-mpopcnt" POPCNT_SUPPORTED)
if (POPCNT_SUPPORTED)
    target_compile_options(demux_core PRIVATE -mpopcnt)
else ()
    message(WARNING "Compiler does not support -mpopcnt")
endif ()


# demux_strand_index_core START
add_default_library(
        demux_strand_index_core
        demux-strand-index/demux-strand-build-pipeline.cpp
        demux-strand-index/strand-build.cpp
)

target_include_directories(demux_strand_index_core PRIVATE .)
target_include_directories(demux_strand_index_core SYSTEM INTERFACE .)

target_link_libraries(demux_strand_index_core PUBLIC demux_shared_core)

target_link_xoos_libraries(demux_strand_index_core
        PRIVATE
        cli
        concurrent
        enum
        error
        io
        log
)
