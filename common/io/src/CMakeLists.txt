add_default_library(xoos_io
        xoos/io/alignment-reader.cpp
        xoos/io/bam/aligner.cpp
        xoos/io/bam/read-group.cpp
        xoos/io/bed-region.cpp
        xoos/io/fast-bed-parser.cpp
        xoos/io/library-type.cpp
        xoos/io/fasta-reader.cpp
        xoos/io/htslib-util/cigar-util.cpp
        xoos/io/htslib-util/hts-idx-internal.cpp
        xoos/io/htslib-util/htslib-bcf-ptr.cpp
        xoos/io/htslib-util/htslib-bcf-util.cpp
        xoos/io/htslib-util/htslib-ptr.cpp
        xoos/io/htslib-util/htslib-util.cpp
        xoos/io/htslib-util/kstring.cpp
        xoos/io/htslib-util/read-util.cpp
        xoos/io/malloc-ptr.cpp
        xoos/io/metadata-util.cpp
        xoos/io/vcf/vcf-header.cpp
        xoos/io/vcf/vcf-reader.cpp
        xoos/io/vcf/vcf-record.cpp
        xoos/io/vcf/vcf-writer.cpp
)
target_include_directories(xoos_io
        SYSTEM
        INTERFACE
        .
)
target_include_directories(xoos_io
        PRIVATE
        .
)

find_package(xoos-htslib REQUIRED)
find_package(fmt REQUIRED)
find_package(xoos-vincentlaucsb-csv-parser REQUIRED)
find_package(Microsoft.GSL REQUIRED)

target_link_libraries(xoos_io
        PUBLIC
        xoos-htslib::xoos-htslib
        Microsoft.GSL::GSL
        xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser
        PRIVATE
        fmt::fmt
)
target_link_xoos_libraries(xoos_io
        PUBLIC
        types
        PRIVATE
        util
        log
        error
        enum
        yc_decode
)
