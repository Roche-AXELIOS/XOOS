add_library(xoos_io
        STATIC
        xoos/io/bed-region.cpp
        xoos/io/fasta-reader.cpp
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
find_package(vincentlaucsb-csv-parser REQUIRED)
find_package(Microsoft.GSL REQUIRED)

target_link_libraries(xoos_io
        PUBLIC
        xoos-htslib::xoos-htslib
        Microsoft.GSL::GSL
        PRIVATE
        vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
        fmt::fmt
)
target_link_xoos_libraries(xoos_io
        PUBLIC
        types
        PRIVATE
        util
        log
        error
        yc_decode
)
