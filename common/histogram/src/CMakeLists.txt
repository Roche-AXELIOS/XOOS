add_default_library(xoos_histogram
        xoos/histogram/histogram-summary.cpp
)

target_include_directories(xoos_histogram PRIVATE .)
target_include_directories(xoos_histogram SYSTEM INTERFACE .)

find_package(xoos-htslib REQUIRED)
find_package(xoos-vincentlaucsb-csv-parser REQUIRED)
find_package(fmt REQUIRED)
find_package(Microsoft.GSL REQUIRED)
target_link_libraries(xoos_histogram
        PUBLIC
        xoos-htslib::xoos-htslib
        xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser
        fmt::fmt
        PRIVATE
        Microsoft.GSL::GSL
)

target_link_xoos_libraries(xoos_histogram
        PUBLIC
        io
        types
        PRIVATE
        util
)
