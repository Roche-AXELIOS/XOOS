add_library(xoos_yc_decode
        STATIC
        xoos/yc-decode/yc-decoder.cpp
        xoos/yc-decode/cigar-utils.cpp
)

target_include_directories(xoos_yc_decode PRIVATE .)
target_include_directories(xoos_yc_decode SYSTEM INTERFACE .)

find_package(xoos-htslib REQUIRED)
target_link_libraries(xoos_yc_decode
        PUBLIC
        xoos-htslib::xoos-htslib
        PRIVATE
)

target_link_xoos_libraries(xoos_yc_decode
        PUBLIC
        error
        log
        PRIVATE
        io
        types
        util
)
