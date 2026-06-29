add_library(xoos_compress
        STATIC
        xoos/compress/compress.cpp
        xoos/compress/gzip.cpp
        xoos/compress/zstandard.cpp
)
target_include_directories(xoos_compress
        PRIVATE
        .
)
target_include_directories(xoos_compress
        SYSTEM
        INTERFACE
        .
)

find_package(fmt REQUIRED)
find_package(ZLIB REQUIRED)
find_package(zstd REQUIRED)

target_link_xoos_libraries(xoos_compress
        PRIVATE
        error
        log
        types
        util
)

target_link_libraries(xoos_compress
        PUBLIC
        fmt::fmt
        ZLIB::ZLIB
        zstd::libzstd_static
)
