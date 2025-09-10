add_library(xoos_util
        STATIC
        xoos/util/container-functions.h
        xoos/util/distance-functions.cpp
        xoos/util/file-functions.cpp
        xoos/util/hash.h
        xoos/util/math.cpp
        xoos/util/parser.cpp
        xoos/util/sequence-functions.cpp
        xoos/util/signal-handler.cpp
        xoos/util/simple-cache.cpp
        xoos/util/string-functions.cpp
        xoos/util/tmp-dir.cpp
        xoos/util/uuid.cpp
)
target_include_directories(xoos_util
        PRIVATE
        .
)
target_include_directories(xoos_util
        SYSTEM
        INTERFACE
        .
)

find_package(fmt REQUIRED)
find_package(stduuid REQUIRED)

target_link_xoos_libraries(xoos_util
        PRIVATE
        error
        types
)

target_link_libraries(xoos_util
        PUBLIC
        fmt::fmt
        stduuid::stduuid
)
