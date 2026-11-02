add_default_library(xoos_util
        xoos/util/adopt.h
        xoos/util/char-functions.h
        xoos/util/container-functions.h
        xoos/util/file-functions.cpp
        xoos/util/format-functions.cpp
        xoos/util/hash.h
        xoos/util/ios-util.cpp
        xoos/util/math.cpp
        xoos/util/missing-value.cpp
        xoos/util/parse-float.cpp
        xoos/util/parse-int.cpp
        xoos/util/sequence-functions.cpp
        xoos/util/signal-handler.cpp
        xoos/util/simple-cache.cpp
        xoos/util/string-functions.cpp
        xoos/util/tmp-dir.cpp
        xoos/util/tmp-file.cpp
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

target_link_xoos_libraries(xoos_util
        PRIVATE
        error
        types
)

target_link_libraries(xoos_util
        PUBLIC
        fmt::fmt
)
