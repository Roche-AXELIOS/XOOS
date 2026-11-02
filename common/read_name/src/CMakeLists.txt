add_default_library(xoos_read_name
        xoos/read-name/consensus-read-name.cpp
        xoos/read-name/read-name-parser.cpp
        xoos/read-name/umi.cpp
)
target_include_directories(xoos_read_name
        PRIVATE
        .
)
target_include_directories(xoos_read_name
        SYSTEM
        INTERFACE
        .
)

find_package(fmt REQUIRED)

target_link_xoos_libraries(xoos_read_name
        PRIVATE
        error
        types
        util
)

target_link_libraries(xoos_read_name
        PUBLIC
        fmt::fmt
)
