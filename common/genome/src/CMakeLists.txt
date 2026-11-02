add_default_library(xoos_genome
        xoos/genome/region.cpp
)

target_include_directories(xoos_genome PRIVATE .)
target_include_directories(xoos_genome SYSTEM INTERFACE .)

find_package(fmt REQUIRED)

target_link_libraries(xoos_genome
        PUBLIC
        fmt::fmt
)

target_link_xoos_libraries(xoos_genome
        PUBLIC
        types
        PRIVATE
        error
        util
)
