find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)

add_default_library(xoos_gc_correct
        xoos/gc_correct/gc-correct.cpp
        xoos/gc_correct/loess.cpp
)
target_include_directories(xoos_gc_correct PRIVATE .)
target_include_directories(xoos_gc_correct SYSTEM INTERFACE .)

target_link_libraries(xoos_gc_correct
        PUBLIC
        Eigen3::Eigen
        PRIVATE
        fmt::fmt
)

target_link_xoos_libraries(xoos_gc_correct
        PUBLIC
        types
        PRIVATE
        log
        util
)
