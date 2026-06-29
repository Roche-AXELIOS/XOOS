add_library(xoos_concurrent
        STATIC
        xoos/concurrent/enumerable-thread-local.h
        xoos/concurrent/enumerable-thread-local.cpp
)
set_target_properties(xoos_concurrent
        PROPERTIES
        LINKER_LANGUAGE CXX
)
target_include_directories(xoos_concurrent
        PRIVATE
        .
)
target_include_directories(xoos_concurrent
        SYSTEM
        INTERFACE
        .
)
