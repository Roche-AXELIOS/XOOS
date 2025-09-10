add_library(xoos_log
        STATIC
        xoos/log/logging.h
        xoos/log/logging.cpp
)
target_include_directories(xoos_log
        PRIVATE
        .
)
target_include_directories(xoos_log
        SYSTEM
        INTERFACE
        .
)

find_package(spdlog REQUIRED)

target_link_libraries(xoos_log
        PUBLIC
        spdlog::spdlog
)
