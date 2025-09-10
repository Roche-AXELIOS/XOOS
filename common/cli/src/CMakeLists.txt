add_library(xoos_cli
        STATIC
        xoos/cli/bam-option-util.cpp
        xoos/cli/cli.cpp
        xoos/cli/fasta-option-util.cpp
        xoos/cli/file-option-util.cpp
        xoos/cli/thread-count-option-util.cpp
        xoos/cli/validators/bam-index-validator.cpp
        xoos/cli/validators/file-extension-validator.cpp
        xoos/cli/validators/file-permission-validator.cpp
)
target_include_directories(xoos_cli
        PRIVATE
        .
)
target_include_directories(xoos_cli
        SYSTEM
        INTERFACE
        .
)
target_link_xoos_libraries(xoos_cli
        PRIVATE
        log
        util
        types
        enum
)

find_package(CLI11 REQUIRED)
target_link_libraries(xoos_cli PUBLIC CLI11::CLI11)

find_package(fmt REQUIRED)
target_link_libraries(xoos_cli PUBLIC fmt::fmt)
