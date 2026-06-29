add_library(xoos_cli
        STATIC
        xoos/cli/bam-option-util.cpp
        xoos/cli/cli.cpp
        xoos/cli/fasta-option-util.cpp
        xoos/cli/file-option-util.cpp
        xoos/cli/thread-count-option-util.cpp
        xoos/cli/validators/bam-index-validator.cpp
        xoos/cli/validators/bed-validator.cpp
        xoos/cli/validators/file-extension-validator.cpp
        xoos/cli/validators/file-permission-validator.cpp
        xoos/cli/validators/indexed-bam-validator.cpp
        xoos/cli/validators/indexed-fasta-validator.cpp
        xoos/cli/validators/indexed-vcf-validator.cpp
        xoos/cli/validators/nonempty-file-validator.cpp
        xoos/cli/validators/validator-util.cpp
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

find_package(CLI11 REQUIRED)
target_link_libraries(xoos_cli PUBLIC CLI11::CLI11)

find_package(fmt REQUIRED)
target_link_libraries(xoos_cli PUBLIC fmt::fmt)

find_package(xoos-htslib REQUIRED)
target_link_libraries(xoos_cli PUBLIC xoos-htslib::xoos-htslib)

target_link_xoos_libraries(xoos_cli
        PRIVATE
        log
        util
        types
        enum
        error
        io
)
