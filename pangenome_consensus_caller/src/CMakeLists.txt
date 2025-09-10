add_library(pangenome_consensus_caller_core
        STATIC
        pangenome-consensus-caller-core.cpp
)

target_include_directories(pangenome_consensus_caller_core SYSTEM INTERFACE .)
target_include_directories(pangenome_consensus_caller_core PRIVATE .)

find_package(xoos-htslib REQUIRED)

target_link_libraries(pangenome_consensus_caller_core
        PRIVATE
        xoos-htslib::xoos-htslib
)
