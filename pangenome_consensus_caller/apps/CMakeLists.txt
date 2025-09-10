add_default_executable(pangenome_consensus_caller
        pangenome-consensus-caller.cpp
)
target_link_libraries(pangenome_consensus_caller
        PRIVATE
        pangenome_consensus_caller_core
)
