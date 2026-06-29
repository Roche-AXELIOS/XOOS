add_default_executable(demux
        demux.cpp
)
target_link_libraries(demux
        PRIVATE
        demux_core
)
target_link_xoos_libraries(demux
        PRIVATE
        util
        enum
        cli
        concurrent
        vfs
        log
        histogram
        io
)


add_executable(demux_strand_index
        demux-strand-index.cpp
)

target_link_libraries(demux_strand_index
        PRIVATE
        demux_strand_index_core
)

target_link_xoos_libraries(demux_strand_index
        PRIVATE
        cli
        concurrent
        enum
        log
        io
        util
)
