add_library(xoos_vfs
        STATIC
        xoos/vfs/archive.cpp
        xoos/vfs/mem.cpp
        xoos/vfs/std.cpp
        xoos/vfs/vfs.cpp
        xoos/vfs/vfs-iterator.cpp
)
target_include_directories(xoos_vfs
        PRIVATE
        .
)
target_include_directories(xoos_vfs
        SYSTEM
        INTERFACE
        .
)

find_package(xoos-libarchive REQUIRED)

target_link_libraries(xoos_vfs
        PRIVATE
        xoos-libarchive::xoos-libarchive
)
