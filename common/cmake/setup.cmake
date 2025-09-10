macro(setup_default_build_type)
    if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
          set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    endif()
endmacro()

macro(setup_default_cxx warning_as_error_disable)
    # clear any currently set compiler flags a future improvement
    # would be to use target specific compiler flags
    unset(CMAKE_CXX_FLAGS)

    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")

    # this is added as an escape hatch, since newer compiler versions can introduce new warnings
    # which could lead to the code no longer building
    # there is native support in this as of cmake 3.24, but we are using 3.22
    # https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_WARNING_AS_ERROR.html
    if (NOT ${warning_as_error_disable})
        option(COMPILE_NO_WARNING_AS_ERROR "Disable compiler warnings as error (-Werror)" OFF)
        if (NOT COMPILE_NO_WARNING_AS_ERROR)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
        endif ()
    endif ()
endmacro()

function(setup_target_link_libraries_util)
    include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/target_link_libraries_util.cmake)
endfunction()

macro(setup_lint)
    unset(CMAKE_CXX_CLANG_TIDY)
    option(LINT_ENABLE "Enable lint" ON)
    if (LINT_ENABLE)
        # to support newer clangd-tidy and older clang-tidy projects
        find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-18" REQUIRED)
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE} -header-filter=. -warnings-as-errors=*)
    endif ()
endmacro()

macro(setup_coverage)
    option(CODE_COVERAGE_ENABLE "Enable code coverage generation" OFF)
    if (CODE_COVERAGE_ENABLE)
        set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
        # turn on atomic profile count updates, this fixes issues with corrupted
        # counts when multi-threading is used. also disable elide constructors and default inline to provide
        # more accurate line uncovered counting
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-update=atomic -fno-elide-constructors -fno-default-inline")
    endif ()
endmacro()

macro(setup_doc name sources)
    option(BUILD_DOC "Enable documentation generation" OFF)
    if (BUILD_DOC)
        find_package(Doxygen)
        if (DOXYGEN_FOUND)
            doxygen_add_docs(${name}_doc ${sources} ALL)
        endif ()
    endif ()
endmacro()

macro(setup_tests sources)
    option(BUILD_TESTING "Enable tests build" ON)
    if (BUILD_TESTING)
        add_subdirectory(${sources})
    endif ()
endmacro()

macro(setup_address_sanitizer)
    option(ADDRESS_SANITIZER_ENABLE "Enable address sanitizer" OFF)
    if (ADDRESS_SANITIZER_ENABLE)
        set(STATIC_LINK_DISABLE ON)
        add_compile_options(-fsanitize=address -fsanitize-recover=address)
        add_link_options(-fsanitize=address -fsanitize-recover=address)
    endif ()
endmacro()

macro(setup_static_link static_link_disable_default)
    option(STATIC_LINK_DISABLE "Disable static link of executable" ${static_link_disable_default})
    function(add_default_executable target sources)
        if (NOT DEFINED VERSION)
            set(VERSION unknown)
        endif ()
        add_executable(${target}
            ${sources}
        )
        target_compile_definitions(${target}
            PRIVATE
            PROGRAM_NAME="${target}"
            VERSION="${VERSION}"
        )
        if (NOT STATIC_LINK_DISABLE)
            target_link_libraries(${target}
                PRIVATE
                -static-libstdc++
                -static-libgcc
            )
            if (NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                target_link_libraries(${target}
                    PRIVATE
                    -static
                )
            endif ()
        endif ()
    endfunction()
endmacro()

# parameters:
#   LINT_DISABLE - turn off clang-tidy
#   WARNING_AS_ERROR - enable warnings, but disable failing with error (-Werror)
#   ADD_SRC_DISABLE - module doesn't have a src directory to add
macro(setup_default_app_project name)
    cmake_parse_arguments(SETUP_DEFAULT_ARG "LINT_DISABLE;WARNING_AS_ERROR_DISABLE;STATIC_LINK_DISABLE" "" "" ${ARGN})
    if (DEFINED SETUP_DEFAULT_ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unrecognized arguments: ${SETUP_DEFAULT_ARG_UNPARSED_ARGUMENTS}")
    endif ()

    setup_default_build_type()
    setup_default_cxx(${SETUP_DEFAULT_ARG_WARNING_AS_ERROR_DISABLE})
    setup_target_link_libraries_util()
    setup_address_sanitizer()
    setup_static_link(${SETUP_DEFAULT_ARG_STATIC_LINK_DISABLE})
    if (NOT ${SETUP_DEFAULT_ARG_LINT_DISABLE})
        setup_lint()
    endif ()
    setup_coverage()
    setup_doc(${name} apps include src)

    add_subdirectory(src)
    add_subdirectory(apps)

    setup_tests(tests)
endmacro()

# parameters:
#   LINT_DISABLE - turn off clang-tidy
#   WARNING_AS_ERROR - enable warnings, but disable failing with error (-Werror)
#   ADD_SRC_DISABLE - module doesn't have a src directory to add
function(setup_default_libxoos_project name)
    cmake_parse_arguments(SETUP_DEFAULT_ARG "LINT_DISABLE;WARNING_AS_ERROR_DISABLE;ADD_SRC_DISABLE" "" "" ${ARGN})
    if (DEFINED SETUP_DEFAULT_ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unrecognized arguments: ${SETUP_DEFAULT_ARG_UNPARSED_ARGUMENTS}")
    endif ()

    setup_default_build_type()
    setup_default_cxx(${SETUP_DEFAULT_ARG_WARNING_AS_ERROR_DISABLE})
    setup_target_link_libraries_util()
    if (NOT ${SETUP_DEFAULT_ARG_LINT_DISABLE})
        setup_lint()
    else ()
        unset(CMAKE_CXX_CLANG_TIDY)
    endif ()
    setup_coverage()
    setup_doc(${name} apps include src)

    if (NOT ${SETUP_DEFAULT_ARG_ADD_SRC_DISABLE})
        add_subdirectory(src)
    endif ()

    if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
        setup_tests(tests)
    endif ()
endfunction()
