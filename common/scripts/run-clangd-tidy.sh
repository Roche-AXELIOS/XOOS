#!/bin/bash

set -eu -o pipefail

build_testing="${build_testing:-OFF}"
lint_enable="${lint_enable:-ON}"

# Run clangd-tidy on a standard C++ module directory, clangd-tidy can run
# clang-tidy checks much faster than clang-tidy (up to 10x) at the cost
# of some unsupported checks. This is worthwhile because clang-tidy
# is extremely slow.
#
# This script must be run from the module directory, and after
# a CMake build has been completed with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON.

if [ "${lint_enable}" = "OFF" ] ; then
  exit 0
fi

dirs=()
if [ -d "include" ] ; then
  dirs+=("include")
fi
if [ -d "src" ] ; then
  dirs+=("src")
fi
if [ -d "apps" ] ; then
  dirs+=("apps")
fi
if [ -d "tests" -a "${build_testing}" = "ON" ] ; then
  dirs+=("tests")
fi

find "${dirs[@]}" \( -name '*.h' -or -name '*.cpp' \) \
  | xargs -n 50 clangd-tidy --compile-commands-dir build -j $(nproc)
