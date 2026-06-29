#!/bin/bash

set -eu -o pipefail

module_dir="${module_dir:-}"
build_testing="${build_testing:-OFF}"
code_coverage_enable="${code_coverage_enable:-OFF}"
cmake_build_type="${cmake_build_type:-Release}"
version="${version:-unknown}"

mkdir -p "${module_dir}/build"
cd "${module_dir}/build"

source "/conan/build-${cmake_build_type}/conanbuild.sh"
cmake .. \
    -G Ninja \
    -DCMAKE_TOOLCHAIN_FILE="/conan/build-${cmake_build_type}/conan_toolchain.cmake" \
    -DBUILD_TESTING="${build_testing}" \
    -DVERSION="${version}" \
    -DCODE_COVERAGE_ENABLE="${code_coverage_enable}" \
    -DCMAKE_BUILD_TYPE="${cmake_build_type}" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DLINT_ENABLE=OFF
cmake --build . --config "${cmake_build_type}" --parallel $(nproc)
source "/conan/build-${cmake_build_type}/deactivate_conanbuild.sh"
