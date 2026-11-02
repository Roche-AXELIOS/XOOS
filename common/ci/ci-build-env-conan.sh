#!/bin/bash

set -eu -o pipefail

# CMake 3.31+ (Ubuntu 26.04) dropped compatibility with cmake_minimum_required
# versions below 3.5. Some ConanCenter packages (e.g. edlib/1.2.7) still use
# older versions. This env var restores compatibility globally.
export CMAKE_POLICY_VERSION_MINIMUM=3.5

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

build_dependencies="${build_dependencies:-missing}"

custom_remote_name="${custom_remote_name:-custom}"
custom_remote_url="${custom_remote_url:-}"
custom_remote_username="${custom_remote_username:-}"
custom_remote_password="${custom_remote_password:-}"
custom_remote_upload_pattern="${custom_remote_upload_pattern:-}"
conancenter_enabled="${conancenter_enabled:-true}"

custom_remote_enabled="False"
if [ -n "${custom_remote_url}" -a -n "${custom_remote_username}" -a -n "${custom_remote_password}" ] ; then
  custom_remote_enabled="True"
fi

remotes=()
if [ "${conancenter_enabled}" == "true" ] ; then
    remotes=(--remote=conancenter)
fi

function cleanup() {
  conan remote remove vendor

  if [ "${custom_remote_enabled}" == "True" ] ; then
    conan remote logout "${custom_remote_name}"
    conan remote remove "${custom_remote_name}"
  fi
}
trap cleanup EXIT

vendor="${script_dir}/../vendor"
conan remote add \
  vendor \
  "${vendor}" \
  --force \
  --index=0 \
  --allowed-packages="xoos-*/*"
remotes+=(--remote=vendor)

if [ "${custom_remote_enabled}" == "True" ] ; then
  conan remote add \
    "${custom_remote_name}" \
    "${custom_remote_url}" \
    --index=1
  conan remote login \
    "${custom_remote_name}" \
    "${custom_remote_username}" \
    -p "${custom_remote_password}"
  remotes+=(--remote="${custom_remote_name}")
fi

cmake_extra_vars="{}"

for build_type in Release Debug RelWithDebInfo ; do
  profile="${vendor}/profile-ubuntu-$(lsb_release -rs)-$(uname -m)"
  conan install \
    --output-folder=/ci/build \
    --build="${build_dependencies}" \
    "${remotes[@]}" \
    --settings=build_type=$build_type \
    --conf=tools.build:skip_test=true \
    --conf=\&:tools.build:skip_test=false \
    --conf=tools.cmake.cmaketoolchain:generator="Ninja Multi-Config" \
    --conf=tools.cmake.cmaketoolchain:extra_variables="${cmake_extra_vars}" \
    --profile:all="${profile}" \
    .
done

if [ -n "${custom_remote_upload_pattern}" -a "${custom_remote_enabled}" == "True" ] ; then
  conan upload \
    --confirm \
    --remote "${custom_remote_name}" \
    "${custom_remote_upload_pattern}"
fi

# Make /ci/build writable and the Conan package cache readable for the
# non-root docker run user. /ci/build is where cmake writes build artifacts.
# /root/.conan2 contains headers and libraries referenced by absolute path
# in the Conan-generated cmake config files; /root itself must be traversable
# to reach .conan2.
chmod -R a+rwX /ci/build
chmod a+rx /root
chmod -R a+rX /root/.conan2
