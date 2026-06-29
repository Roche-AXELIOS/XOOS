#!/bin/bash

set -eu -o pipefail

build_types="${build_type:-Release Debug}"
build_dependencies="${build_dependencies:-missing}"
custom_remote_name="${custom_remote_name:-custom}"
skip_test="${skip_test:-True}"
include_vendor="${include_vendor:-True}"
deploy_direct_dependencies="${deploy_direct_dependencies:-False}"
use_conancenter="${use_conancenter:-True}"

function main() {
  remotes=()
  if [ "${use_conancenter}" == "True" ] ; then
    remotes+=(--remote=conancenter)
  fi

  if [ "${include_vendor}" == "True" ] ; then
    conan remote add \
      vendor \
      /conan/vendor \
      --allowed-packages="xoos-*/*" \
      --index=0
    remotes+=(--remote=vendor)
  fi

  if [ -n "${custom_remote_url:-}" ] ; then
    conan remote add \
      "${custom_remote_name}" \
      "${custom_remote_url}" \
      --index=1
    conan remote login \
      "${custom_remote_name}" \
      "$(cat /run/secrets/CONAN_LOGIN_USERNAME)" \
      -p "$(cat /run/secrets/CONAN_PASSWORD)"
    remotes+=(--remote="${custom_remote_name}")
  fi

  if [ -n "${upload_package_pattern:-}" ] ; then
    build_types="Release Debug RelWithDebInfo"
  fi

  local conan_install_cmd=(
    conan install
    /conan
    "--profile:all=/conan/vendor/profile-ubuntu-$(lsb_release -rs)-$(uname -m)"
    "--build=${build_dependencies}"
    "--conf=tools.build:skip_test=${skip_test}"
    "--conf=tools.cmake.cmaketoolchain:generator=Ninja"
  )
  conan_install_cmd+=("${remotes[@]}")

  for build_type in ${build_types} ; do
    local current_cmd=("${conan_install_cmd[@]}")
    current_cmd+=(
      "--output-folder=/conan/build-${build_type}"
      "--settings=build_type=${build_type}"
    )

    if [ "${deploy_direct_dependencies}" == "True" ] ; then
      current_cmd+=(
        "--deployer=direct_deploy"
        "--deployer-folder=/conan/build-${build_type}"
      )
    fi

    echo "--- Installing for build_type: ${build_type} ---"
    "${current_cmd[@]}"
  done

  find /conan -type d -exec chmod a+wx {} +

  if [ -n "${upload_package_pattern:-}" -a -n "${custom_remote_url:-}" ] ; then
    conan upload \
      --confirm \
      --remote "${custom_remote_name}" \
      "${upload_package_pattern}"
  else
    echo "upload_package_pattern or custom_remote_url not set, skipping upload"
  fi

  if [ -n "${custom_remote_url:-}" ] ; then
    conan remote logout "${custom_remote_name}"
  fi
}

main
