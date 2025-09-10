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

    remotes+=(--remote=${custom_remote_name})
  fi

  if [ -n "${upload_package_pattern:-}" ] ; then
    build_types="Release Debug RelWithDebInfo"
  fi

  for build_type in ${build_types} ; do
    deployer=""

    if [ "${deploy_direct_dependencies}" == "True" ] ; then
      deployer="--deployer=direct_deploy --deployer-folder=/conan/build-${build_type}"
    fi

    conan install \
      /conan \
      --profile:all="/conan/vendor/profile-ubuntu-$(lsb_release -rs)-$(uname -m)" \
      --output-folder="/conan/build-${build_type}" \
      --build=${build_dependencies} \
      -c tools.build:skip_test=${skip_test} \
      --settings=build_type=${build_type} \
      ${remotes[@]} \
      ${deployer}
  done
  find /conan -type d -exec chmod a+wx {} +

  if [ -n "${upload_package_pattern:-}" -a -n "${custom_remote_url:-}" ] ; then
    conan upload \
      --confirm \
      --remote ${custom_remote_name} \
      "${upload_package_pattern}"
  else
    echo "upload_package_pattern or custom_remote_url not set, skipping upload"
  fi

  if [ -n "${custom_remote_url:-}" ] ; then
    conan remote logout "${custom_remote_name}"
  fi
}

main
