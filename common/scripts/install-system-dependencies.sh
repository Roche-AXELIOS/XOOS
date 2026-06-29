#!/bin/bash

set -eu -o pipefail

python_venv_dir="${python_venv_dir:-/opt/venv}"

apt-get update
apt-get upgrade -y
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  cmake \
  curl \
  doxygen \
  g++ \
  gcc \
  gcovr \
  gdb \
  git-lfs \
  git \
  lsb-release \
  make \
  ninja-build \
  python3-pip \
  python3-venv \
  unzip

apt-get purge -y zlib1g-dev
rm -rf /var/lib/apt/lists/*

python3 -m venv "${python_venv_dir}"
source "${python_venv_dir}/bin/activate"
python3 -m pip install -q conan~=2.0 clangd~=20.0 clangd-tidy~=1.0
