#!/bin/bash

set -eu -o pipefail

cat > /etc/apt/apt.conf.d/99-defaults <<'EOF'
# 1. MAKE LOGS LESS NOISY
# Equivalent to -q=2. Implies -y.
# 0=quiet, 1=errors only, 2=quiet (no progress bars)
quiet "2";

# Disable the "progress" bar entirely (reduces noise in CI/logs)
Dpkg::Use-Pty "0";

# 2. DON'T INSTALL RECOMMENDS/SUGGESTS
# Only install exactly what is requested and its hard dependencies
APT::Install-Recommends "0";
APT::Install-Suggests "0";

# 3. MAKE IT NON-INTERACTIVE
# Automatically answer "yes" to prompts
APT::Get::Assume-Yes "true";

# FIX DPKG CONFIGURATION PROMPTS
# If a config file changes, keep the old version without asking
# (or use "force-confnew" to always overwrite)
Dpkg::Options {
   "--force-confdef";
   "--force-confold";
};

APT::Keep-Downloaded-Packages "false";
EOF

export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get upgrade
apt-get install \
  autoconf \
  automake \
  cmake \
  curl \
  doxygen \
  g++ \
  gcc \
  gcovr \
  gdb \
  git \
  git-lfs \
  lsb-release \
  make \
  ninja-build \
  python3-pip \
  python3-venv \
  unzip

python_venv_dir="${python_venv_dir:-/opt/venv}"

python3 -m venv "${python_venv_dir}"
source "${python_venv_dir}/bin/activate"

python3 -m pip install -q --upgrade pip
python3 -m pip install -q conan~=2.27 clangd~=20.0 clangd-tidy~=1.1

# sccache (opt-in, set sccache_enable=true to install)
sccache_enable="${sccache_enable:-false}"
if [ "${sccache_enable}" = "true" ] ; then
  mkdir -p /ci/sccache
  # -f: fail on HTTP errors instead of saving the error body; --retry: tolerate transient failures
  curl -fsSL --retry 5 --retry-all-errors --retry-delay 2 -o sccache-v0.14.0-x86_64-unknown-linux-musl.tar.gz https://github.com/mozilla/sccache/releases/download/v0.14.0/sccache-v0.14.0-x86_64-unknown-linux-musl.tar.gz
  echo '8424b38cda4ecce616a1557d81328f3d7c96503a171eab79942fad618b42af44 sccache-v0.14.0-x86_64-unknown-linux-musl.tar.gz' | sha256sum -c
  tar xz -C /ci/sccache --strip-components=1 -f sccache-v0.14.0-x86_64-unknown-linux-musl.tar.gz
  rm sccache-v0.14.0-x86_64-unknown-linux-musl.tar.gz
fi

# Rust toolchain (opt-in, set rust_enable=true for modules that need it)
rust_enable="${rust_enable:-false}"
if [ "${rust_enable}" = "true" ] ; then
  rust_toolchain="${rust_toolchain:-1.91}"
  rust_nightly="${rust_nightly:-nightly-2026-03-08}"
  cbindgen_version="${cbindgen_version:-0.29.2}"

  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain "${rust_toolchain}" --profile minimal
  export PATH="${HOME}/.cargo/bin:${PATH}"

  rustup toolchain install "${rust_nightly}"
  cargo install cbindgen --version "${cbindgen_version}" --locked

  rustc --version && cargo --version && cbindgen --version
fi
