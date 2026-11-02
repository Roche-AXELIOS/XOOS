# Contributing

We welcome contributions to XOOS!
Please follow these guidelines to ensure a smooth contribution process.

## Pull Request Process

1. Fork the repository and create a feature branch from `main`
2. Test your changes thoroughly
3. Submit a pull request with a clear description of your changes
4. Address any feedback from code review

## Development Environment

This project uses:

- **Conan** - Package manager
- **CMake** - Build system configuration
- **clang-format** - Code formatting (run before committing)
- **clang-tidy** - Static code analysis and linting
- **Docker** - Development and release target

Please ensure your code passes formatting and linting checks before submitting.

## Building a Module

Each C++ module can be built with Docker or natively with CMake.

### Docker Build

Every module uses a multi-stage Dockerfile that requires two external build
contexts: `common` (shared libraries and CI scripts) and `vendor` (Conan
recipes). The build installs the entire toolchain inside the image, so the
only requirement on your machine is Docker. Run from the module directory:

```bash
cd <module>
docker build \
  --build-arg build_type=Release \
  --build-context common=../common \
  --build-context vendor=../vendor \
  -t <module>:local \
  .
```

Add `--build-arg BUILD_TESTING=ON` to also build the unit tests. Run the
built image:

```bash
docker run --rm --user $(id -u):$(id -g) <module>:local <binary> --help
```

See the [Docker Guide](docs/docker-guide.md) for running images, including
Podman and Singularity/Apptainer.

### Native Build

A native build requires the following toolchain (Ubuntu 24.04 shown):

```bash
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  autoconf automake cmake clang-tidy-18 g++ gcc git git-lfs libtool make \
  ninja-build python3-pip python3-venv

# Ubuntu 24.04 marks the system Python as externally managed (PEP 668),
# so install the Python tooling into a virtual environment.
python3 -m venv /opt/venv
source /opt/venv/bin/activate
pip install 'conan~=2.27' 'clangd~=20.0' 'clangd-tidy~=1.1'
```

| Dependency | Minimum Version |
|------------|-----------------|
| CMake | 3.26+ |
| Conan | 2.x |
| GCC (g++) | 15+ (C++23) |
| Ninja | any recent |
| autoconf / automake / libtool / make | any recent |
| clang-tidy | 18+ (only when linting is enabled) |

One-time Conan setup:

```bash
conan profile detect
git_root="$(git rev-parse --show-toplevel)"
conan remote add vendor "${git_root}/vendor" \
  --force --index=0 --allowed-packages="xoos-*/*"
conan remote add conancenter https://center2.conan.io --force
```

Install dependencies, configure, and build a module:

```bash
cd <module>
conan install \
  --output-folder=build \
  --build=missing \
  --remote=conancenter \
  --remote=vendor \
  --settings=build_type=Release \
  --conf=tools.build:skip_test=true \
  --conf='&:tools.build:skip_test=false' \
  --conf=tools.cmake.cmaketoolchain:generator="Ninja Multi-Config" \
  --profile:all=../vendor/profile-ubuntu-24.04-x86_64 \
  .

cmake --preset conan-default -DBUILD_TESTING=ON
cmake --build --preset conan-release
```

Run the module binary and the tests to verify the build:

```bash
module=<module>
./build/apps/Release/"${module}" --help
./build/tests/Release/"${module}_tests"
```

## Contributor License Agreement (CLA)

All contributors must sign a Contributor License Agreement before their contributions can be accepted.
You will be prompted to sign the CLA when you submit your first pull request.

## Contribution Capacity

Please note that we have limited capacity to review and accept contributions.
We prioritize:

- Documentation improvements
- Bug fixes
- Performance improvements
- Small, focused feature additions

Large feature requests or architectural changes may take longer to review or may not be accepted due to resource constraints.

Thank you for your interest in contributing to XOOS!
