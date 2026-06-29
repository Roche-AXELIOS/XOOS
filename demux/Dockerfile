# syntax=docker/dockerfile:1
ARG CMAKE_BUILD_TYPE=Release

FROM ubuntu:24.04 AS build_env_stage

ENV VENV=/opt/venv

COPY --from=common scripts/install-system-dependencies.sh /common/scripts/install-system-dependencies.sh

RUN python_venv_dir=${VENV} /common/scripts/install-system-dependencies.sh

ENV PATH=${VENV}/bin:${PATH}

ENV CONAN_HOME=/conan/.conan2

COPY --from=vendor . /conan/vendor
COPY --from=common scripts/conan-project-setup.sh /common/scripts/conan-project-setup.sh
COPY conanfile.txt /conan/conanfile.txt

ARG CMAKE_BUILD_TYPE
ARG CONAN_REMOTE_URL
ARG CONAN_BUILD_DEPENDENCIES=missing
ARG CONAN_INCLUDE_VENDOR=True
ARG CONAN_USE_CONANCENTER=True
ARG CONAN_UPLOAD_PACKAGE_PATTERN
RUN --mount=type=secret,id=CONAN_LOGIN_USERNAME \
    --mount=type=secret,id=CONAN_PASSWORD \
    custom_remote_url="${CONAN_REMOTE_URL}" \
    build_type="${CMAKE_BUILD_TYPE}" \
    build_dependencies="${CONAN_BUILD_DEPENDENCIES}" \
    include_vendor="${CONAN_INCLUDE_VENDOR}" \
    use_conancenter="${CONAN_USE_CONANCENTER}" \
    upload_package_pattern="${CONAN_UPLOAD_PACKAGE_PATTERN}" \
    /common/scripts/conan-project-setup.sh

FROM build_env_stage AS build_bin_stage

COPY --from=common . /common
COPY . /demux

ARG BUILD_TESTING=OFF
ARG CODE_COVERAGE_ENABLE=OFF
ARG CMAKE_BUILD_TYPE
ARG VERSION=unknown
ARG LINT_ENABLE=ON

RUN module_dir=/demux \
    build_testing=${BUILD_TESTING} \
    code_coverage_enable=${CODE_COVERAGE_ENABLE} \
    cmake_build_type="${CMAKE_BUILD_TYPE}" \
    version="${VERSION}" \
    /common/scripts/cmake-build.sh

RUN cd /demux \
  && build_testing=${BUILD_TESTING} lint_enable=${LINT_ENABLE} /common/scripts/run-clangd-tidy.sh

FROM ubuntu:24.04 AS adapter_design_bundle_stage

RUN apt-get -qq update && apt-get -qq upgrade && apt-get -qq install --no-install-recommends -y \
      zip \
  && rm -rf /var/lib/apt/lists/*

COPY resources /resources

WORKDIR /resources/adapter-design-bundle
RUN zip -r /adapter-design-bundle.zip *

FROM ubuntu:24.04

RUN apt-get -qq update && apt-get -qq upgrade && rm -rf /var/lib/apt/lists/*

ENV PATH=/apps:${PATH}

COPY --from=adapter_design_bundle_stage /adapter-design-bundle.zip /resources/adapter-design-bundle.zip
COPY --from=build_bin_stage /demux/build/apps/demux /apps/demux
COPY --from=build_bin_stage /demux/build/apps/demux_strand_index /apps/demux_strand_index

WORKDIR /apps

COPY LICENSE.md /LICENSE.md
