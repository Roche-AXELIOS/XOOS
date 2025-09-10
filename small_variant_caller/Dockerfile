# syntax=docker/dockerfile:1
ARG CMAKE_BUILD_TYPE=Release

FROM ubuntu:24.04 AS build_env_stage

ENV VENV=/opt/venv

RUN apt-get update && apt-get upgrade -y \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
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
      unzip \
  && apt-get purge -y zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && python3 -m venv "${VENV}" \
  && . "${VENV}/bin/activate" \
  && python3 -m pip install -q conan~=2.0 clangd~=20.0 clangd-tidy~=1.0

ENV PATH=${VENV}/bin:${PATH}

ENV CONAN_HOME=/conan/.conan2

COPY --from=vendor . /conan/vendor

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
    /conan/vendor/conan-project-setup.sh

FROM build_env_stage AS build_bin_stage

COPY --from=common / /common
COPY . /small_variant_caller

ARG BUILD_TESTING=OFF
ARG BUILD_DOC=OFF
ARG CODE_COVERAGE_ENABLE=OFF
ARG CMAKE_BUILD_TYPE
ARG VERSION=unknown
ARG LINT_ENABLE=ON

RUN mkdir -p /small_variant_caller/build \
    && cd /small_variant_caller/build \
    && . /conan/build-${CMAKE_BUILD_TYPE}/conanbuild.sh \
    && cmake .. \
        -G Ninja \
        -DCMAKE_TOOLCHAIN_FILE=/conan/build-${CMAKE_BUILD_TYPE}/conan_toolchain.cmake \
        -DBUILD_TESTING=${BUILD_TESTING} \
        -DBUILD_DOC=${BUILD_DOC} \
        -DVERSION=${VERSION} \
        -DCODE_COVERAGE_ENABLE=${CODE_COVERAGE_ENABLE} \
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DLINT_ENABLE=OFF \
    && cmake --build . --config ${CMAKE_BUILD_TYPE} --parallel $(nproc) \
    && . /conan/build-${CMAKE_BUILD_TYPE}/deactivate_conanbuild.sh

RUN cd /small_variant_caller \
  && build_testing=${BUILD_TESTING} lint_enable=${LINT_ENABLE} /common/scripts/run-clangd-tidy.sh

FROM ubuntu:24.04 AS minimal

RUN apt-get -qq update && apt-get -qq upgrade && rm -rf /var/lib/apt/lists/*

ENV PATH=/apps:${PATH}

COPY --from=build_bin_stage /small_variant_caller/build/apps/compute_bam_features /apps/compute_bam_features
COPY --from=build_bin_stage /small_variant_caller/build/apps/compute_vcf_features /apps/compute_vcf_features
COPY --from=build_bin_stage /small_variant_caller/build/apps/filter_variants /apps/filter_variants
COPY --from=build_bin_stage /small_variant_caller/build/apps/train_model /apps/train_model
COPY --from=build_bin_stage /small_variant_caller/build/apps/vcf_to_bed /apps/vcf_to_bed
COPY --from=build_bin_stage /small_variant_caller/resources/ /resources/
