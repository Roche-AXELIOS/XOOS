ARG build_type=Release

FROM ubuntu:26.04 AS builder-base

ENV PATH=/ci/bin:/ci/sccache:/ci/.venv/bin:${PATH}

COPY --from=common /ci/ci-build-env-base.sh /ci/bin/

ARG SCCACHE_ENABLE=false
RUN --mount=type=cache,target=/root/.cache/pip \
  --mount=type=cache,target=/var/cache/apt,sharing=locked \
  --mount=type=cache,target=/var/lib/apt,sharing=locked \
  sccache_enable=${SCCACHE_ENABLE} python_venv_dir=/ci/.venv ci-build-env-base.sh

COPY --from=common /ci/ci-build-env-conan.sh /ci/bin/
COPY --from=vendor . /ci/vendor
COPY --from=common . /ci/common

WORKDIR /ci/app
COPY conanfile.txt .
COPY CMakeLists.txt .

ARG SCCACHE_BUCKET SCCACHE_REGION
ARG BUILD_DEPENDENCIES=missing CONANCENTER_ENABLED=true CUSTOM_REMOTE_URL CUSTOM_REMOTE_UPLOAD_PATTERN
RUN --mount=type=secret,id=AWS_ACCESS_KEY_ID \
  --mount=type=secret,id=AWS_SECRET_ACCESS_KEY \
  --mount=type=secret,id=AWS_SESSION_TOKEN \
  --mount=type=secret,id=CUSTOM_REMOTE_USERNAME \
  --mount=type=secret,id=CUSTOM_REMOTE_PASSWORD \
  SCCACHE_BUCKET="${SCCACHE_BUCKET}" \
  SCCACHE_REGION="${SCCACHE_REGION}" \
  AWS_ACCESS_KEY_ID="$([ -f /run/secrets/AWS_ACCESS_KEY_ID ] && cat /run/secrets/AWS_ACCESS_KEY_ID)" \
  AWS_SECRET_ACCESS_KEY="$([ -f /run/secrets/AWS_SECRET_ACCESS_KEY ] && cat /run/secrets/AWS_SECRET_ACCESS_KEY)" \
  AWS_SESSION_TOKEN="$([ -f /run/secrets/AWS_SESSION_TOKEN ] && cat /run/secrets/AWS_SESSION_TOKEN)" \
  root_dir=/ci \
  build_dependencies="${BUILD_DEPENDENCIES}" \
  custom_remote_url="${CUSTOM_REMOTE_URL}" \
  custom_remote_username="$([ -f /run/secrets/CUSTOM_REMOTE_USERNAME ] && cat /run/secrets/CUSTOM_REMOTE_USERNAME)" \
  custom_remote_password="$([ -f /run/secrets/CUSTOM_REMOTE_PASSWORD ] && cat /run/secrets/CUSTOM_REMOTE_PASSWORD)" \
  custom_remote_upload_pattern="${CUSTOM_REMOTE_UPLOAD_PATTERN}" \
  conancenter_enabled="${CONANCENTER_ENABLED}" \
  ci-build-env-conan.sh

FROM builder-base AS builder-app
ARG SCCACHE_BUCKET SCCACHE_REGION

COPY --from=common /ci/ci-app.sh /ci/bin/
COPY . /ci/app/

ARG BUILD_TESTING=OFF CODE_COVERAGE_ENABLE=OFF build_type VERSION=unknown
RUN --mount=type=secret,id=AWS_ACCESS_KEY_ID \
  --mount=type=secret,id=AWS_SECRET_ACCESS_KEY \
  --mount=type=secret,id=AWS_SESSION_TOKEN \
  SCCACHE_BUCKET="${SCCACHE_BUCKET}" \
  SCCACHE_REGION="${SCCACHE_REGION}" \
  AWS_ACCESS_KEY_ID="$([ -f /run/secrets/AWS_ACCESS_KEY_ID ] && cat /run/secrets/AWS_ACCESS_KEY_ID)" \
  AWS_SECRET_ACCESS_KEY="$([ -f /run/secrets/AWS_SECRET_ACCESS_KEY ] && cat /run/secrets/AWS_SECRET_ACCESS_KEY)" \
  AWS_SESSION_TOKEN="$([ -f /run/secrets/AWS_SESSION_TOKEN ] && cat /run/secrets/AWS_SESSION_TOKEN)" \
  build_testing=${BUILD_TESTING} \
  code_coverage_enable=${CODE_COVERAGE_ENABLE} \
  build_type=${build_type} \
  version=${VERSION} \
  ci-app.sh

ARG LINT_ENABLE=ON
COPY --from=common /ci/ci-clangd-tidy.sh /ci/bin/
RUN build_testing=${BUILD_TESTING} lint_enable=${LINT_ENABLE} ci-clangd-tidy.sh

FROM cgr.dev/chainguard/wolfi-base:latest

RUN apk add --no-cache bash procps

ENV PATH="/apps:${PATH}"

ARG build_type
COPY --from=builder-app /ci/build/apps/${build_type}/compute_bam_features /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/compute_vcf_features /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/filter_variants /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/train_model /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/vcf_to_bed /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/generate_pon /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/phase_variants /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/detect_sample_type /apps/
COPY --from=builder-app /ci/build/apps/${build_type}/features_to_vcf /apps/
COPY --from=builder-app /ci/app/resources/ /resources/
COPY LICENSE.md /LICENSE.md
