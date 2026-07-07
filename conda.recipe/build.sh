#!/usr/bin/env bash
set -ex

export CC=mpicc
export FC=mpifort

OPENMP_ARGS=()
OPENMP_VARIANT="${openmp:-${OPENMP:-false}}"
if [[ "${OPENMP_VARIANT,,}" == "true" ]]; then
    OPENMP_ARGS+=(--config-settings=setup-args="-Dopenmp=true")
fi

$PYTHON -m pip install . -v --no-build-isolation \
    --config-settings=setup-args="-Dlink_args_blas=-L"$PREFIX"/lib,-lopenblas" \
    --config-settings=setup-args="-Dincdirs="$PREFIX"/include" \
    --config-settings=setup-args="-Dlibdirs="$PREFIX"/lib" \
    "${OPENMP_ARGS[@]}" \
    --config-settings=setup-args="-Dbuildtype=release"
