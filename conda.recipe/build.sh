#!/usr/bin/env bash
set -ex

export CC=mpicc
export FC=mpifort

OPENMP_ARGS=()
OPENMP_VARIANT="${openmp:-${OPENMP:-false}}"
OPENMP_VARIANT_LOWER=$(printf '%s' "$OPENMP_VARIANT" | tr '[:upper:]' '[:lower:]')
if [[ "$OPENMP_VARIANT_LOWER" == "true" ]]; then
    OPENMP_ARGS+=(--config-settings=setup-args="-Dopenmp=true")
fi

$PYTHON -m pip install . -v --no-build-isolation \
    --config-settings=setup-args="-Dlink_args_blas=-L"$PREFIX"/lib,-lopenblas" \
    --config-settings=setup-args="-Dincdir_blas="$PREFIX"/include" \
    "${OPENMP_ARGS[@]}" \
    --config-settings=setup-args="-Dbuildtype=release"
