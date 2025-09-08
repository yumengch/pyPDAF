#!/usr/bin/env bash
set -ex

CC=mpicc
FC=mpifort

$PYTHON -m pip install . -v --no-build-isolation \
    --config-settings=setup-args="-Dblas_lib=['openblas',]" \
    --config-settings=setup-args="-Dincdirs="$PREFIX"/include" \
    --config-settings=setup-args="-Dlibdirs="$PREFIX"/lib" \
    --config-settings=setup-args="-Dbuildtype=release"