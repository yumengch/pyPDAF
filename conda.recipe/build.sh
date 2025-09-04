#!/usr/bin/env bash
set -ex

CC=mpicc
FC=mpifort

$PYTHON -m pip install . --no-build-isolation \
    --config-settings=setup-args="-Dblas_lib=['openblas','blas','lapack']"