#!/usr/bin/env bash
set -ex

CC=mpicc
FC=mpifort
# Install the Python package, but without dependencies,
# because Conda takes care of that
$PYTHON -m pip install . --no-build-isolation \
--config-settings=setup-args="-Dblas_lib=['openblas','blas','lapack']"