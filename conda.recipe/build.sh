#!/usr/bin/env bash
set -ex

sed -i "s^use_MKL=^use_MKL=True^g" setup.cfg
sed -i "s^MKLROOT=^MKLROOT=${CONDA_PREFIX}/lib^g" setup.cfg


# Install the Python package, but without dependencies,
# because Conda takes care of that
$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
sed -i "s^MKLROOT=${CONDA_PREFIX}/lib^MKLROOT=^g" setup.cfg