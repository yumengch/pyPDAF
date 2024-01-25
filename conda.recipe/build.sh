#!/usr/bin/env bash
set -ex


# Install the Python package, but without dependencies,
# because Conda takes care of that
$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv