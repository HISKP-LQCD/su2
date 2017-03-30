#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the GNU Public License version 3

set -e
set -u

if [[ -d build ]]; then
    echo "Directory build exists, this script will do nothing."
    exit 1
fi

mkdir -p build/debug
mkdir -p build/release

pushd build/debug
cmake ../.. -DCMAKE_BUILD_TYPE=Debug
popd

pushd build/release
cmake ../.. -DCMAKE_BUILD_TYPE=Release
popd
