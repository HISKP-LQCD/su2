#!/bin/bash
# Installing xtensor from source inside a build/ directory


git clone https://github.com/jbeder/yaml-cpp
cd yaml-cpp

mkdir build/
cd build

sdir=$(realpath ../)
bdir=$(realpath ./)
mkdir install_dir/
idir=$(realpath ./install_dir/)

cmake \
  -S ${sdir} \
  -B ${bdir} \
  -D CMAKE_INSTALL_PREFIX=${idir}

make -j$(nproc) install