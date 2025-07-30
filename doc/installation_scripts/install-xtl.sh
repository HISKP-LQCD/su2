#!/bin/bash
# Installing xtl from source inside a build/ directory


git clone https://github.com/xtensor-stack/xtl
cd xtl/
git checkout d11fb6b


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