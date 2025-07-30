#!/bin/bash
# Installing xtensor from source inside a build/ directory


git clone https://github.com/xtensor-stack/xtensor
cd xtensor
git checkout 8c0a484f

mkdir build/
cd build

xtl_dir=/path/to/your/xtl/installation/directory/

sdir=$(realpath ../)
bdir=$(realpath ./)
mkdir install_dir/
idir=$(realpath ./install_dir/)

cmake \
  -S ${sdir} \
  -B ${bdir} \
  -D CMAKE_INSTALL_PREFIX=${idir} \
  -D CMAKE_PREFIX_PATH=${xtl}

make -j$(nproc) install