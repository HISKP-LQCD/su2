#!/bin/bash

echo "Generating the executables"

bdir=$(pwd -P) # build directory

for btype in debug release; do
  echo "Build type: ${btype}"
  
  d1=${bdir}/${btype}
  cd $d1
  
  make -j$1
  
  cd ${bdir}
done