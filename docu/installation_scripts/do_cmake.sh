#!/bin/bash

echo "Configuring the project with cmake"

bash ./clean.sh # clean old cmake configuration

bdir=$(pwd -P) # build directory

yaml_cpp="/path/to/yaml-cpp/installation/directory/"
xtensor="/path/to/xtensor/installation/directory/"

for btype in debug release; do
  echo "Build type: ${btype}"
  
  d1=${bdir}/${btype}
  mkdir -p $d1
  
  cmake  \
  -D CMAKE_BUILD_TYPE=${btype} \
  -S /path/to/source/code/ \
  -B ${d1} \
  -D CMAKE_PREFIX_PATH="$yaml_cpp;$xtensor"
  
  cd ${bdir}
done
