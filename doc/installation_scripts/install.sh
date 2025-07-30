#!/bin/bash

echo "Configuring the project with cmake"

# it's a good idea to first clean up the old cmake configuration
rm -r CMakeCache.txt CMakeFiles/ cmake_install.cmake generated/ Makefile

bdir=/path/to/your/build/directory/

# specify these dependencies paths only if you are not loading them already as modules
yaml_cpp=/path/to/yaml-cpp/installation/directory/
xtl=/path/to/xtl/installation/directory/
xtensor=/path/to/xtensor/installation/directory/
eigen=/path/to/eigen/installation/directory/


cmake  \
  -D CMAKE_BUILD_TYPE=${btype} \
  -S /path/to/source/code/ \
  -B ${d1} \
  -D CMAKE_PREFIX_PATH="$yaml_cpp;$xtl;$xtensor;$eigen"
  

make -j$1 install
