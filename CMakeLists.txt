# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the GNU Public License version 3

cmake_minimum_required(VERSION 2.8.11)
project(su2 C CXX)

add_definitions(--std=c++11)
# There is one header file which is in a subdirectory. This add the `-i` to the
# compiler.
#include_directories("include")

# We could just have two `add_executable` blocks listing almost all source
# files. However, CMake would then compile the source files twice, once for
# each executable. The rationale is that different flags can be used for
# different executables. Here we want to re-use this. One can either track the
# generated objects files or create a simple static library. The latter is done
# here.
add_library(su2
    gauge_energy.cc
    gaugeconfig.cc
    get_staples.cc
    random_gauge_trafo.cc
    su2.cc
    sweep.cc
)

# Create the main program.
add_executable(su2-main main.cc)
target_link_libraries(su2-main su2)

# Create the test driver.
add_executable(su2-test test.cc)
target_link_libraries(su2-test su2)
