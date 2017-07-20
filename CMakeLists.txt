# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
#             2017 Carsten Urbach <urbach@hiskp.uni-bonn.de>
# Licensed under the GNU Public License version 3

cmake_minimum_required(VERSION 2.8.11)
project(su2 C CXX)
project(su2hmc C CXX)

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
    expsu2.cc
    wilsonloop.cc
    energy_density.cc
    gradient_flow.cc
)

# Create the main program.
add_executable(su2-metropolis main.cc)
target_link_libraries(su2-metropolis su2)

add_executable(su2-hmc hmc.cc)
target_link_libraries(su2-hmc su2)

add_executable(su2-kramers kramers.cc)
target_link_libraries(su2-kramers su2)

add_executable(measure measure.cc)
target_link_libraries(measure su2)


