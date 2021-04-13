# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
#             2017,2021 Carsten Urbach <urbach@hiskp.uni-bonn.de>
# Licensed under the GNU Public License version 3

cmake_minimum_required(VERSION 2.8.11)
project(su2 C CXX)
project(su2hmc C CXX)

add_definitions(--std=c++11)

## here we generate version.hh
## to track the git version of the executables
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
endif(EXISTS "${CMAKE_SOURCE_DIR}/.git")

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hh")
configure_file(
  ${CMAKE_SOURCE_DIR}/version.hh.in
  ${CMAKE_BINARY_DIR}/generated/version.hh
)

include_directories(${CMAKE_SOURCE_DIR}/)
include_directories(${CMAKE_BINARY_DIR}/generated/)

find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -D_USE_OMP_")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D_USE_OMP_")

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
  u1.cc
  su2.cc
  genzsu2.cc
  linearsu2.cc
  print_program_options.cc
  parse_commandline.cc
)

# Create the main program.
find_package(Boost COMPONENTS program_options filesystem system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(su2-metropolis main.cc)
add_executable(genz genz.cc)
add_executable(linear linear.cc)
add_executable(hysteresis hysteresis.cc)
add_executable(u1-metropolis main-u1.cc)
add_executable(test-groups test.cc)
#add_executable(measure measure.cc)

#foreach(target su2-metropolis genz measure)
foreach(target su2-metropolis genz linear test-groups u1-metropolis hysteresis)
  if(Boost_FOUND)
    target_link_libraries(${target} su2 ${Boost_LIBRARIES})
  else()
    target_link_libraries(${target} su2 "-lboost_program_options -lboost_filesystem -lboost_system")
  endif()
endforeach(target)

