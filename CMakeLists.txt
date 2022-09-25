# Copyright © 2022 Carsten Urbach <urbach@hiskp.uni-bonn.de>
#             2022 Simone Romiti <simone.romiti@uni-bonn.de>
# Licensed under the GNU Public License version 3

cmake_minimum_required(VERSION 3.22.1)
project(su2 C CXX)

set (CMAKE_CXX_STANDARD 11) # c++ standard

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
else()
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
endif()

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hh")
configure_file(
  ${CMAKE_SOURCE_DIR}/version.hh.in
  ${CMAKE_BINARY_DIR}/generated/version.hh
)

include_directories(${CMAKE_SOURCE_DIR}/)
include_directories(${CMAKE_SOURCE_DIR}/include/)
include_directories(${CMAKE_BINARY_DIR}/generated/)


# openMP threads parallelization
find_package(OpenMP)


# setting C and CXX compilation flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -D_USE_OMP_")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D_USE_OMP_")

include_directories(${YAML_SRC_DIR}/include/ ${YAML_SRC_DIR}/include/yaml-cpp/)


# We could just have two `add_executable` blocks listing almost all source
# files. However, CMake would then compile the source files twice, once for
# each executable. The rationale is that different flags can be used for
# different executables. Here we want to re-use this. One can either track the
# generated objects files or create a simple static library. The latter is done
# here.
add_library(${PROJECT_NAME}
    exp_gauge.cc
    print_program_options.cc
    parse_commandline.cc
    parse_input_file.cc
)

# Create the main program.
find_package(Boost COMPONENTS program_options filesystem system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(su2-metropolis main.cc)
add_executable(u1-metropolis metropolis-u1.cc)
add_executable(su2-hmc hmc.cc)
add_executable(u1-hmc hmc-u1.cc)
add_executable(su2-kramers kramers.cc)
add_executable(su2-measure measure.cc)
add_executable(u1-measure measure-u1.cc)
add_executable(test-groups test.cc)
add_executable(scaling scaling.cc)
add_executable(try tryreduce.cc)

set(link_dir  ${YAML_BINARY_DIR})

target_link_directories(${CMAKE_PROJECT_NAME} PUBLIC ${link_dir})

foreach(target su2-hmc u1-hmc su2-metropolis u1-metropolis su2-kramers su2-measure u1-measure test-groups scaling try)
  target_link_libraries(${target} su2 libyaml-cpp.a)
  target_link_directories(${target} PUBLIC ${link_dir})

  if(Boost_FOUND)
    target_link_libraries(${target} su2 ${Boost_LIBRARIES})
  else()
    target_link_libraries(${target} su2 "-lboost_program_options -lboost_filesystem -lboost_system")
  endif()
endforeach(target)

