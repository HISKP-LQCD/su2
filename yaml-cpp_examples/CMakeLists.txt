cmake_minimum_required(VERSION 3.10)

# set the project name
project(yce VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
#set(yaml_src_dir /home/simone/Documents/bonn_work/U1_project/yaml-cpp/)

set(yaml_src_dir ${YAML_SRC})
set(
    inc_dir  
    ${yaml_src_dir}/include/
    ${yaml_src_dir}/include/yaml-cpp/
    )
include_directories(${inc_dir})

set(
    link_dir  
    ${yaml_src_dir}/build/
    ${yaml_src_dir}/build/lib/
)
link_directories(${link_dir})
add_executable(yce main.cpp)
target_link_libraries(yce libyaml-cpp.a)


# configured header file
configure_file(yce.h.in yce.h)

target_include_directories(yce PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )