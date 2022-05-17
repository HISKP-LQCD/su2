cmake_minimum_required(VERSION 3.10)

# set the project name
project(test VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executables

# conjugate gradient algorithms
add_executable(CG.exe CG.cpp)
add_executable(BiCGStab.exe BiCGStab.cpp)

#yaml input file parsing

set(yaml_src_dir ${YAML_SRC})
set(
    inc_dir  
    ${yaml_src_dir}/include/
    ${yaml_src_dir}/include/yaml-cpp/
    )

add_executable(yaml.exe yaml.cpp)
target_include_directories(yaml.exe PUBLIC ${inc_dir})

set(
    link_dir  
    ${yaml_src_dir}/build/
    ${yaml_src_dir}/build/lib/
)


target_link_directories(yaml.exe PUBLIC ${link_dir})
target_link_libraries(yaml.exe libyaml-cpp.a)

