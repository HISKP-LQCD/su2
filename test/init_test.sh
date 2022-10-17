# init_test.sh

cd ..
d1=$(pwd -P)
cd test

YAML_SRC_PATH=$d1/external/yaml-cpp/

rm -rf build # removing the previous build/ directory
mkdir build # creating a new one

cp config.yaml build

cd build
cmake ../ -D YAML_SRC=$YAML_SRC_PATH
