# Makefile

all: a

CXX_FLAGS = -I../../yaml-cpp/include/ -I../../yaml-cpp/include/yaml-cpp/
LD_FLAGS = -L../../yaml-cpp/build/ -lyaml-cpp  -L../../yaml-cpp/build/lib/ -lgmock_main  -lgtest_main

a:
	g++ -c $(CXX_FLAGS) a.cpp -o a.o
	g++ $(LD_FLAGS) a.o -o a.exe


clean:
	rm *.exe
