all: su2

CXX=g++ -std=c++11
CXXFLAGS=-Wpedantic

su2: main.cc include/su2.hh include/geometry.hh Makefile
	$(CXX) $(CXXFLAGS) main.cc -o su2
