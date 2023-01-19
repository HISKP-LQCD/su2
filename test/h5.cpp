// h5.cpp
// example taken from here: https://github.com/BlueBrain/HighFive

#include <iostream>
#include <highfive/H5Easy.hpp>

int main() {
    H5Easy::File file("example.h5", H5Easy::File::Overwrite);

    int A = 2;
    H5Easy::dump(file, "/path/to/A", A);

    A = H5Easy::load<int>(file, "/path/to/A");
    std::cout << "A = " << A << "\n";
}