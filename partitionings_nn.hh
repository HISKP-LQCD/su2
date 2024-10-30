/**
 * @file partitionings_nn.hh
 * @author Sebastian Müller (s6sbmuel@uni-bonn.de)
 * @brief class representation of partition elements for the neirest neighbor algorithm
 * @version 0.1
 * @date 2024-10-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once 
#include<iomanip>
#include<stdexcept>
#include<complex>

#include<cmath>
#include<cassert>
#include<iostream>
#include<vector>
#include "accum_type.hh"
#include "dagger.hh"
#include<algorithm>
#include "su2.hh"
using Complex = std::complex<double>;


/**
 * @brief the partitionings are digitised in the fundamental SU(2) representation
 * following arXiv:2201.09625v1 [hep-lat]
 * 
 */

class _partitioning_nn {
    public:
    const size_t N_c = 2;
    inline static std::vector<double> point0; // vector to save x0
    inline static std::vector<double> point1; // vector to save x1
    inline static std::vector<double> point2; // vector to save x2
    inline static std::vector<double> point3; // vector to save x3
    inline static std::vector<double> weights; // vector to save the weights of each element
    inline static std::vector< std::vector <size_t> > nn_lookup; // lookup table for the neighrest neighbor of each element (length of contained vectors may vary)
    size_t index; // index identifiying the elment of the partitioning
    explicit _partitioning_nn(): index(0){} // create a new instance with index 0
    explicit _partitioning_nn(size_t i): index(i){} // create a new instance with given index
    _partitioning_nn(const _partitioning_nn &U): index(U.index){} // copy an instance by copiny the index
    inline double getweight() {return weights[index];} // get the weight of the element
    inline std::vector <size_t> getneigborindeces() const {return _partitioning_nn::nn_lookup[index];} // get a vector with the neighrest neighbors
    inline su2 getsu2() const {return su2(Complex(point0[index], point1[index]), Complex(point2[index], point3[index]));} // get the su2 matrix of the partitioning
    friend inline su2 operator+(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline su2 operator-(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline  su2 operator*(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline su2 operator*(const Complex &U1, const _partitioning_nn &U2);
    inline size_t getindex() const {return (index);} // get the index of the partitioning
    inline void operator=(const _partitioning_nn &U){ index = U.getindex();}
    explicit _partitioning_nn( const su2 &U){ // get the closest partitioning to an su2 matrix
        su2 helpmatrix = U;
        helpmatrix.restoreSU();
         
        double min_distance = 1000000000000;
        size_t min_index = 0; 
        size_t i = 0;
        // this is probably the fastest way without assuming any other information about the partitioning 
        while  ( i < point0.size()){
            double distance = pow(point0[i] - std::real(helpmatrix.geta()), 2.) + pow(point1[i] - std::imag(helpmatrix.geta()), 2.) + pow(point2[i] - std::real(helpmatrix.getb()), 2.) + pow(point3[i] - std::imag(helpmatrix.getb()), 2.);
            
            if (distance <= min_distance){
                min_index = i;
                min_distance = distance;
                
            }
            i = i +1;
        } 

        index = min_index;       
    }
    inline Complex geta() const { // get a from the su2 matrix
        Complex a(point0[index], point1[index]);
        return a;
    }
    inline su2 operator*=(const su2 &U1){
        su2 helpmatrix = getsu2();
        return helpmatrix*U1;
    }
    inline su2 operator*=(const _partitioning_nn &U1){
        su2 helpmatrix1 = getsu2();
        su2 helpmatrix2 = U1.getsu2();
        return helpmatrix1*helpmatrix2;
    }
    inline su2 operator+=(const su2 &U1){
        su2 helpmatrix = getsu2();
        return helpmatrix+U1;
    }
    inline su2 operator+=(const _partitioning_nn &U1){su2 helpmatrix1 = getsu2();
        su2 helpmatrix2 = U1.getsu2();
        return helpmatrix1+helpmatrix2;
        }
    inline Complex getb() const { //get b from the su2 matrix
        Complex b(point2[index], point3[index]);
        return b;
    }
    su2 round(size_t n){ // uses su2 round to get the rounded partioning
        su2 helpmatrix(Complex(point0[index], point1[index]), Complex(point2[index], point3[index]));
        return helpmatrix.round(n);
    }
    inline su2 operator+=(_partitioning_nn &U){
        return getsu2()*U.getsu2();
    }
    inline su2 set_to_identity(){ //returns an su2 identity matrix since the identity is not necessarily part of the partitioning (eg. for Fibonacci)
        su2 helpmatrix;
        helpmatrix.set_to_identity();
        return helpmatrix;
    }
    void set(const size_t _i){index = _i;} // set the partitioning to a specific element
    inline su2 dagger() const {su2 helpmatrix = getsu2();
    return helpmatrix.dagger();}
    inline double retrace() {return 2*point0[index];}
    Complex det() {return Complex(1, 0);} //the determinant of a partioning is 1 by construction
    void restoreSU() {} // a partitioning is an SU(2) element by construction
    void print(){
        std::cout << "--------------------\n";
        std::cout << index << "\n";
        std::cout << point0[index] << " " << point1[index] << " " << point2[index] << " " << point3[index] << "\n";
        std::cout << "--------------------\n";
    }
};

inline double retrace( _partitioning_nn const &U){
    return 2*_partitioning_nn::point0[U.getindex()];
}
inline Complex trace(_partitioning_nn const &U){
    return Complex(2*_partitioning_nn::point0[U.getindex()], 0);
}
inline su2 dagger(const _partitioning_nn &u){ //dagger is not necessarily an element of the partitioning
    su2 helpmatrix = u.getsu2();
    return helpmatrix.dagger();
}
inline su2 traceless_antiherm(const _partitioning_nn &x) { //not necessarily an element of the partitioning. It thus returns an su2 element
    su2 helpmatrix = x.getsu2();
    return traceless_antiherm(helpmatrix);
}

inline su2 operator*(const _partitioning_nn &U1, const _partitioning_nn &U2){
    return U1.getsu2()*U2.getsu2();
}
inline su2 operator*(const su2 &U1, const _partitioning_nn &U2){
    return U1 * U2.getsu2();
}
inline su2 operator*(const _partitioning_nn &U1, const su2 &U2){
    return U1.getsu2() * U2;
}
inline su2 operator*(const Complex &U1, const _partitioning_nn &U2){
    return U1*U2.getsu2();
}
inline su2 operator*(const _partitioning_nn &U1, const Complex &U2){
    return U1.getsu2()*U2;
}
inline su2 operator+(const _partitioning_nn &U1, const _partitioning_nn &U2){
    return U1.getsu2() + U2.getsu2();
}
inline su2 operator+(const _partitioning_nn &U1, const su2 &U2){
    return U1.getsu2()*U2;
}
inline su2 operator+(const su2 &U1, const _partitioning_nn &U2){
    return U1 * U2.getsu2();
}

inline su2 operator-(const _partitioning_nn &U1, const _partitioning_nn &U2){
    return U1.getsu2() - U2.getsu2();
}
inline su2 operator-(const _partitioning_nn &U1, const su2 &U2){
    return U1.getsu2() - U2;
}
inline su2 operator-(const su2 &U1, const _partitioning_nn &U2){
    return U1 - U2.getsu2();
}

template <> struct accum_type<_partitioning_nn> {
    typedef su2 type;
};
using partitioning_nn = _partitioning_nn;