/**
 * @file partitionings.hh
 * @author Sebastian Mueller (s6sbmuel@uni-bonn.de)
 * @brief class representing the SU(2)-partitioning (in a fundamental representation)
 * @version 0.1
 * @date 2024-10-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */




#pragma once

#include<fstream>
#include<sstream>
#include<cstdio>
#include<iomanip>
#include<stdexcept>
#include<complex>

//#include<math.h>
#include<cmath>
#include<cassert>
#include<iostream>
#include<vector>
#include "accum_type.hh"
#include "dagger.hh"
#include<algorithm>
using Complex = std::complex<double>;

class _partitioning {

    public:
    const size_t N_c = 2;
    inline static std::vector<size_t> dagger_vector; // index of element closest to daggered
    inline static std::vector<double> point0 ; // 0 coordinate of point on hypersphere
    inline static std::vector<double> point1 ; // 1 coordinate of point on hypersphere
    inline static std::vector<double> point2 ; // 2 coordinate of point on hypersphere
    inline static std::vector<double> point3 ; // 3 coordinate of point on hypersphere
    inline static std::vector<double> weights; //the weights of each gauge configuration
    inline static std::vector<double> distance_to_identity ; // the euclidean distance of every point on the hypersphere to the idenity
    inline static std::vector< std::vector<size_t> > multiplication_lookup_table; // a lookup table for the indeces of the results of multiplication
    inline static std::vector< std::vector<size_t> > addition_lookup_table ; // a lookup table for the indeces of the results of addition
    size_t index;
    Complex multiplicator;
    explicit _partitioning(): index(0), multiplicator(Complex(1, 0)) {}
    explicit _partitioning(size_t i): index(i) , multiplicator(Complex(1,0)){}
    explicit _partitioning(size_t i, Complex lambda): index(i), multiplicator(lambda) {}
    _partitioning(const _partitioning &U) : index(U.index), multiplicator(U.multiplicator)  {}

    friend inline _partitioning operator+(const _partitioning &U1, const _partitioning &U2);
    friend inline _partitioning operator-(const _partitioning &U1, const _partitioning &U2);
    friend inline _partitioning operator*(const _partitioning &U1, const _partitioning &U2);
    friend inline _partitioning operator*(const Complex &U1, const _partitioning &U2);
    inline size_t getindex() const {return (index);}
    inline Complex getmultiplicator() const {return (multiplicator);}
    inline void operator=(const _partitioning &U){
        index = U.getindex();
        multiplicator = U.getmultiplicator();
    }
    _partitioning round(size_t n) {
        double dn = n;
        return 
    _partitioning(index, 
    Complex(std::round(std::real(multiplicator)*dn), std::round(std::imag(multiplicator)*dn))/dn);}
    _partitioning &operator*=(const _partitioning &U1){
        this -> index = multiplication_lookup_table[this->index][U1.index];
        this -> multiplicator = (this->multiplicator)*U1.getmultiplicator();
        return *this; 
    }
    inline Complex geta() const {
        Complex a(point0[index], point1[index]);
        return a*multiplicator;
    }
    inline Complex getb() const {
        Complex b(point2[index], point3[index]);
        return b*multiplicator;
    }
    inline void operator+=(const _partitioning &U){
        Complex new_determinant = (geta() + U.geta())*std::conj(geta() + U.geta()) + (getb() + U.getb())*std::conj((getb() + U.getb()));
        multiplicator = sqrt(new_determinant);
        index = addition_lookup_table[index][U.index];

    }
    void set(const size_t _index){index = _index;}
    void set(const size_t _index, const Complex _multiplicator){
    index = _index;
    multiplicator = _multiplicator;
    }
    void set_to_identity(){
        std::vector<double>::iterator result = std::min_element(distance_to_identity.begin(), distance_to_identity.end());
        index = std::distance(distance_to_identity.begin(), result);
        multiplicator = Complex(1,0);
    }
    inline _partitioning dagger() const { return (_partitioning(dagger_vector[index], std::conj(multiplicator)));}
    Complex det(){return (multiplicator);}
    void restoreSU() {multiplicator = Complex(1, 0);}
    void print(){
        std::cout << "--------------------\n";
        std::cout << index << " \n";
        std::cout << multiplicator << " \n";
        std::cout << point0[index] << " " << point1[index] << " " << point2[index] << " " << point3[index] << " \n";
        std::cout << "--------------------\n";
    }

};

inline double retrace (_partitioning const &U) {
    double a = std::real(U.geta());
    return (2*a);
}
inline Complex trace(_partitioning const &U) {
    double a = std::real(U.geta());
    return (Complex(2*a, 0));
}

template <> inline _partitioning dagger(const _partitioning &u) {
    size_t resindex = _partitioning::dagger_vector[u.getindex()];
    _partitioning dag(resindex, std::conj(u.getmultiplicator()));
    return dag;
    //return (_partitioning(_partitioning.dagger_vector[u.getindex]));
}

inline _partitioning operator*(const _partitioning &U1, const _partitioning &U2){
    _partitioning res;
    res.set(_partitioning::multiplication_lookup_table[U1.getindex()][U2.getindex()], U1.getmultiplicator()*U2.getmultiplicator());
    return (res);
}
inline _partitioning operator+(const _partitioning &U1, const _partitioning &U2){
    _partitioning res;
    Complex new_determinant = (U1.geta() + U2.geta())*std::conj((U1.geta() + U2.geta())) + (U1.getb() + U2.getb())*std::conj((U1.getb() + U2.getb())); 
    res.set(_partitioning::addition_lookup_table[U1.getindex()][U2.getindex()], sqrt(new_determinant));
    return (res);
}

inline _partitioning operator-(const _partitioning &U1, const _partitioning &U2){
    _partitioning res;
    
    Complex new_determinant = (U1.geta() - U2.geta())*std::conj((U1.geta() - U2.geta())) + (U1.getb() - U2.getb())*std::conj((U1.getb() - U2.getb())); 
    std::vector<size_t>::iterator helpvalue = std::find(_partitioning::addition_lookup_table[U2.getindex()].begin(), _partitioning::addition_lookup_table[U2.getindex()].end(), U1.index);
    res.set(std::distance(_partitioning::addition_lookup_table[U2.getindex()].begin(), helpvalue), sqrt(new_determinant));
    return (res); 
}

inline _partitioning operator*(const Complex &U1, const _partitioning &U2){
    _partitioning res;
    res.set(U2.getindex(), U2.getmultiplicator()*U1);
    return (U2);
}

inline _partitioning operator*(const _partitioning &U1, const Complex &U2){
    _partitioning res;
    res.set(U1.getindex(), U1.getmultiplicator()*U2);
    return (U1);
}

using partitioning_accum = _partitioning;

inline _partitioning accum_to_Group(const partitioning_accum &x){
    _partitioning res;
    res.set(x.getindex(), Complex(1, 0));
    return x;
}

template <> struct accum_type<_partitioning> {
  typedef _partitioning type;
};

/*static std::vector<size_t> dagger_vector = {0}; // index of element closest to daggered
    static std::vector<double> point0 = {0}; // 0 coordinate of point on hypersphere
    static std::vector<double> point1 =  {0} ; // 1 coordinate of point on hypersphere
    static std::vector<double> point2 = {0 }; // 2 coordinate of point on hypersphere
    static std::vector<double> point3 = {0}; // 3 coordinate of point on hypersphere
    static std::vector<double> distance_to_identity  = {0 }; // the euclidean distance of every point on the hypersphere to the idenity
    static std::vector< std::vector<size_t> > multiplication_lookup_table = {{0}}; // a lookup table for the indeces of the results of multiplication
    static std::vector< std::vector<size_t> > addition_lookup_table = {{0}}; // a lookup table for the indeces of the results of addition
  */  
using partitioning = _partitioning;
