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

class _partitioning_nn {
    public:
    const size_t N_c = 2;
    inline static std::vector<double> point0;
    inline static std::vector<double> point1;
    inline static std::vector<double> point2;
    inline static std::vector<double> point3;
    inline static std::vector<double> weights;
    inline static std::vector< std::vector <size_t> > nn_lookup;
    size_t index;
    explicit _partitioning_nn(): index(0){}
    explicit _partitioning_nn(size_t i): index(i){}
    _partitioning_nn(const _partitioning_nn &U): index(U.index){}
    inline double getweight() {return weights[index];}
    inline std::vector <size_t> getneigborindeces() const {return _partitioning_nn::nn_lookup[index];}
    inline su2 getsu2() const {return su2(Complex(point0[index], point1[index]), Complex(point2[index], point3[index]));}
    friend inline su2 operator+(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline su2 operator-(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline  su2 operator*(const _partitioning_nn &U1, const _partitioning_nn &U2);
    friend inline su2 operator*(const Complex &U1, const _partitioning_nn &U2);
    inline size_t getindex() const {return (index);}
    inline void operator=(const _partitioning_nn &U){ index = U.getindex();}
    explicit _partitioning_nn( const su2 &U){
        su2 helpmatrix = U;
        helpmatrix.restoreSU();
        //std::cout << "problem ran \n"; 
        double min_distance = 1000000000000;
        size_t min_index = 0; //TODO: Do this without looping over everything...
        size_t i = 0;
        //std::cout << point0.size() << "point 0 size \n";
        while  ( i < point0.size()){
            double distance = pow(point0[i] - std::real(helpmatrix.geta()), 2.) + pow(point1[i] - std::imag(helpmatrix.geta()), 2.) + pow(point2[i] - std::real(helpmatrix.getb()), 2.) + pow(point3[i] - std::imag(helpmatrix.getb()), 2.);
            //std::cout << "distance " << distance << "\n";
            if (distance <= min_distance){
                min_index = i;
                min_distance = distance;
                //std::cout << "got here \n";
            }
            i = i +1;
        } 

        index = min_index;
        //std::cout << "min index " << min_index << "\n";
        //std::cout << "min distance " << min_distance << "\n";
        //return _partitioning_nn(min_index);
        //std::vector<double> helpvector;
        //std::transform(point0.begin(), point0.end(), std::back:inserter(helpvector), &)       
    }
    inline Complex geta() const {
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
    inline Complex getb() const {
        Complex b(point2[index], point3[index]);
        return b;
    }
    su2 round(size_t n){
        su2 helpmatrix(Complex(point0[index], point1[index]), Complex(point2[index], point3[index]));
        return helpmatrix.round(n);
    }
    inline su2 operator+=(_partitioning_nn &U){
        return getsu2()*U.getsu2();
    }
    inline su2 set_to_identity(){
        su2 helpmatrix;
        helpmatrix.set_to_identity();
        return helpmatrix;
    }
    void set(const size_t _i){index = _i;}
    inline su2 dagger() const {su2 helpmatrix = getsu2();
    return helpmatrix.dagger();}
    inline double retrace() {return 2*point0[index];}
    Complex det() {return Complex(1, 0);}
    void restoreSU() {}
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
inline su2 dagger(const _partitioning_nn &u){
    su2 helpmatrix = u.getsu2();
    return helpmatrix.dagger();
}
inline su2 traceless_antiherm(const _partitioning_nn &x) {
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
//inline su2 &su2::operator*=(const _partitioning_nn &U1){
//    *this = *this * U1.getsu2();
//    return *this;
 // }
template <> struct accum_type<_partitioning_nn> {
    typedef su2 type;
};
using partitioning_nn = _partitioning_nn;