// staggered.hpp
/**
 * Staggered fermions 
 */

#pragma once

#include <array>
#include <vector>
#include <complex>
#include <random>

#include "adjointfield.hh"
#include "gauge_energy.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"

#include "cg.hpp"


namespace staggered {

// eta_{\mu}(x) as in eq. (16) of https://arxiv.org/pdf/2112.14640.pdf
int eta(const std::vector<size_t>& x, const size_t& mu){
  int e=1;
  for(int i=0; i<mu; i++) { e *= std::pow(-1, x[i]); }
  return e;
}

// index of the lattice point given the dimensions
size_t index_from_coords_3d(const std::vector<size_t>& x, const std::vector<size_t>& D) {
  return ((x[0]*D[1] + x[1])*D[2] +x[2])*D[3] + x[3];
}


// vector of staggered "spinors" (no Dirac structure) for all the points of the lattice
template<class Float, class Type>
class spinor_lat_3d {
private:
std::vector<Type> Psi;

public:

spinor_lat_3d(){}
~spinor_lat_3d(){}

Type& operator () (const size_t& i){ return Psi[i]; }
Type operator () (const size_t& i) const { return Psi[i]; }


Type& operator () (const std::vector<size_t>& x, const std::vector<size_t>& dims){
  const size_t i = index_from_coords_3d(x, dims);
  return (*this)[i];
}

Type operator () (const std::vector<size_t>& x, const std::vector<size_t>& dims) const {
  const size_t i = index_from_coords_3d(x, dims);
  return (*this)[i];
}

Type& operator [] (const size_t& i){ return Psi[i];}
Type operator [] (const size_t& i) const { return Psi[i]; }

spinor_lat_3d(const size_t& n){ Psi.resize(n); }

size_t size() const { return Psi.size();}

spinor_lat_3d<Float, Type> operator /(const Type& lambda){
  spinor_lat_3d<Float, Type> phi;
  const int N = (*this).size();
  for (size_t i = 0; i < N; i++) { phi[i] = Psi[i]/lambda; }
  return phi;
}


void operator +=(const spinor_lat_3d<Float, Type>& psi){ (*this) = (*this) + psi; }

Float norm_squared() const{ 
  return complex_dot_product((*this),(*this)).real(); // A^{\dagger} * A is real
}

Float norm() const { return sqrt(this->norm_squared()); }

}; // class spinor_lat_3d


// a + lambda*b
template<class Float, class Type, class Type_lambda>
spinor_lat_3d<Float, Type> a_plus_lambda_b(const spinor_lat_3d<Float, Type>& a, const Type_lambda& lambda, const spinor_lat_3d<Float, Type>& b) 
{
  const int N = a.size();
  spinor_lat_3d<Float, Type> c(N);
  for (size_t i = 0; i < N; i++) { c[i] = a[i] + lambda*b[i]; }
  return c;
}


template<class Float, class Type>
spinor_lat_3d<Float, Type> operator +(const spinor_lat_3d<Float, Type>& a, const spinor_lat_3d<Float, Type>& b){
  return a_plus_lambda_b(a, 1.0, b);
}

template<class Float, class Type>
spinor_lat_3d<Float, Type> operator -(const spinor_lat_3d<Float, Type>& a, const spinor_lat_3d<Float, Type>& b) {
  return a_plus_lambda_b(a, -1.0, b);
}

// change of sign : psi --> -psi
template<class Float, class Type>
spinor_lat_3d<Float, Type> operator -(const spinor_lat_3d<Float, Type>& psi){
  const spinor_lat_3d<Float, Type> v( psi.size() );
  return (v-psi);
}


template<class Float, class Type >
spinor_lat_3d<Float, Type> gaussian_spinor_normalized(const size_t& n, const size_t& seed)
{
  std::normal_distribution<Float> dis{0.0, 1.0 / sqrt(2)};   // e^{-x^2} has sigma=1/sqrt(2)
  std::mt19937 gen_re(seed), gen_im(seed+1);
  spinor_lat_3d<Float, Type> psi_gauss(n);
  Type x, y;
  Type norm2 = 0.0;
  for (size_t i = 0; i < n; i++) {   // lattice points
    x = (dis(gen_re), dis(gen_im));
    psi_gauss[i] = x;
    norm2 += (conj(x)*x).real(); // x^{\dagger}*x is real
  }
  return psi_gauss / sqrt(norm2);
}



// \sum_{i} A_i^{\dagger}*B_i
template<class Float, class Type>
Type complex_dot_product(const spinor_lat_3d<Float, Type>&A, const spinor_lat_3d<Float, Type>&B)
{
  const int N = A.size();
  Type sum;
  for (size_t i = 0; i < N; i++) { sum += conj(A[i]) * B[i]; }
  return sum;
}

template<class Float, class Type>
Type operator *(const spinor_lat_3d<Float, Type>&A, const spinor_lat_3d<Float, Type>&B)
{ return complex_dot_product(A, B); }

template<class Float, class Type>
spinor_lat_3d<Float, Type> operator *(const Type& lambda, const spinor_lat_3d<Float, Type>& psi)
{
  const spinor_lat_3d<Float, Type> v(psi.size());
  return a_plus_lambda_b(v, lambda, psi);
}

template<class Float, class Type>
Float norm(const spinor_lat_3d<Float, Type>& psi){ return psi.norm(); }


template<class Float, class Type>
class matrix_lat_3d {
private:
size_t Lt, Lx, Ly;
std::vector<spinor_lat_3d<Float, Type>> M;

public:
// M^{-1}*psi , found withthe conjugate gradient algorithm
spinor_lat_3d<Float, Type> inv(const spinor_lat_3d<Float, Type> & psi, const Float& tol, const size_t& verb, const size_t& seed)
{
  typedef spinor_lat_3d<Float, Type> LAvector;
  typedef matrix_lat_3d<Float, Type> LAmatrix;
  cg::LinearCG<Float, Type, LAmatrix, LAvector>  LCG((*this), psi);
  const size_t N = psi.size();
  const LAvector phi0 = gaussian_spinor_normalized<Float, Type>(N, seed);

  LCG.solve(phi0, tol, verb);
  return LCG.get_solution();
}

size_t getVolume() const {
  return Lt*Lx*Ly;
}

size_t rows() const { return M.size(); }
size_t cols() const { return M.size(); }

void add_row(const spinor_lat_3d<Float, Type>& psi){ M.push_back(psi); }

Type& operator () (const size_t& i, const size_t& j){ return M[i][j]; }

Type operator () (const size_t& i, const size_t& j) const { return M[i][j]; }


Type& operator () (const std::vector<size_t>& x, const std::vector<size_t>& y){
  const int i = index_from_coords_3d(x, {Lt, Lx, Ly});
  const int j = index_from_coords_3d(y, {Lt, Lx, Ly});
  return M[i][j];
}

Type operator () (const std::vector<size_t>& x, const std::vector<size_t>& y) const {
  const int i = index_from_coords_3d(x, {Lt, Lx, Ly});
  const int j = index_from_coords_3d(y, {Lt, Lx, Ly});
  return M[i][j];
}


};


template<class Float, class Type, class Group>
matrix_lat_3d<Float, Type> get_Qp_3d(gaugeconfig<Group>* U, const Type& m0){  
  const size_t nd = 3;
  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy();//, Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly}; // vector of dimensions

  matrix_lat_3d<Float, Type> Qp;
  std::vector<size_t> x(nd), xp(nd), xm(nd); // updated at each loop step

#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        x = {x0, x1, x2};
        xm = x; xp = x;
        for(size_t mu = 0; mu < nd; mu++) {
          const Float eta_x_mu = eta(x, mu);
          if (eta_x_mu > 0)
          {// only even sites of the lattice
          xm[mu] -= 1; // x - mu
          xp[mu] += 1; // x + mu
          
          if(xp[mu]<dims[mu]){ // not over the right border of the matrix 
            Qp(x, xp) = (1.0/2.0)*eta_x_mu*((*U)(x, mu)); 
          }
          if(xm[mu]>0){ // not over the left border of the matrix
            Qp(x, xm) = (1.0/2.0)*eta_x_mu*((*U)(xm, mu).dagger()); // note the minus sign       
          }
          Qp(x, x)  = (1.0/2.0)*m0;
        
          xm[mu] += 1; // =x again
          xp[mu] -= 1; // =x again
          }

        }
      }
    }
  }
  return Qp;
}


// To be optimized: the matrix is sparse
template<class Float, class Type>
spinor_lat_3d<Float, Type> operator * (const matrix_lat_3d<Float, Type>& M, const spinor_lat_3d<Float, Type>& psi)
{
  const int N = M.getVolume();
  spinor_lat_3d<Float, Type> R(N);
  for (size_t i = 1; i < N; i++) {
          for (size_t j = 0; j < N; j++) {
                  R(i) += M(i, j)*psi(j);
          }
  }
  return R;
}




} // namespace staggered
