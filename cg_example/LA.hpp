// LA.hpp
// containes for linear algebra: vectors and matrices

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

// typedef std::vector<std::vector<T> > matrix;
typedef std::complex<double> Real;

template<class T>
class matrix
{

private:
std::vector<std::vector<T> > M;

public:

matrix(){}
~matrix(){}

matrix(const std::vector<std::vector<T>>& m)
{
    const int nr = m.size();
    const int nc = m[0].size();

    M.resize(nr); // saving the space for enough rows
    for(int i=0; i<nr; ++i)
    {
        M[i].resize(nc); // saving the space for enough column entries
        for(int j=0; j<nc; ++j)
        {
            M[i][j] = m[i][j];
        }
    }
}

matrix(const size_t& nr, const size_t& nc){ this->resize(nr, nc); };

// storing enough space for the matrix elements
void resize(const size_t& nr, const size_t& nc)
{
    M.resize(nr);
    for(int i=0; i<nr; ++i){ M[i].resize(nc); }
}

size_t rows() const { return M.size(); }
size_t cols() const { return M[0].size(); }

void add_row (const std::vector<T>& v){ M.push_back(v); }

std::vector<T>& operator [] (const int& i){ return M[i];}
std::vector<T> operator [] (const int& i) const { return M[i];}

std::vector<std::vector<T>> get_M() const { return M; }
void operator =(const matrix& m){ (*this).M = m.get_M(); }

void operator +=(const matrix& m){ 
  const size_t nr = this->rows(), nc = this->cols();  
  for(size_t i=0; i<nr; ++i){ 
    for (size_t j = 0; j < nc; j++)
    {
          M[i][j] += m[i][j];
    }  
  }
  return;
}

};


template<class T>
std::vector<T> operator *(const matrix<T>& A, const std::vector<T>& v)
{
        const int n = A.rows(), m = v.size();
        std::vector<T> w(n, 0.0);
        for (size_t i = 0; i < n; i++) {
                for (size_t k = 0; k < m; k++) {
                        w[i] += A[i][k]*v[k];
                }
        }
        return w;
}

template<class T, class P>
std::vector<T> a_plus_lambda_b(const std::vector<T>& a, const P& lambda, const std::vector<T>& b)
{
        const int n = a.size();
        std::vector<T> c(n, 0.0);
        for (size_t i = 0; i < n; i++) {
                c[i] = a[i] + lambda*b[i];
        }
        return c;

}

template<class T>
std::vector<T> operator -(const std::vector<T>& a, const std::vector<T>& b)
{
        return a_plus_lambda_b(a, -1.0, b);
}

template<class T>
std::vector<T> operator +(const std::vector<T>& a, const std::vector<T>& b)
{
        return a_plus_lambda_b(a, +1.0, b);
}

template<class T>
std::vector<T> operator -(const std::vector<T>& a)
{
        std::vector<T> foo(a.size(), 0.0);
        return foo-a; // \vec{0} - \vec{a}
}

template<class T>
std::vector<T> operator *(const T& lambda, const std::vector<T>& v)
{
        std::vector<T> foo(v.size(), 0.0);
        return a_plus_lambda_b(foo, lambda, v);// \vec{0} + lambda*\vec{v}
}


template<class T>
T dot(const std::vector<T>& a, const std::vector<T>& b)
{
        T x=0.0;
        const int n = a.size();
        for (size_t i = 0; i < n; i++) {
                x += std::conj(a[i]) * b[i];
        }
        return x;
}

template<class T>
T operator*(const std::vector<T>& a, const std::vector<T>& b)
{
        return dot(a,b);
}


template<class Float, class T>
Float norm(const std::vector<T>& v)
{
        return (Float) sqrt(dot(v,v).real());
}

// division by a scalar
template<class T> matrix<T> operator /(const matrix<T>& M, const T& lambda)
{
  const size_t nr = M.rows(), nc = M.cols();
  matrix<T> Q(nr,nc);

  for (size_t i = 0; i < nr; i++)
  {
    for (size_t j = 0; j < nc; j++)
    {
        Q[i][j] = M[i][j]/lambda;
    }
  }
  return Q;        
} 