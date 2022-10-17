// LA.hpp
// containes for linear algebra: vectors and matrices

#pragma once

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

namespace LA {

  // vector for Linear Algebra operations
  template <class Float, class T> class LAvector {
  private:
    std::vector<T> v;

  public:
    LAvector() {}
    ~LAvector() {}

    LAvector(const std::vector<T> &_v) { v = _v; }
    LAvector(const size_t &n) { this->resize(n); }
    LAvector(const size_t &n, const T &val) {
      std::vector<T> w(n, val);
      v = w;
    }

    // storing enough space for the LAvector elements
    void resize(const size_t &n) { v.resize(n); }
    size_t size() const { return v.size(); }

    T &operator[](const int &i) { return v[i]; }
    T operator[](const int &i) const { return v[i]; }

    std::vector<T> get_v() const { return v; }

    void operator=(const LAvector &_V) { (*this).v = _V.get_v(); }
    void operator=(const std::vector<T> &_v) { (*this).v = _v; }

    // dot procuct \sum_{i} v^{\dagger}_i * w_i
    T dot(const LAvector &w) const {
      T x = 0.0;
      const int n = v.size();
      for (size_t i = 0; i < n; i++) {
        x += std::conj(v[i]) * w[i];
      }
      return x;
    }
    Float norm() { return (Float)sqrt((this->dot(v)).real()); }

    // real dot  procuct \sum_{i} v_i * w_i
    T real_dot(const LAvector &w) const {
      T x = 0.0;
      const int n = v.size();
      for (size_t i = 0; i < n; i++) {
        x += v[i] * w[i];
      }
      return x;
    }


  };

  template <class T> class LAmatrix {
  private:
    std::vector<std::vector<T>> M;

  public:
    LAmatrix() {}
    ~LAmatrix() {}

    LAmatrix(const std::vector<std::vector<T>> &m) {
      const int nr = m.size();
      const int nc = m[0].size();

      M.resize(nr); // saving the space for enough rows
      for (int i = 0; i < nr; ++i) {
        M[i].resize(nc); // saving the space for enough column entries
        for (int j = 0; j < nc; ++j) {
          M[i][j] = m[i][j];
        }
      }
    }

    LAmatrix(const size_t &nr, const size_t &nc) { this->resize(nr, nc); };

    // storing enough space for the LAmatrix elements
    void resize(const size_t &nr, const size_t &nc) {
      M.resize(nr);
      for (int i = 0; i < nr; ++i) {
        M[i].resize(nc);
      }
    }

    size_t rows() const { return M.size(); }
    size_t cols() const { return M[0].size(); }

    void add_row(const std::vector<T> &v) { M.push_back(v); }

    std::vector<T> &operator[](const int &i) { return M[i]; }
    std::vector<T> operator[](const int &i) const { return M[i]; }

    std::vector<std::vector<T>> get_M() const { return M; }
    void operator=(const LAmatrix &m) { (*this).M = m.get_M(); }

    void operator+=(const LAmatrix &m) {
      const size_t nr = this->rows(), nc = this->cols();
      for (size_t i = 0; i < nr; ++i) {
        for (size_t j = 0; j < nc; j++) {
          M[i][j] += m[i][j];
        }
      }
      return;
    }
  };

  template <class Float, class T>
  LAvector<Float, T> operator*(const LAmatrix<T> &A, const LAvector<Float, T> &v) {
    const int n = A.rows(), m = v.size();
    LAvector<Float, T> w(n, 0.0);
    for (size_t i = 0; i < n; i++) {
      for (size_t k = 0; k < m; k++) {
        w[i] += A[i][k] * v[k];
      }
    }
    return w;
  }

  template <class Float, class T, class P>
  LAvector<Float, T> a_plus_lambda_b(const LAvector<Float, T> &a,
                                     const P &lambda,
                                     const LAvector<Float, T> &b) {
    const int n = a.size();
    LAvector<Float, T> c(n, 0.0);
    for (size_t i = 0; i < n; i++) {
      c[i] = a[i] + lambda * b[i];
    }
    return c;
  }

  template <class Float, class T>
  LAvector<Float, T> operator-(const LAvector<Float, T> &a, const LAvector<Float, T> &b) {
    return a_plus_lambda_b(a, -1.0, b);
  }

  template <class Float, class T>
  LAvector<Float, T> operator+(const LAvector<Float, T> &a, const LAvector<Float, T> &b) {
    return a_plus_lambda_b(a, +1.0, b);
  }

  template <class Float, class T>
  LAvector<Float, T> operator-(const LAvector<Float, T> &a) {
    LAvector<Float, T> foo(a.size(), 0.0);
    return foo - a; // \vec{0} - \vec{a}
  }

  template <class Float, class T>
  LAvector<Float, T> operator*(const T &lambda, const LAvector<Float, T> &v) {
    LAvector<Float, T> foo(v.size(), 0.0);
    return a_plus_lambda_b(foo, lambda, v); // \vec{0} + lambda*\vec{v}
  }

  //   template <class Float, class T>
  //   T dot(const LAvector<Float, T> &a, const LAvector<Float, T> &b) {
  //     T x = 0.0;
  //     const int n = a.size();
  //     for (size_t i = 0; i < n; i++) {
  //       x += std::conj(a[i]) * b[i];
  //     }
  //     return x;
  //   }

  template <class Float, class T>
  T operator*(const LAvector<Float, T> &a, const LAvector<Float, T> &b) {
    return a.dot(b);
  }

  // division by a scalar
  template <class T> LAmatrix<T> operator/(const LAmatrix<T> &M, const T &lambda) {
    const size_t nr = M.rows(), nc = M.cols();
    LAmatrix<T> Q(nr, nc);

    for (size_t i = 0; i < nr; i++) {
      for (size_t j = 0; j < nc; j++) {
        Q[i][j] = M[i][j] / lambda;
      }
    }
    return Q;
  }

} // namespace LA
