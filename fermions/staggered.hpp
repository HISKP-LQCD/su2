// staggered Dirac operator

#include <fstream>
#include <iostream>
#include <string>

#include "fermions/memory.hpp"
#include "include/gaugeconfig.hh"

namespace staggered {

  /**
   * Dirac operator components as numbers. The index is the checkerboard for
   * spacetime+color.
   *
   * The elements are stored as for a sparse matrix.
   * D = pointer to values
   * col_idx = pointer to column components
   * The row indices are implicit. They can be thought as 3*N triplets:
   * {0,0,0, 1,1,1, 2,2,2, ...}.
   * Like so, the non-vanishing matrix elements are:
   * (row_idx[j], col_idx[j]) : D[j] for all "0 <= j < 3*N"
   *
   * Ref: eq. 6.47 of Degrand's book:
   * https://www.worldscientific.com/worldscibooks/10.1142/6065?srsltid=AfmBOopLT5inz0mksWaaeXm_5YTYw8dI3gLXDYNTZfVZw39SDivh9opG#t=aboutBook
   */
  template <class Float, class Group>
  std::pair<Float *, size_t *> *dirac_op_ptr(const gaugeconfig<Group> &U,
                                             const double &m);

  template <class Float>
  std::pair<Float *, size_t *> dirac_op_ptr(const gaugeconfig<u1> &U, const double &m) {
    /**
     */
    typedef typename accum_type<u1>::type accum;
    const size_t N_c = 1; // U(1) theory
    const size_t N = U.getSize() * N_c;

    // checking that we have enough space in the memory
    const double mem_elements = 3 * N * sizeof(Float) / 1e+9; // elements of D
    const double mem_col_idx = 3 * N * sizeof(size_t) / 1e+9; // column indices
    check_RAM_allocation(mem_elements + mem_col_idx);

    Float *D;
    D = (Float *)(calloc(3 * N, sizeof(Float)));
    size_t *col_idx; // column indices
    col_idx = (size_t *)(calloc(3 * N, sizeof(size_t)));

    geometry Geom = U.get_geometry(); // geometry of the lattice
    std::vector<size_t> L = Geom.get_L(); // lattice sizes
    size_t n_dims = Geom.get_n_dims(); // number of dimensions
    size_t N_pts = Geom.get_N_pts(); // number of points

    // std::vector<int> eta_exp(n_dims, 0); // exponents of \eta_\mu

    for (size_t i = 0; i < N_pts; i++) {
      std::vector<int> x = spacetime_lattice::index_to_x<int>(i, L);
      const size_t i_x = n_dims * i;
      int eta_exp_mu = 0;
      for (size_t mu = 0; mu < n_dims; mu++) {
        size_t i_g = i_x + mu; // global index
        // (-1)^{x_0+x_1+...+x_{\mu-1}}
        eta_exp_mu += (x[mu] % 2);
        const int eta_mu = std::pow(-1.0, eta_exp_mu);

        std::cout << x[0] << " " << x[1] << " " << x[2] << " | n_dims=" << n_dims << "\n";
        std::cout << eta_exp_mu << " " << eta_mu << "\n";

        // hopping terms
        col_idx[3 * i_g] = (i_g + 1) % N;
        D[3 * i_g] = eta_mu * accum(U(x, mu));

        x[mu] = (x[mu] - 1 + N) % N; // x - mu
        col_idx[3 * i_g + 1] = (i_g - 1 + N) % N;
        D[3 * i_g + 1] = -eta_mu * accum(U(x, mu));
        x[mu] = (x[mu] + 1) % N; // x

        // mass term
        col_idx[3 * i_g + 2] = i_g;
        D[3 * i_g + 2] = +m * accum(U(x, mu));

        // i_g++;
      }
    }
    return std::pair(D, col_idx);
  }

} // namespace staggered