// staggered Dirac operator

#include "include/gaugeconfig.hh"
#include <fstream>
#include <iostream>
#include <string>

double get_available_RAM_Gb() {
  // ACHTUNG: it is assumed you are running on a Linux machine!
  std::ifstream meminfo("/proc/meminfo");

  if (!meminfo.is_open()) {
    std::cerr << "Error opening /proc/meminfo" << std::endl;
    std::abort();
  }

  std::string line;
  double g = 0;

  // first check from the physical RAM
  while (getline(meminfo, line)) {
    if (line.find("MemTotal:") != std::string::npos) {
      // Found the line containing total memory information
      const std::string totalMemoryString = line.substr(10);
      const double k = std::stoi(totalMemoryString); // RAM in kb
      g += k * (1e-6);
      break;
    }
  }

  // now check for the possible SWAP partition
  while (getline(meminfo, line)) {
    if (line.find("SwapTotal:") != std::string::npos) {
      // Found the line containing total memory information
      const std::string swapMemoryString = line.substr(10);
      const double k = std::stoi(swapMemoryString); // RAM in kb
      g += k * (1e-6);
      break;
    }
  }

  meminfo.close();

  return g;
}

void check_RAM_allocation(const size_t &size_in_Gb) {
  if (size_in_Gb > get_available_RAM_Gb()) {
    std::cerr << "## Available RAM " << get_available_RAM_Gb() << " Gb\n";
    std::cerr << "## Requested RAM " << size_in_Gb << " Gb\n";
    std::cerr << "## Error from " << __func__ << ":\n";
    std::cerr << "## Aborting. \n";
    std::abort();
  }
}

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
   * Ref: eq. 6.47 of Degrand's book
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

    int i = 0;
    std::vector<int> eta = {1, 1, 1, 1};
    for (size_t x0 = 0; x0 < U.getLt(); x0++) {
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            std::vector<size_t> x = {x0, x1, x2, x3};
            for (size_t mu = 0; mu < U.getndims(); mu++) {
              const int eta_mu =
                std::accumulate(eta.begin(), eta.begin() + mu, 1, std::multiplies<int>());
              // hopping terms
              col_idx[3 * i] = (i + 1) % N;
              D[3 * i] = eta_mu * accum(U(x, mu));

              x[mu] = (x[mu] - 1 + N) % N; // x - mu
              col_idx[3 * i + 1] = (i - 1 + N) % N;
              D[3 * i + 1] = -eta_mu * accum(U(x, mu));
              x[mu] = (x[mu] + 1) % N; // x

              // mass term
              col_idx[3 * i + 2] = i;
              D[3 * i + 2] = +m * accum(U(x, mu));

              i++;
            }
            eta[3] *= -1;
          }
          eta[2] *= -1;
        }
        eta[1] *= -1;
      }
      eta[0] *= -1;
    }
    return std::pair(D, col_idx);
  }

} // namespace staggered