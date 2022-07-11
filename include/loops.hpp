/**
 * @file loops.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-07-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <iostream>
#include <vector>

#include "gaugeconfig.hh"
#include "geometry.hh"

namespace links {
  template <class T> using nd_max_arr = typename spacetime_lattice::nd_max_arr<T>;

  nd_max_arr<int> operator+(const nd_max_arr<int> &a, const nd_max_arr<int> &b) {
    nd_max_arr<int> c;
    for (size_t i = 0; i < spacetime_lattice::nd_max; i++) {
      c[i] = a[i] + b[i];
    }
    return c;
  }

  /**
   * @brief path from one point (the origin) to another according to a sequence of steps.
   * Example: d=2; (0,0)->(0,1)->(0,2)->(-1,2)
   *
   */
  class path {
  private:
    std::vector<nd_max_arr<int>> steps; // unit vectors of each step
    std::vector<size_t> dirs; // directions
    std::vector<bool> pos; // whether the unit vectors are along the positive direction

  public:
    path() {}
    path(const std::vector<nd_max_arr<int>> &_steps,
         const std::vector<size_t> &_dirs,
         const std::vector<bool> &_positive) {
      (*this).steps = _steps;
      (*this).dirs = _dirs;
      (*this).pos = _positive;
    }
    ~path() {}

    bool is_closed() const {
      const int n = (*this).steps.size();
      nd_max_arr<int> res;
      res.fill(0);
      for (size_t i = 0; i < n; i++) {
        res = res + steps[i];
      }
      nd_max_arr<int> zero_arr;
      zero_arr.fill(0);
      return (res == zero_arr);
    }

    std::vector<nd_max_arr<int>> get_steps() const { return (*this).steps; }
    std::vector<size_t> get_dirs() const { return (*this).dirs; }
    std::vector<bool> get_pos() const { return (*this).pos; }

    void append(const path &p) {
      std::vector<nd_max_arr<int>> steps1 = p.get_steps();
      std::vector<size_t> dirs1 = p.get_dirs();
      std::vector<bool> pos1 = p.get_pos();

      (*this).steps.insert((*this).steps.end(), steps1.begin(), steps1.end());
      (*this).dirs.insert((*this).dirs.end(), dirs1.begin(), dirs1.end());
      (*this).pos.insert((*this).pos.end(), pos1.begin(), pos1.end());
      return;
    }

    /**
     * @brief product of links along the path
     *
     * @tparam Group gauge group U belongs to
     * @param p path
     * @param U gauge configuration
     * @param x starting point
     * @param apply_P true when applying spatial parity operator
     * @return Group
     */
    template <class Group>
    Group to_gauge(const gaugeconfig<Group> &U,
                   const nd_max_arr<int> &x,
                   const bool &apply_P = false) const {
      const std::vector<nd_max_arr<int>> steps = this->get_steps();
      const std::vector<size_t> dirs = this->get_dirs();
      const std::vector<bool> pos = this->get_pos();

      size_t n_steps = steps.size();
      nd_max_arr<int> y = x;
      if (apply_P) {
        for (size_t i = 1; i < spacetime_lattice::nd_max; i++) {
          y[i] *= -1; // (y0, yi) --> (y0, -yi)
        }
      }

      Group L;
      for (size_t i = 0; i < n_steps; i++) {
        const Group Ui = U(y, dirs[i], pos[i] ^ apply_P);

        L = L * Ui;

        y[dirs[i]] += std::pow(-1, apply_P && (dirs[i] != 0)) *
                      steps[i][dirs[i]]; // updating 'y' after one step
      }
      return L;
    }
  };

  /**
   * @brief On a d-dimensional lattice, find all closed path of a given lenght
   *
   */
  class closed_paths {
  private:
    size_t d = spacetime_lattice::nd_max; // number of dimensions
    size_t first_dim = 0; // index of 1st dimension
    size_t length; // number of steps
    std::vector<path> P; // vector of closed paths

    void find_all(path &pt, const size_t &l) {
      if (l == 0) {
        std::vector<bool> pos = pt.get_pos();
        // considering only closed loops
        bool cond = pt.is_closed();
        // avoid loops equivalent up to translations
        cond = cond && (pos[0] == true) && (pos.back() == false);
        if (cond) {
          P.push_back(pt);
        }
      } else {
        for (size_t i = (*this).first_dim; i < (*this).d; i++) {
          for (int s = -1; s < 2; s += 2) {
            nd_max_arr<int> step;
            step.fill(0);
            step[i] = s;

            // we don't want loops were after one step, the following is the same but
            // reversed.
            nd_max_arr<int> zero_step;
            zero_step.fill(0);
            nd_max_arr<int> prev_step;
            prev_step.fill(0);
            if (pt.get_steps().size() > 0) {
              prev_step = pt.get_steps().back();
            }
            if ((step + prev_step) == zero_step) {
              continue;
            }

            const path pi({step}, {i}, {(s == 1)}); // elementary path of just one step

            path pt2 = pt; // copy of the present path
            pt2.append(pi);
            this->find_all(pt2, l - 1);
          }
        }
      }
      return;
    }

    void find_all(const size_t &l) {
      path pt;
      this->find_all(pt, l);
    }

  public:
    closed_paths(const size_t &_d, const size_t &_first_dim, const size_t &_length) {
      d = _d;
      first_dim = _first_dim;
      if ((first_dim + d) > spacetime_lattice::nd_max) {
        std::cerr << "Error: invalid number of dimensions and/or 1st dimension index. "
                     "Aborting.\n";
        std::abort();
      }
      length = _length;
      this->find_all(length);
    }

    ~closed_paths() {}

    size_t size() const { return P.size(); }

    std::vector<path> get_paths() const { return P; }

    void print() const {
      const size_t n = this->size();
      for (size_t i = 0; i < n; i++) {
        std::cout << "#path: " << i << "\n";
        std::vector<nd_max_arr<int>> steps = P[i].get_steps();
        for (size_t s = 0; s < steps.size(); s++) {
          const size_t ns = steps[s].size();
          for (size_t j = 0; j < ns; j++) {
            std::cout << steps[s][j] << " ";
          }
          std::cout << "\n";
        }
      }
    }
  };

  /**
   * @brief class for generating all closed loops of a given lenght
   *
   * @tparam Group
   */
  template <class Group> class links_loops {
  private:
    std::vector<Group> loops;
    gaugeconfig<Group> U;
    nd_max_arr<int> x;
    size_t d; // number of dimensions
    size_t first_dim; // index of 1st dimension
    size_t length; // lenght of loops

  public:
    links_loops(const gaugeconfig<Group> &_U,
                const nd_max_arr<int> _x,
                const size_t &_first_dim,
                const size_t &_length) {
      if (_length < 4) {
        std::cerr << "Error: cannot create closed loops of gauge links with length < 4. "
                     "Aborting.\n";
        std::abort();
      }
      length = _length;
      U = _U;
      x = _x;
      first_dim = _first_dim;
      d = U.getndims();
    };
    ~links_loops() {}

    Group from_path(const path &p) { return p.to_gauge<Group>(U, x); }

    void find_all(const size_t &li) {
      closed_paths CP(d, first_dim, li);
      std::vector<path> P = CP.get_paths();
      for (size_t i = 0; i < P.size(); i++) {
        loops.push_back(this->from_path(P[i]));
      }
    }

    void find_all_shorter() {
      for (size_t i = 1; i <= length; i++) { /* note: i==0 is meaningless */
        this->find_all(i);
      }
      return;
    }

    std::vector<Group> get_loops() const { return (*this).loops; }
  };

} // namespace links