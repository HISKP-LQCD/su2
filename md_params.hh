#pragma once

class md_params {

public:
  md_params(size_t _nsteps, double _tau) : nsteps(_nsteps), tau(_tau) {}
  size_t getnsteps() const {
    return nsteps;
  }
  double gettau() const {
    return tau;
  }
  void setnsteps(const size_t _nsteps) {
    nsteps = _nsteps;
  }
  void settau(const double _tau) {
    tau = _tau;
  }

private:
  size_t nsteps;
  double tau;
};

