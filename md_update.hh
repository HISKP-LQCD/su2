#pragma once

class md_params {
public:
  size_t getnsteps() const {
    return nsteps;
  }
  double gettau() const {
    return tau;
  }

private:
  size_t nsteps;
  double tau;
};
