#pragma once

class md_params {

public:
  md_params(size_t _nsteps, double _tau) : 
    nsteps(_nsteps), tau(_tau), H(0.),
    deltaH(0.), deltadeltaH(0.), deltadeltaU(0.),
    revtest(false), accept(true) {}
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
  void enablerevtest() {
    revtest = true;
  }
  void disablerevtest() {
    revtest = false;
  }
  bool getrevtest() const {
    return revtest;
  }
  void setdeltaH(const double h) {
    deltaH = h;
  }
  double getdeltaH() const {
    return deltaH;
  }
  void setdeltadeltaH(const double h) {
    deltadeltaH = h;
  }
  double getdeltadeltaH() const {
    return deltadeltaH;
  }
  void setdeltadeltaU(const double h) {
    deltadeltaU = h;
  }
  double getdeltadeltaU() const {
    return deltadeltaU;
  }
  void setH(const double h) {
    H = h;
  }
  double getH() const {
    return H;
  }
  void setaccept(const bool a) {
    accept = a;
  }
  bool getaccept() const {
    return accept;
  }

private:
  size_t nsteps;
  double tau;
  double H, deltaH, deltadeltaH, deltadeltaU;
  bool revtest, accept;
};

