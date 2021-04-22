#pragma once

template<class URNG> class md_params {

public:
  md_params(size_t _nsteps, double _tau, size_t n = 1000000000) : 
    nsteps(_nsteps), n_prec(n), kmax(5), exponent(16), tau(_tau), H(0.),
    deltaH(0.), deltadeltaH(0.), deltadeltaU(0.),
    gamma(2.0), omflambda(0.1938), temperedwidth(0.06),
    revtest(false), accept(true), acceptreject(true) {}
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
  void setgamma(const double g) {
    gamma = g;
  }
  double getgamma() const {
    return gamma;
  }
  void setaccept(const bool a) {
    accept = a;
  }
  bool getaccept() const {
    return accept;
  }
  void disableacceptreject() {
    acceptreject = false;
  }
  void enableacceptreject() {
    acceptreject = true;
  }
  bool getacceptreject() const {
    return acceptreject;
  }
  void setn_prec(size_t n) {
    n_prec = n;
  }
  size_t getn_prec() const {
    return n_prec;
  }
  void setkmax(size_t k) {
    kmax = k;
    return;
  }
  size_t getkmax() const {
    return kmax;
  }
  void setexponent(size_t k) {
    exponent = k;
  }
  size_t getexponent() const {
    return exponent;
  }
  void setengine(URNG * _engine) {
    engine = _engine;
  }
  URNG * getengine() {
    return(engine);
  }
  double gettemperedwidth() {
    return temperedwidth;
  }
private:
  size_t nsteps;
  size_t n_prec;
  size_t kmax;
  size_t exponent;
  double tau;
  double H, deltaH, deltadeltaH, deltadeltaU;
  double gamma;
  double omflambda;
  double temperedwidth;
  bool revtest, accept, acceptreject;
  URNG * engine;
};

