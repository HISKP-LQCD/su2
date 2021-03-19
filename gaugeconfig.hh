#pragma once

#include"su2.hh"
#include"u1.hh"
#include"random_element.hh"
#include<random>
#include<vector>
#include<cmath>
#include<fstream>
#include<complex>
#include<iostream>
#include<cassert>

using std::vector;

template<class T> class gaugeconfig {
public:
  using value_type = T;

  gaugeconfig(const size_t Ls, const size_t Lt, const double beta=0) : 
    Ls(Ls), Lt(Lt), volume(Ls*Ls*Ls*Lt), beta(beta) {
    data.resize(volume*4);
  }
  gaugeconfig(const gaugeconfig &U) :
    Ls(U.getLs()), Lt(U.getLt()), volume(U.getVolume()), beta(U.getBeta()) {
    data.resize(volume*4);
#pragma omp parallel for
    for(size_t i = 0; i < getSize(); i++) {
      data[i] = U[i];
    }
  }

  size_t storage_size() const { return data.size() * sizeof(value_type); };
  size_t getLs() const {
    return(Ls);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getVolume() const {
    return(volume);
  }
  size_t getSize() const {
    return(volume*4);
  }
  double getBeta() const {
    return beta;
  }
  void setBeta(const double _beta){
    beta = _beta;
  }
  int getNc() {
    return(data[0].N_c);
  }
  void restoreSU() {
#pragma omp parallel for
    for(size_t i = 0; i < getSize(); i++) {
      data[i].restoreSU();
    }
  }

  void operator=(const gaugeconfig &U) {
    Ls = U.getLs();
    Lt = U.getLt();
    volume = U.getVolume();
    data.resize(U.getSize());
#pragma omp parallel for
    for(size_t i = 0; i < U.getSize(); i++) {
      data[i] = U[i];
    }
  }

  value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  const value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) const {
    return data[ getIndex(t, x, y, z, mu) ];
  }

  value_type &operator()(std::vector<size_t> const &coords, size_t const mu) {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  const value_type &operator()(std::vector<size_t> const &coords, size_t const mu) const {
    return data[ getIndex(coords[0], coords[1], coords[2], coords[3], mu) ];
  }

  value_type &operator[](size_t const index) {
    return data[ index ];
  }

  const value_type &operator[](size_t const index) const {
    return data[index];
  }

  void save(std::string const &path) const;
  int load(std::string const &path);

private:
  size_t Ls, Lt, volume;
  double beta;

  vector<value_type> data;

  size_t getIndex(const size_t t, const size_t x, const size_t y, const size_t z, const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Ls) % Ls;
    size_t y2 = (y + Ls) % Ls;
    size_t y3 = (z + Ls) % Ls;
    size_t _mu = (mu + 4) % 4;
    return( (((y0*Ls + y1)*Ls + y2)*Ls + y3)*4 + _mu );
  }
};

template<class T> void gaugeconfig<T>::save(std::string const &path) const {
  std::ofstream ofs(path, std::ios::out | std::ios::binary);
  ofs.write(reinterpret_cast<char const *>(data.data()), storage_size());
  return;
}

template<class T> int gaugeconfig<T>::load(std::string const &path) {
  std::cout << "## Reading config from file " << path << std::endl;
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  if(ifs) {
    ifs.read(reinterpret_cast<char *>(data.data()), storage_size());
    return 0;
  }
  else
    std::cerr << "Error: could not read file from " << path << std::endl;
  return 1;
}


template<class T> void coldstart(gaugeconfig<T> &config) {

#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = T(1., 0.);
  }
}

template<class T> void coldstart(gaugeconfig<_u1> &config) {

#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = _u1(0.);
  }
}

template<class T> void  hotstart(gaugeconfig<T> & config,
                                 const int seed, const double _delta) {

  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);

  for(size_t i = 0; i < config.getSize(); i++) {
    random_element(config[i], engine, delta);
  }
}


