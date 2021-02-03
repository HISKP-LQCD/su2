#pragma once

#include"su2.hh"
#include"genzsu2.hh"
#include"random_su2.hh"
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

  gaugeconfig(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, const size_t ndims=4, const double beta=0) : 
    Lx(Lx), Ly(Ly), Lz(Lz), Lt(Lt), volume(Lx*Ly*Lz*Lt), beta(beta), ndims(ndims) {
    data.resize(volume*ndims);
  }
  gaugeconfig(const gaugeconfig<T> &U) :
    Lx(U.getLx()), Ly(U.getLy()), Lz(U.getLz()), Lt(U.getLt()), volume(U.getVolume()), beta(U.getBeta()), ndims(U.getndims()) {
    data.resize(volume*ndims);
#pragma omp parallel for
    for(size_t i = 0; i < getSize(); i++) {
      data[i] = U[i];
    }
  }

  size_t storage_size() const { return data.size() * sizeof(value_type); };
  size_t getLx() const {
    return(Lx);
  }
  size_t getLy() const {
    return(Ly);
  }
  size_t getLz() const {
    return(Lz);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getndims() const {
    return(ndims);
  }
  size_t getVolume() const {
    return(volume);
  }
  size_t getSize() const {
    return(volume*ndims);
  }
  double getBeta() const {
    return beta;
  }
  void setBeta(const double _beta){
    beta = _beta;
  }

  void restoreSU() {
#pragma omp parallel for
    for(size_t i = 0; i < getSize(); i++) {
      data[i].restoreSU();
    }
  }

  void operator=(const gaugeconfig<T> &U) {
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
  void loadEigen(std::string const &path);

private:
  size_t Lx, Ly, Lz, Lt, volume, ndims;
  double beta;

  vector<value_type> data;

  size_t getIndex(const size_t t, const size_t x, const size_t y, const size_t z, const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    size_t _mu = (mu + ndims) % ndims;
    return( (((y0*Lx + y1)*Ly + y2)*Lz + y3)*ndims + _mu );
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

gaugeconfig<Gsu2> coldstart(const size_t Lx, const size_t Ly,
                            const size_t Lz, const size_t Lt,
                            const size_t ndims, const size_t m=10) {
  gaugeconfig<Gsu2> config(Lx, Ly, Lz, Lt, ndims);
#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = Gsu2(m);
  }
  return(config);
}

gaugeconfig<su2> coldstart(const size_t Lx, const size_t Ly,
                           const size_t Lz, const size_t Lt,
                           const size_t ndims) {
  
  gaugeconfig<su2> config(Lx, Ly, Lz, Lt, ndims);
#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = su2(1., 0.);
  }
  return(config);
}

gaugeconfig<Gsu2> hotstart(const size_t Lx, const size_t Ly,
                           const size_t Lz, const size_t Lt, 
                           const int seed, const size_t m,
                           const double _delta, const size_t ndims) {
  gaugeconfig<Gsu2> config(Lx, Ly, Lz, Lt, ndims);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);

  for(size_t i = 0; i < config.getSize(); i++) {
    random_su2(config[i], engine, m, delta);
  }
  return(config);
}

gaugeconfig<Osu2> Ohotstart(const size_t Lx, const size_t Ly,
                            const size_t Lz, const size_t Lt, 
                            const int seed, const size_t m,
                            const double _delta, const size_t ndims) {
  gaugeconfig<Osu2> config(Lx, Ly, Lz, Lt, ndims);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);

  for(size_t i = 0; i < config.getSize(); i++) {
    random_su2(config[i], engine, m, delta);
  }
  return(config);
}


gaugeconfig<su2> hotstart(const size_t Lx, const size_t Ly,
                          const size_t Lz, const size_t Lt, 
                          const int seed, const double _delta,
                          const size_t ndims) {
  
  gaugeconfig<su2> config(Lx, Ly, Lz, Lt, ndims);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);
  
  for(size_t i = 0; i < config.getSize(); i++) {
    random_su2(config[i], engine, delta);
  }
  return(config);
}

