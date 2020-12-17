#pragma once

#include"su2.hh"
#include<vector>

using std::vector;

class gaugeconfig {
public:
  using value_type = su2;

  gaugeconfig(const size_t Ls, const size_t Lt, const double beta=0) : 
    Lx(Ls), Ly(Ls), Lz(Ls), Lt(Lt), volume(Lx*Ly*Lz*Lt), beta(beta), ndims(4) {
    data.resize(volume*4);
  }
  gaugeconfig(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, const size_t ndims=4, const double beta=0) : 
    Lx(Lx), Ly(Ly), Lz(Lz), Lt(Lt), volume(Lx*Ly*Lz*Lt), beta(beta), ndims(ndims) {
    data.resize(volume*4);
  }
  gaugeconfig(const gaugeconfig &U) :
    Lx(U.getLx()), Ly(U.getLy), Lz(U.getLz), Lt(U.getLt()), volume(U.getVolume()), beta(U.getBeta()), ndims(U.getndims()) {
    data.resize(volume*4);
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
    return(volume*4);
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

  void operator=(const gaugeconfig &U) {
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
    size_t _mu = (mu + 4) % 4;
    return( (((y0*Lx + y1)*Ly + y2)*Lz + y3)*4 + _mu );
  }
};

gaugeconfig hotstart(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, 
                     const int seed, const double _delta, sonst size_t ndims);
gaugeconfig coldstart(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, const size_t ndims);
