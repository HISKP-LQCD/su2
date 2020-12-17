#include"su2.hh"
#include"random_su2.hh"
#include"gaugeconfig.hh"
#include<random>
#include<vector>
#include<cmath>
#include<fstream>
#include<complex>
#include<iostream>
#include<cassert>

void gaugeconfig::save(std::string const &path) const {
  std::ofstream ofs(path, std::ios::out | std::ios::binary);
  ofs.write(reinterpret_cast<char const *>(data.data()), storage_size());
  return;
}

int gaugeconfig::load(std::string const &path) {
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

void gaugeconfig::loadEigen(std::string const &path) {
  std::cout << "## Reading config from file " << path << " in Eigen format" << std::endl;
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  double X[8];
  for(size_t t = 0; t < Lt; t++) {
    for(size_t x = 0; x < Lx; x++) {
      for(size_t y = 0; y < Ly; y++) {
        for(size_t z = 0; z < Lz; z++) {
          //size_t coord[4] = {t, x, y, z};
          for(size_t mu = 0; mu < ndims; mu++) {
            if(ifs.good()) {
              ifs.read(reinterpret_cast<char *>(X), 8*sizeof(double));
              su2 U(std::complex<double>(X[0], X[1]), std::complex<double>(X[4], X[5]));
              data[ getIndex(t, x, y, z, mu) ] = U;
              //std::cout << X[0] << " " << coord[mu] << std::endl;
              //std::cout << U.det() << " " << trace(U) <<  " " << ifs.gcount() << std::endl;
            }
            else {
              std::cerr << "Could not read all data from file " << path << std::endl;
            }
          }
        }
      }
    }
  }
  return;
}

gaugeconfig coldstart(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, const size_t ndims) {

  gaugeconfig config(Lx, Ly, Lz, Lt, ndims);
#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = su2(1., 0.);
  }
  return(config);
}

gaugeconfig hotstart(const size_t Lx, const size_t Ly, const size_t Lz, const size_t Lt, 
                     const int seed, const double _delta, const size_t ndims) {

  gaugeconfig config(Lx, Ly, Lz, Lt, ndims);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);

  for(size_t i = 0; i < config.getSize(); i++) {
    random_su2(config[i], engine, delta);
  }
  return(config);
}

