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

void gaugeconfig::load(std::string const &path) {
  std::cout << "## Reading config from file " << path << std::endl;
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  ifs.read(reinterpret_cast<char *>(data.data()), storage_size());
  return;
}

void gaugeconfig::loadEigen(std::string const &path) {
  std::cout << "## Reading config from file " << path << " in Eigen format" << std::endl;
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  double X[8];
  for(size_t t = 0; t < Lt; t++) {
    for(size_t x = 0; x < Ls; x++) {
      for(size_t y = 0; y < Ls; y++) {
        for(size_t z = 0; z < Ls; z++) {
          //size_t coord[4] = {t, x, y, z};
          for(size_t mu = 0; mu < 4; mu++) {
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

gaugeconfig coldstart(size_t Ls, size_t Lt) {

  gaugeconfig config(Ls, Lt);
#pragma omp parallel for
  for(size_t i = 0; i < config.getSize(); i++) {
    config[i] = su2(1., 0.);
  }
  return(config);
}
gaugeconfig hotstart(size_t Ls, size_t Lt, 
                     const int seed, const double _delta) {

  gaugeconfig config(Ls, Lt);
  double delta = _delta;
  if(delta < 0.) delta = 0;
  if(delta > 1.) delta = 1.;
  std::mt19937 engine(seed);

  for(size_t i = 0; i < config.getSize(); i++) {
    random_su2(config[i], engine, delta);
  }
  return(config);
}

