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

