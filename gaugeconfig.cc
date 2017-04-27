#include"su2.hh"
#include"random_su2.hh"
#include"gaugeconfig.hh"
#include<random>
#include<vector>
#include<cmath>
#include <fstream>

void gaugeconfig::save(std::string const &path) const {
  std::ofstream os(path, std::ios::out | std::ios::binary);
  os.write(reinterpret_cast<char const *>(data.data()), storage_size());
}

void gaugeconfig::load(std::string const &path) {
  std::ifstream ifs(path, std::ios::out | std::ios::binary);
  ifs.read(reinterpret_cast<char *>(data.data()), storage_size());
}


gaugeconfig coldstart(size_t Ls, size_t Lt) {

  gaugeconfig config(Ls, Lt);
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
