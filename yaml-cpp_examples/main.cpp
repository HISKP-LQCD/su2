// a.cpp

//#include "yce.h"
#include <iostream>
#include "yaml.h"


int main(int argc, char const *argv[]) {

  const YAML::Node config = YAML::LoadFile("./config.yaml");
  std::cout << "p1 : " << config["p1"].as<double>() << "\np2 : " << config["p2"].as<std::string>() << "\n";

  return 0;
}
