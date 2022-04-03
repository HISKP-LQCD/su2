// a.cpp

//#include "yce.h"
#include <iostream>
#include "yaml.h"


int main(int argc, char const *argv[]) {

  const YAML::Node config = YAML::LoadFile("./config.yaml");
  std::cout << "p1 : " << config["p1"].as<double>() << "\np2 : " << config["p2"].as<std::string>() << "\n";

  std::cout << "n1 : x : " << config["n1"]["x"].as<double>() << "\n";
  std::cout << "n1 : y : " << config["n1"]["y"].as<double>() << "\n";

  return 0;
}
