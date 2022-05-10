// yaml.cpp

#include <iostream>

#include "yaml.h"

void find_all(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    std::string k = it->first.as<std::string>(); // key
    std::cout << "First: " << k << "\n";
    if (node[k]) {
      find_all(node[k]);
    }
  }
}

int main(int argc, char const *argv[]) {
  const YAML::Node config = YAML::LoadFile("./config.yaml");
  std::cout << "p1 : " << config["p1"].as<double>()
            << "\np2 : " << config["p2"].as<std::string>() << "\n";

  std::cout << "n1 : x : " << config["n1"]["x"].as<double>() << "\n";
  std::cout << "n1 : y : " << config["n1"]["y"].as<double>() << "\n";

  const YAML::Node node = YAML::LoadFile("./config.yaml");
  std::cout << node.Type() << " " << YAML::NodeType::Null << "\n";

  find_all(node);

  return 0;
}
