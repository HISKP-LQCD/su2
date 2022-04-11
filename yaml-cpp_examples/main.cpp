// a.cpp

//#include "yce.h"
#include <iostream>
#include "yaml.h"


void find_all(const YAML::Node& node){
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it){
    std::string k = it->first.as<std::string>(); // key
    std::cout << "First: " << k << "\n";
    if(node[k]){
      find_all(node[k]);
    }
    
    // it->second.as<std::string>(); // can't do this until it's type is checked!!
  }
}

int main(int argc, char const *argv[]) {

  const YAML::Node config = YAML::LoadFile("./config.yaml");
  std::cout << "p1 : " << config["p1"].as<double>() << "\np2 : " << config["p2"].as<std::string>() << "\n";

  std::cout << "n1 : x : " << config["n1"]["x"].as<double>() << "\n";
  std::cout << "n1 : y : " << config["n1"]["y"].as<double>() << "\n";

  const YAML::Node node = YAML::LoadFile("./config.yaml");
  std::cout << node.Type() << " " << YAML::NodeType::Null << "\n";

  find_all(node);

  switch (node.Type()) {
    case YAML::NodeType::Null: // ...
    case YAML::NodeType::Scalar: // ...
    // case YAML::NodeType::Sequence:
    //   for (auto it = node.begin(); it != node.end(); ++it) {
    //     auto element = *it;
    //     // recurse on "element"
    //   }
    //   break;
    // case Map:
    //   for (auto it = node.begin(); it != node.end(); ++it) {
    //     auto key = it->first;
    //     auto value = it->second;
    //     // recurse on "key" and "value"
    //     // if you're sure that "key" is a string, just grab it here
    //   }
    break;
    case YAML::NodeType::Undefined: std::cout << "ciao\n";
  }
  return 0;
}
