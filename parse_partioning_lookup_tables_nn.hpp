/**
 * @file parse_partioning_lookup_tables_nn.hpp
 * @author Sebastian Müller (s6sbmuel@uni-bonn.de)
 * @brief loads the lookup table for the neighrest neighbor implementation of the partitioning
 * @version 0.1
 * @date 2024-10-30
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once
#include "partitionings_nn.hh"
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sstream>
#include <vector>


namespace read_partitionings_nn {
  
  std::string readFileIntoString(const std::string& path) {
  auto ss = std::ostringstream{};
  std::ifstream input_file(path);
  std::cout << input_file.is_open() << std::endl;
  if (!input_file.is_open()) {
    std::cerr << "Could not open the file - '" << path << "'" << std::endl;
    exit(EXIT_FAILURE);
  }
  ss << input_file.rdbuf();
  input_file.close();
  return ss.str();
}

void load_tables(){
    
    
    char delimiter = ',';
    std::string file_contents = readFileIntoString( "lookuptable_nn.csv"); //filename fixed since this is loaded before the YAML file
    std::istringstream sstream(file_contents);
    std::string record;
    int linecounter = 0;
    while (std::getline(sstream, record)) {
      std::istringstream line(record);
      std::cout << "### reading lookup talbe ### " << "\n";
      int wordcounter = 0;
      std::vector<size_t> neighbors;
      while (std::getline(line, record, delimiter)){
          if (wordcounter == 1){
            partitioning_nn::point0.push_back(std::stod(record));
          }
          else if (wordcounter == 2){
            partitioning_nn::point1.push_back(std::stod(record)); 
          }
          else if (wordcounter == 3) {
            partitioning_nn::point2.push_back(std::stod(record));
          }
          else if (wordcounter == 4) {
            partitioning_nn::point3.push_back(std::stod(record));
          }
          else if (wordcounter == 5){
            partitioning_nn::weights.push_back(std::stod(record));
          }
          else if (wordcounter > 5){
            
            if (record != ""){ // "" is the represnetation of NaN and thus the end of the neighrest neighbor list
            neighbors.push_back(std::stol(record));}
          else { // breaks the loop if end of neighrest neighbor list
            break;
          }
          }
          
          wordcounter += 1;
      }
      partitioning_nn::nn_lookup.push_back(neighbors);
      neighbors.clear();
      linecounter += 1;
    }
  }

}
 //end namespace read_partitionings