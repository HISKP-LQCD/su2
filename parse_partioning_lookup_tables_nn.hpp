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
    std::string file_contents = readFileIntoString( "lookuptable_nn.csv");
    //std::string file_contents_multiplication = readFileIntoString("lookup_table_multiplication.csv");
    //std::string file_contents_addition = readFileIntoString("lookup_table_addition.csv");
    
    std::istringstream sstream(file_contents);
    std::string record;

    int linecounter = 0;
    while (std::getline(sstream, record)) {
      std::istringstream line(record);
      std::cout << "reading main " << "\n";
      int wordcounter = 0;
      std::vector<size_t> neighbors;
      while (std::getline(line, record, delimiter)){
          std::cout << "linecounter"  << linecounter << "\n";
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
            
            if (record != ""){
            neighbors.push_back(std::stol(record));}
          else {
            break;
          }
          }
          
         std::cout << "wordcounter " << wordcounter << "\n"; 
          wordcounter += 1;
      }
      std::cout << "hello there \n";
      partitioning_nn::nn_lookup.push_back(neighbors);
      std::cout << "i am still here \n";
      neighbors.clear();
      linecounter += 1;
    }
  std::cout << "stain alive \n";
  }

}
 //end namespace read_partitionings