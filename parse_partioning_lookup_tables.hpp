/**
 * @file parse_partioning_lookup_tables.hpp
 * @author Sebastian Mueller (s6sbmuel@uni-bonn.de)
 * @brief loads the lookup tables for the partitioning.
 * Their names are fixed since they are loaded before the yaml file
 * @version 0.1
 * @date 2024-10-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */


#pragma once
#include "partitionings.hh"
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sstream>
#include <vector>


namespace read_partitionings {
  

  std::string readFileIntoString(const std::string& path) {
  auto ss = std::ostringstream{};
  std::ifstream input_file(path);
  //input_file.open(path);
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
    std::string file_contents = readFileIntoString( "lookuptable1.csv");
    //std::cout << file_contents << "\n";
    std::string file_contents_multiplication = readFileIntoString("lookup_table_multiplication.csv");
    std::string file_contents_addition = readFileIntoString("lookup_table_addition.csv");
    
    std::istringstream sstream(file_contents);
    std::istringstream multiplicationsstream(file_contents_multiplication);
    std::istringstream additionsstream(file_contents_addition);
    std::string record;

    int linecounter = 0;
    int header = 1;
    while (std::getline(sstream, record)) {
      std::istringstream line(record);
      std::cout << "reading main " << "\n";
      int wordcounter = 0;
      
      while (std::getline(line, record, delimiter)){
          if (linecounter >=header){
          if (wordcounter == 1){
            partitioning::point0.push_back(std::stod(record));
          }
          else if (wordcounter == 2){
            partitioning::point1.push_back(std::stod(record)); 
          }
          else if (wordcounter == 3) {
            partitioning::point2.push_back(std::stod(record));
          }
          else if (wordcounter == 4) {
            partitioning::point3.push_back(std::stod(record));
          }
          else if (wordcounter == 5){
            partitioning::weights.push_back(std::stod(record));
          }
          else if (wordcounter == 6){
            partitioning::distance_to_identity.push_back(std::stod(record));
          }
          else if (wordcounter == 7){
            //std::cout << record << " \n";
            partitioning::dagger_vector.push_back(std::stol(record));
          }
          }
          wordcounter += 1;
      }
      
      linecounter += 1;
    }
//std::cout << "got here" << "\n";
  linecounter = 0;
  while (std::getline(multiplicationsstream, record)){
    std::istringstream multiplicationline(record);
    //std::cout << record << "n";
    int wordcounter = 0;
    while(std::getline(multiplicationline, record, delimiter)){
      //std::cout <<linecounter << "\n";
      if (linecounter == 0){
      
        std::vector<size_t> create_word_vector; 
        partitioning::multiplication_lookup_table.push_back(create_word_vector);
        
      }
      else if ((linecounter < header) && (linecounter != 0)){
        //std::cout << "did not stay dead multi" << "\n";
      }
      else {
        if (wordcounter > 0){
        partitioning::multiplication_lookup_table[wordcounter-1].push_back(std::stol(record));
        //std::cout << record << " \n"; 
        }
      }
      //std::cout << record << " \n";
      wordcounter += 1;

    }
    linecounter += 1;
  }
  //std::cout <<header << "\n";
 std::cout << "stain alive" << "\n"; 
  linecounter = 0;
  while (std::getline(additionsstream, record)){
    std::istringstream additionline(record);
    int wordcounter = 0;
    while (std::getline(additionline, record, delimiter)){
      if (linecounter == 0){
        std::vector<size_t> create_word_vector;
        partitioning::addition_lookup_table.push_back(create_word_vector);
      }
      else if ((linecounter < header) && (linecounter != 0)){ 
        //std::cout << "did not stay dead addition" << "\n";
      }
      else{
        if (wordcounter > 0){
        //std::cout << record << "\n";
        partitioning::addition_lookup_table[wordcounter-1].push_back(std::stol(record));
        }
      }
      wordcounter += 1;
    }
    linecounter += 1;
  }

}
} //end namespace read_partitionings