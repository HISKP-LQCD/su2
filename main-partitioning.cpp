/**
 * @file main-partitioning.cpp
 * @author Sebastian Müller (s6sbmuel@uni-bonn.de)
 * @brief main programm running any simulations of this library for the partitionings
 *          More details about them can be found in https://doi.org/10.1140/epjc/s10052-022-10192-5
 * @version 0.1
 * @date 2024-10-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */
 #define parti
#include "partitionings.hh"
#include "parse_partioning_lookup_tables.hpp"

 #include "run_program.hpp"
 int main(int argc, char *argv[]){
      read_partitionings::load_tables();
    typedef _partitioning Group;

    run_program<Group>(argc, argv);
    return(0);
 }