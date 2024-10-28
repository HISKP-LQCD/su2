#define partinn
#include "partitionings_nn.hh"
#include "parse_partioning_lookup_tables_nn.hpp"

 #include "run_program.hpp"
 int main(int argc, char *argv[]){
    read_partitionings_nn::load_tables();
    typedef _partitioning_nn Group;

    run_program<Group>(argc, argv);
    return(0);
 }