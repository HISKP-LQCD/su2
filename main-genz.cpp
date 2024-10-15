
/** 
* @file main-genz.cpp
* @author Sebastian Mueller (s6sbmuel@uni-bonn.de)
* @brief main program running any simulation of this libary with the Genz Point partitioning
* @version  0.1
* @date 2024-04-09


*/



#define Genz
#include "run_program.hpp"
int main(int argc, char *argv[]){
    typedef _Gsu2 Group;
    run_program<Group>(argc, argv);
    
    return (0);
}