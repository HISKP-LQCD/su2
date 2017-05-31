#pragma once

#include"gaugeconfig.hh"
#include<vector>

double compute_hamiltonian(std::vector<double> const &momenta,
                           gaugeconfig & U);
