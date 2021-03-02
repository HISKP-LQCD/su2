#pragma once
#include"su2.hh"
#include"gaugeconfig.hh"
#include<string>

void gradient_flow(gaugeconfig<su2> &U, std::string const &path, const double tmax = 9.99);
