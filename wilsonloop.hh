// Copyright (C) 2017 C. Urbach

#pragma once

#include"su2.hh"
#include"gaugeconfig.hh"
#include<string>

double planar_wilsonloop_dir(gaugeconfig &U, const size_t r, const size_t t, 
                      const size_t mu, const size_t nu);
double wilsonloop(gaugeconfig &U, const size_t r, const size_t t);
void compute_all_loops(gaugeconfig &U, std::string const &path);
void compute_spacial_loops(gaugeconfig &U, std::string const &path);
