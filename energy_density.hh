#pragma once
#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"

void energy_density(gaugeconfig<su2> &U, double &res, double &Q);
void energy_density(gaugeconfig<_u1> &U, double &res, double &Q);
