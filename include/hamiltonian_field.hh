#pragma once

#include "adjointfield.hh"
#include "gaugeconfig.hh"
#include <vector>

template <typename Float, class Group> struct hamiltonian_field {
  adjointfield<Float, Group> *momenta;
  gaugeconfig<Group> *U;

  // There's a difference when the the theory is abelian.
  // In this case the factor in front of the force changes:
  // the RHS of eq. 8.42 of Gattringer&Lang is multiplied by 2, because the only generator
  // is the identity.
  double Fact_Nc_force = 1.0 / 2.0; // default, non-abelian case

  // Note:
  // - for the U(1) case, get_deriv gives: Im(S) = -(i/2) * (S - S^\dagger)
  // - for the SU(2) case, get_deriv gives: (-i)*(S - S^\dagger)
  double Fact_get_deriv = 1.0 / 2.0; // default: SU(2)

  hamiltonian_field(adjointfield<Float, Group> &momenta, gaugeconfig<Group> &U)
    : momenta(&momenta), U(&U) {
    if (U.getNc() == 1) {
      Fact_Nc_force = 1.0;
      Fact_get_deriv = 1.0;
    }
  }
};
