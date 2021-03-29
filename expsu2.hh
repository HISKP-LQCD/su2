#pragma once

#include "su2.hh"
#include "u1.hh"
#include "adjointfield.hh"

_su2 exp(adjointsu2<double> const & x);

inline _u1 exp(adjointu1<double> const & x) {
  return _u1(x.geta());
}
