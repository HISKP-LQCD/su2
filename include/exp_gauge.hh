#ifndef Genz
#pragma once

#include "su2.hh"
#include "partitionings.hh"
#include "su3.hh"
#include "u1.hh"
#include "adjointfield.hh"

inline _u1 exp(adjointu1<double> const & x) {
  return _u1(x.geta());
}

_su2 exp(adjointsu2<double> const & x);
_su3 exp(adjointsu3<double> const & x);
#endif