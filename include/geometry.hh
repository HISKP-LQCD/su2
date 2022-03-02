#pragma once

#include <cstddef>

class geometry { 
public:
  explicit geometry(const size_t _Lx, const size_t _Ly,
                    const size_t _Lz, const size_t _Lt) :
    Lx(_Lx), Ly(_Ly), Lz(_Lz), Lt(_Lt) {}
  size_t getLx() const {
    return(Lx);
  }
  size_t getLy() const {
    return(Ly);
  }
  size_t getLz() const {
    return(Lz);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getIndex(const int t, const int x, const int y, const int z) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    return( ((y0*Lx + y1)*Ly + y2)*Lz + y3 );
  }

  // *** Ls variable is not defined ***
  // void getCoordinate(size_t c[], const size_t index) const {
  //   size_t x0 = index, x = index/Ls;
  //   for(int i = 3; i >= 0; i--) {
  //     c[i] = x0-x*Ls;
  //     x0 = x;
  //     x /= Ls;
  //   }
  // }

private:
  size_t Lx, Ly, Lz, Lt;

};

