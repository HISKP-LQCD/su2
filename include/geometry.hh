#pragma once


class geometry { 
public:
  explicit geometry(size_t _Ls, size_t _Lt) : Ls(_Ls), Lt(_Lt) {}
  size_t getLs() const {
    return(Ls);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getIndex(const int t, const int x, const int y, const int z) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Ls) % Ls;
    size_t y2 = (y + Ls) % Ls;
    size_t y3 = (z + Ls) % Ls;
    return( ((t*Ls + x)*Ls + y)*Ls + z );
  }
  void getCoordinate(size_t c[], const size_t index) const {
    size_t x0 = index, x = index/Ls;
    for(int i = 3; i >= 0; i--) {
      c[i] = x0-x*Ls;
      x0 = x;
      x /= Ls;
    }
  }

private:
  size_t Ls, Lt;

};

