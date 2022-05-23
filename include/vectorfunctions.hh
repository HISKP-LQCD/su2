#pragma once
#include<vector>

//elementwise addition, if vectors are of the same length
template<class numbertype=double> void operator+=(std::vector<numbertype> &v1, std::vector<numbertype> v2) {
  if(v1.size() != v2.size()){
      std::cerr << "vectors do not match in size! No addition took place." << std::endl;
      return;
  }
  for(size_t i = 0; i < v1.size(); i++){
    v1[i] += v2[i];
  }
}

//returns a vector of length length filled with zeros
template<class numbertype=double> std::vector<numbertype> zerovector(size_t length){
    std::vector<numbertype> res(length); //zero is implicit
    return res;
}
