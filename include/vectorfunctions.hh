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
