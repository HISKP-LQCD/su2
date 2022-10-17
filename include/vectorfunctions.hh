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

//interprets the vector as a space-time coordinate and inverts the spatial cordinates
template<class numbertype=size_t> std::vector<numbertype> invertspace(const std::vector<numbertype> &x){
    std::vector<numbertype> res = x;
    res[1] = -x[1];
    res[2] = -x[2];
    res[3] = -x[3];
    return res;
}

template<class numbertype=size_t> std::vector<numbertype> xminusmu(const std::vector<numbertype> &x, const size_t mu){
    std::vector<numbertype> res=x;
    res[mu]--;
    return res;
}



//This seems to be already implemented in the standard library and unneeded
//~ //compares two vectors
//~ template<class numbertype=double> bool operator==(const std::vector<numbertype> &v1, const std::vector<numbertype> &v2){
  //~ if(v1.size() != v2.size()){
    //~ std::cerr << "vectors do not match in size! Comparison is not possible." << std::endl;
    //~ abort();
  //~ }  
  //~ for(size_t i = 0; i < v1.size(); i++){
    //~ if(v1[i] != v2[i]){
      //~ return false;
    //~ }
  //~ }
  //~ return true;
//~ }
//~ //compares two vectors
//~ template<class numbertype=double> bool operator!=(const std::vector<numbertype> &v1, const std::vector<numbertype> &v2){
  //~ return (!(v1==v2));
//~ }

//~ //multiply a vector by a scalar
//~ template<class numbertype=double> bool operator!=(const std::vector<numbertype> &v1, const std::vector<numbertype> &v2){
  //~ return (!(v1==v2));
//~ }
