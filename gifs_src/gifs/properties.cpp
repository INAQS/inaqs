#include "properties.hpp"
#include <algorithm>
#include <armadillo>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
ArmaWrap PropMap::get(QMProperty key) {
  auto itr = prop.find(key);
  if(itr == prop.end()){
    return nullptr;
  }
  else return itr->second;
}

const arma::uvec* PropMap::get_idx(QMProperty key) const { 
  auto itr = prop_vec.find(key);
  if (itr == prop_vec.end()) {
    return nullptr; 
  }
  return &itr->second;
}

bool PropMap::has(QMProperty key) const{
  return prop.end() != prop.find(key); 
}

bool PropMap::has_idx(QMProperty key) const{
  return prop_vec.end() != prop_vec.find(key); 
}

const std::vector<QMProperty> PropMap::keys(void) const {
  std::vector<QMProperty> qv;
  for (auto p: prop){
    qv.push_back(p.first);
  }
  std::sort(qv.begin(), qv.end());
  return qv;
}

void PropMap::emplace(QMProperty p, ArmaWrap vec) {
  if (has(p) || has_idx(p)){
    throw std::invalid_argument("Property already in Map!");
  }
  else{
    prop.emplace(p, vec);
  }
}

void PropMap::emplace(QMProperty p, arma::uvec iv, ArmaWrap vec) {
  if (has(p) || has_idx(p)){
    throw std::invalid_argument("Property already in Map!");
  }
  else{
    prop.emplace(p, vec);
    prop_vec.emplace(p, iv);
  }
}
