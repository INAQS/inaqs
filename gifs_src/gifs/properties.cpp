#include "properties.hpp"

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
std::vector<double>* PropMap::get(QMProperty key) {
  auto itr = prop.find(key);
  if(itr == prop.end()){
    throw std::invalid_argument("Bad Key!");
  }
  else return itr->second;
}

const std::vector<int>* PropMap::get_idx(QMProperty key) const { 
  auto itr = prop_vec.find(key);
  if (itr == prop_vec.end()) {
    return nullptr; 
  }
  return &itr->second;
}

void PropMap::emplace(QMProperty p, std::vector<double>* vec) { prop.emplace(p, vec); }
void PropMap::emplace(QMProperty p, std::vector<int> iv, std::vector<double>* vec) { prop.emplace(p, vec); prop_vec.emplace(p, iv); }
