#include "properties.hpp"
#include <algorithm>
#include <armadillo>

/*
  Using switch to get complaints from compiler when we miss one
*/
std::ostream& operator<<( std::ostream& oss, const QMProperty p)
{
  std::string name;
  switch (p){
  case QMProperty::nacvector:
    name = "nacvector";
    break;
  case QMProperty::nacvector_imag:
    name = "nacvector_imag";
    break;
  case QMProperty::wfoverlap:
    name = "wfoverlap";
    break;
  case QMProperty::qmcharge:
    name = "qmcharge";
    break;
  case QMProperty::diabatic_rot_mat:
    name = "diabatic_rot_mat";
    break;
  case QMProperty::diabatic_H:
    name = "diabatic_H";
    break;
  case QMProperty::diabatic_gradients:
    name = "diabatic_gradients";
    break;
  case QMProperty::qmgradient:
    name = "qmgradient";
    break;
  case QMProperty::mmgradient:
    name = "mmgradient";
    break;
  case QMProperty::qmgradient_multi:
    name = "qmgradient_multi";
    break;
  case QMProperty::mmgradient_multi:
    name = "mmgradient_multi";
    break;
  case QMProperty::energies:
    name = "energies";
    break;
  }
  
  oss << name;
  return oss;
}


ArmaWrap PropMap::get(QMProperty key) {
  auto itr = prop.find(key);
  if(itr == prop.end()){
    return nullptr;
  }
  else{ return itr->second;}
}

const arma::uvec* PropMap::get_idx(QMProperty key) const { 
  auto itr = prop_vec.find(key);
  if (itr == prop_vec.end()) {
    return nullptr; 
  }
  return &itr->second;
}

arma::uvec & PropMap::get_writable_idx(QMProperty key) {
  auto itr = prop_vec.find(key);
  if (itr == prop_vec.end()) {
    throw std::logic_error("get_writable_idx() is dangerous and must be used with has_idx()!");
  }
  return itr->second;
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
