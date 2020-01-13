#ifndef __PROPERTIES_BASE_HPP
#define __PROPERTIES_BASE_HPP

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/

std::vector<double>& PropMap::get(QMProperty key) {
    auto itr = find(key);
    if(itr == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return *(itr->second);
};

#endif
