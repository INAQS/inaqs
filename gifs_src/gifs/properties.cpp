#include "properties.hpp"
#include <vector>

std::vector<double> PropMap::get(QMProperty key);
{
    auto itr = find(key);
    if(itr == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return itr->second;
}
