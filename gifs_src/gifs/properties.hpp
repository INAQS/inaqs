#ifndef __PROPERTIES_HPP
#define __PROPERTIES_HPP

#include <unordered_map>
#include <vector>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
enum class QMProperty{
    qmgradient,
    mmgradient,
    energy,
};

class PropMap:
  public std::unordered_map<QMProperty, std::vector<double>*>{
public:
  PropMap() : std::unordered_map<QMProperty, std::vector<double>*>{} {}
  std::vector<double>& get(QMProperty key);
};

#endif
