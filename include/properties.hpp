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
    energies,
    nacvector,
    nacvector_imag,
    wfoverlap,
    
};


class PropMap
{
public:
  PropMap() : prop{} {};

  std::vector<double>* get(QMProperty key);
  const std::vector<int>* get_idx(QMProperty key) const;
  inline std::vector<double>& operator[](QMProperty key) { return *prop[key]; }  //get(key)
  
  void emplace(QMProperty p, std::vector<double>* vec);
  void emplace(QMProperty p, std::vector<int> iv, std::vector<double>* vec);
  
private:
  std::unordered_map<QMProperty, std::vector<double>*> prop{};
  std::unordered_map<QMProperty, std::vector<int>> prop_vec{};
};

#endif
