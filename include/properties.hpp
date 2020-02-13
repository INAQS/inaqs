#ifndef __PROPERTIES_HPP
#define __PROPERTIES_HPP

#include <unordered_map>
#include <armadillo>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
enum class QMProperty{
    nacvector,
    nacvector_imag,
    wfoverlap,
    qmgradient,
    mmgradient,
    energies,
};


class PropMap
{
public:
  PropMap() : prop{} {};

  arma::cube* get(QMProperty key);
  const arma::uvec* get_idx(QMProperty key) const;
  inline arma::cube& operator[](QMProperty key) { return *get(key); }  //get(key)
  
  void emplace(QMProperty p, arma::cube* vec);
  void emplace(QMProperty p, arma::uvec iv, arma::cube* vec);

  bool has(QMProperty key) const;
  bool has_idx(QMProperty key) const;

  const std::vector<QMProperty> keys(void) const;
  
private:
  std::unordered_map<QMProperty, arma::cube*> prop{};
  std::unordered_map<QMProperty, arma::uvec> prop_vec{};
};

#endif
