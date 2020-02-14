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

/*
  A wrapper around pointers to various armadillo base types so that
  PropMap can store tensors of different rank. Implicitly converts
  back to the base type on assignement; will raise runtime exceptions
  if invalid conversions are attempted.

  To add support for another type:
  * Add a member
  * Add a constructor
  * Add an operator
  * Update the function std::string kind(void) with the new type
*/
class ArmaWrap{
public:
  ArmaWrap(std::nullptr_t x) {(void) x;};
  ArmaWrap(arma::vec  * x): v(x) {};
  ArmaWrap(arma::mat  * x): m(x) {};
  ArmaWrap(arma::cube * x): c(x) {};

  operator std::nullptr_t() const {return nullptr;};
  operator arma::vec  *() const {return notnull(v);};
  operator arma::mat  *() const {return notnull(m);};
  operator arma::cube *() const {return notnull(c);};

private:
  // Member types
  arma::vec  * v = nullptr;
  arma::mat  * m = nullptr;
  arma::cube * c = nullptr;

  template <typename T>
  T notnull(T p) const{
    if (!p){
      std::string msg = "Invalid conversion to " + (std::string) typeid(T).name() + "; I am a " + kind() + "!";
      throw std::logic_error(msg);
    }
    return p;
  }
  
  std::string kind(void) const{
    std::string kind;
    if      (v){kind = "vec";}
    else if (m){kind = "mat";}
    else if (c){kind = "cube";}
    else       {kind = "uninitialized wrapper";}
    return kind; 
  }
};


class PropMap
{
public:
  PropMap() : prop{} {};
  ArmaWrap get(QMProperty key);
  const arma::uvec* get_idx(QMProperty key) const;
  //inline arma::cube& operator[](QMProperty key) { return *get(key); }  //get(key)
  
  void emplace(QMProperty p, ArmaWrap vec);
  void emplace(QMProperty p, arma::uvec iv, ArmaWrap vec);

  bool has(QMProperty key) const;
  bool has_idx(QMProperty key) const;

  const std::vector<QMProperty> keys(void) const;
  
private:
  std::unordered_map<QMProperty, ArmaWrap> prop{};
  std::unordered_map<QMProperty, arma::uvec> prop_vec{};
};

#endif
