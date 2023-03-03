#ifndef __PROPERTIES_HPP
#define __PROPERTIES_HPP

#include <unordered_map>
#include <armadillo>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
enum class QMProperty{
   nacvector,          // mat
   nacvector_imag,     // mat
   wfoverlap,          // mat
   qmcharge,           // vec
   diabatic_rot_mat,   // mat
   diabatic_H,         // mat
   diabatic_gradients, // cube(3, N, 2)
   qmgradient,         // mat
   mmgradient,         // mat
   qmgradient_multi,   // cube
   mmgradient_multi,   // cube
   energies,           // vec
};

std::ostream& operator<<( std::ostream& oss, const QMProperty p);

/*
  A wrapper around pointers to various armadillo base types so that
  PropMap can store tensors of different rank. Implicitly converts
  back to the base type on assignment; will raise runtime exceptions
  if invalid conversions are attempted.

  To add support for another type:
  * Add a pointer member
  * Add a constructor
  * Add operators for & and *
  * Update the function std::string kind(void) with the new type
*/
class ArmaWrap{
public:
  ArmaWrap(std::nullptr_t x) {(void) x;};
  ArmaWrap(arma::vec  * x): v(x) {};
  ArmaWrap(arma::mat  * x): m(x) {};
  ArmaWrap(arma::cube * x): c(x) {};

  // FIXME: should the nullptr_t return be conditional?
  operator std::nullptr_t() const {return nullptr;};
  operator arma::vec  *() const {return notnull(v);};
  operator arma::mat  *() const {return notnull(m);};
  operator arma::cube *() const {return notnull(c);};

  operator arma::vec  &() const {return *notnull(v);};
  operator arma::mat  &() const {return *notnull(m);};
  operator arma::cube &() const {return *notnull(c);};

  ArmaWrap& operator= (const arma::vec&  other) {*notnull(v) = other; return *this;}
  ArmaWrap& operator= (const arma::mat&  other) {*notnull(m) = other; return *this;}
  ArmaWrap& operator= (const arma::cube& other) {*notnull(c) = other; return *this;}

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
    if      (v){return "vec";}
    else if (m){return "mat";}
    else if (c){return "cube";}
    else       {return "uninitialized wrapper";}
  }
};


class PropMap
{
public:
  PropMap() : prop{} {};
  ArmaWrap get(QMProperty key);
  const arma::uvec* get_idx(QMProperty key) const;
  arma::uvec & get_writable_idx(QMProperty key);
  
  void emplace(QMProperty p, ArmaWrap vec);
  void emplace(QMProperty p, arma::uvec iv, ArmaWrap vec);

  bool has(QMProperty key) const;
  bool has_idx(QMProperty key) const;

  const std::vector<QMProperty> keys(void) const;

  friend std::ostream& operator<<(std::ostream& oss, const PropMap &pm){
    oss << "{ ";
    for (const auto& e: pm.prop){
      oss << "{" << e.first;
      if (pm.has_idx(e.first)){
        auto idx = pm.get_idx(e.first);
        oss << ": ";
        for (const auto& v: *idx){
          oss << v << " ";
        }
      }
      oss << "}, ";
    }
    oss << "}";
    return oss;
  }

  
private:
  std::unordered_map<QMProperty, ArmaWrap> prop{};
  std::unordered_map<QMProperty, arma::uvec> prop_vec{};
};


#endif
