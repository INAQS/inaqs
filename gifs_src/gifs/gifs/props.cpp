#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
enum class QMProperty {
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
  
  std::vector<double>* get(QMProperty key) {
    auto itr = prop.find(key);
    if(itr == prop.end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return itr->second;
  }


  void emplace(QMProperty p, std::vector<double>* vec) { prop.emplace(p, vec); }
  void emplace(QMProperty p, std::vector<int> iv, std::vector<double>* vec) { prop.emplace(p, vec); prop_vec.emplace(p, iv); }
  //
  inline std::vector<double>& operator[](QMProperty key) { return *prop[key]; } 
  //
  const std::vector<int>* get_idx(QMProperty key) const { 
      auto itr = prop_vec.find(key);
      if (itr == prop_vec.end()) {
          return nullptr; 
      }
      return &itr->second;
  }
private:
    std::unordered_map<QMProperty, std::vector<double>*> prop{};
    std::unordered_map<QMProperty, std::vector<int>> prop_vec{}; 
};


int main(void){
  std::vector<double> u{1, 2, 3}; 
  std::vector<double> v{4, 5, 6};
  std::vector<int> grads{1, 2, 3};
  std::vector<int> index{1, 2, 3};
  
  PropMap props{};
  props.emplace(QMProperty::qmgradient, grads, &u);
 
  std::cout << props[QMProperty::qmgradient][1] << '\n';
  std::cout << props.get_idx(QMProperty::qmgradient)->size() << '\n';

  /*

  const std::vector<int>* idx = pv.get_idx();
  std::cout << (idx == nullptr) << "\n";
  //std::cout << "size: " << (*idx).size() << std::endl;
  
  PropMap props{};
  props.emplace(mkprop(QMProperty::qmgradient, grads), &u);
  props.emplace(QMProperty::energies, &v);             // Comment this line and uncomment the below
  //props.emplace(mkprop(QMProperty::qmgradient), &v); // a 'Bad Key!' abort
  
  std::vector<double> *out = props.get(QMProperty::energies);
  std::cout << v[1] << "\n";
  (*out)[1] = 3;
  std::cout << v[1] << "\n";
  std::cout << (*out)[1] << "\n";

  auto s = props.find(QMProperty::qmgradient);
  if (s != props.end()){
    std::cout << "found it!" << std::endl;
    //const QMPropertyVector qmpv = dynamic_cast<QMPropertyVector> (s->first);
    const std::vector<int> *id = s->first.get_idx();
    if ( id == nullptr){
      std::cout << "but didn't find idx!" << std::endl;
    }
    else{
//      std::cout << "size: " << (*id).size() << std::endl;
    }
  }
  
  */
  return 0;
}
