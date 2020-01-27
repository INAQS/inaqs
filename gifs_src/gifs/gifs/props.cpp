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

class  QMPropertyWrap
{
public:
    QMPropertyWrap(QMProperty in_prop) : qmprop{in_prop} {}

    inline bool operator==(QMProperty rhs) const { return rhs == qmprop; }
    inline bool operator<(QMProperty rhs) const { return rhs < qmprop; }
    inline bool operator>(QMProperty rhs) const { return rhs > qmprop; }
    inline bool operator==(QMPropertyWrap rhs) const { return rhs.qmprop == qmprop; }
    inline bool operator<(QMPropertyWrap rhs) const { return rhs.qmprop < qmprop; }
    inline bool operator>(QMPropertyWrap rhs) const { return rhs.qmprop > qmprop; }

    QMProperty qmprop{};
    virtual const std::vector<int>* get_idx(void) const {return nullptr;} 
};

namespace std
{
    template<> struct hash<QMPropertyWrap>
    {
        std::size_t operator()(QMPropertyWrap const& prop) const noexcept
        {
            return std::hash<QMProperty>{}(prop.qmprop);
        }

    };
}


class QMPropertyVector:
    public QMPropertyWrap
{
public:
  QMPropertyVector(QMProperty p, std::vector<int> in_idx) : QMPropertyWrap{p}, idx{in_idx} {}
  virtual const std::vector<int>* get_idx(void) const {return &idx;} 
private:
  const std::vector<int> idx;
};


class PropMap : 
    public std::unordered_map<QMPropertyWrap, std::vector<double>*>
{
public:
  PropMap() : std::unordered_map<QMPropertyWrap, std::vector<double>*>{} {};
  
  std::vector<double>* get(QMPropertyWrap key) {
    auto itr = find(key);
    if(itr == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return itr->second;
  }
};

QMPropertyWrap mkprop(QMProperty prop) { return QMPropertyWrap{prop}; }
QMPropertyWrap mkprop(QMProperty prop, std::vector<int> idx) { return QMPropertyVector(prop, idx); }

int main(void){
  std::vector<double> u{1, 2, 3}; 
  std::vector<double> v{4, 5, 6};
  std::vector<int> grads{1, 2, 3};
  
  PropMap props{};
  props.emplace(QMProperty::energies, &v);
  props.emplace(mkprop(QMProperty::qmgradient, grads), &u);
  //props.emplace(mkprop(QMProperty::qmgradient), &v);
  
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
      std::cout << "size: " << (*id).size() << std::endl;
    }
  }
  
  return 0;
}
