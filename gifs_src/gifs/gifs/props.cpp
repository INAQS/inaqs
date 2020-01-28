#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

/*
  List of all possible properties that could be produced by or
  requested of a QMInterface.
*/
enum class QMPropertyEnum {
    qmgradient,
    mmgradient,
    energies,
    nacvector,
    nacvector_imag,
    wfoverlap,
};

class  QMProperty 
{
public:
    QMProperty(QMPropertyEnum in_prop) : qmprop{in_prop} {}

    inline bool operator==(QMPropertyEnum rhs) const { return rhs == qmprop; }
    inline bool operator<(QMPropertyEnum rhs) const { return rhs < qmprop; }
    inline bool operator>(QMPropertyEnum rhs) const { return rhs > qmprop; }
    inline bool operator==(QMProperty rhs) const { return rhs.qmprop == qmprop; }
    inline bool operator<(QMProperty rhs) const { return rhs.qmprop < qmprop; }
    inline bool operator>(QMProperty rhs) const { return rhs.qmprop > qmprop; }

    QMPropertyEnum qmprop{};
};

namespace std
{
    template<> struct hash<QMProperty>
    {
        std::size_t operator()(QMProperty const& prop) const noexcept
        {
            return std::hash<QMPropertyEnum>{}(prop.qmprop);
        }
    };
}


class QMPropertyVector:
    public QMProperty
{
public:
  QMPropertyVector(QMPropertyEnum p, std::vector<int> in_idx) : QMProperty{p}, idx{in_idx} {}
private:
  const std::vector<int> idx;
};


class PropMap : 
    public std::unordered_map<QMProperty, std::vector<double>*>
{
public:
  PropMap() : std::unordered_map<QMProperty, std::vector<double>*>{} {};
  
  std::vector<double>* get(QMProperty key) {
    auto itr = find(key);
    if(itr == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return itr->second;
  }
};

QMProperty get_property(QMPropertyEnum prop) { return QMProperty{prop}; }
QMProperty get_property(QMPropertyEnum prop, std::vector<int> idx) { return QMPropertyVector(prop, idx); }

int main(void){
  std::vector<double> u{1, 2, 3}; 
  std::vector<double> v{4, 5, 6};
  std::vector<int> grads{1, 2, 3};
  
  PropMap props{};
  props.emplace(QMPropertyEnum::energies, &v);
  props.emplace(get_property(QMPropertyEnum::qmgradient, grads), &u);

  std::vector<double>* out = props.get(QMPropertyEnum::energies);
  std::cout << v[1] << "\n";
  (*out)[1] = 3;
  std::cout << v[1] << "\n";
  std::cout << (*out)[1] << "\n";
  
  return 0;
}
