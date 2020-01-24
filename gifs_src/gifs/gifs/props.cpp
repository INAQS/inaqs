#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iostream>

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

class QMReq{
public:
  QMReq(QMProperty p, std::vector<double> &r): prop(p), r(r), idx({}) {}
  QMReq(QMProperty p, std::vector<double> &r, const std::vector<int> idx): prop(p), r(r), idx(idx) {}
  
  QMProperty prop;
  std::vector<double> &r;
  const std::vector<int> idx; 
};

class PropMap : public std::unordered_map<QMProperty, std::vector<double>*>
{
public:
  PropMap() : std::unordered_map<QMProperty, std::vector<double>*>{} {};
  
  std::vector<double>& get(QMProperty key) {
    auto itr = find(key);
    if(itr == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return *(itr->second);
  }
};

class PropVec : public std::vector<QMReq>{
public:
  PropVec() : std::vector<QMReq> {} {};
  
  bool includes(QMProperty key){
    auto it = std::find_if(begin(), end(),
			   [&] (const QMReq &p)->bool{
			     return p.prop == key;});
    return it != end();
  }
  
  std::vector<double>& getmem(QMProperty key){
    auto it = std::find_if(begin(), end(),
			   [&] (const QMReq &p)->bool{
			     return p.prop == key;});
    if(it == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return it->r;
  }

  const std::vector<int>& getidx(QMProperty key){
    auto it = std::find_if(begin(), end(),
			   [&] (const QMReq &p)->bool{
			     return p.prop == key;});
    if(it == end()){
      throw std::invalid_argument("Bad Key!");
    }
    else return it->idx;
  }
};

int main(void){
  std::vector<double> u, v, w;
  
  PropMap props{};
  props.emplace(QMProperty::qmgradient, &v);
  props.emplace(QMProperty::mmgradient, &u);

  PropVec props_v {};
  props_v.push_back(QMReq(QMProperty::qmgradient, v));
  props_v.push_back(QMReq(QMProperty::mmgradient, u, {1, 2, 3}));
  props_v.push_back(QMReq(QMProperty::energies,   w, {0, 2, 3}));

  const std::vector<QMProperty> list_of_props= {QMProperty::qmgradient,
						QMProperty::mmgradient,
						QMProperty::energies,
						QMProperty::nacvector,
						QMProperty::nacvector_imag,
						QMProperty::wfoverlap};
  for (const QMProperty p: list_of_props){
    std::cout << (int) p << ": ";
    if (props_v.includes(p)){
      std::cout << props_v.includes(p) << " :: ";
      for (const auto &i: props_v.getidx(p)){
	std::cout << i << " " ;
      }
    }
    std::cout << std::endl;
  }
  
  return 0;
}
