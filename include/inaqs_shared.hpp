#ifndef INAQS_SHARED_HPP
#define INAQS_SHARED_HPP

#include <armadillo>
#include <typeinfo>
#include <cxxabi.h>

#define INAQS_LOG(SHARED, MSG) do{SHARED->log(this, MSG);}while(0)

class INAQSShared{
public:
  INAQSShared(int mdTimeStep, double classicalTimeStep, std::string h5file_name):
    initial_step{mdTimeStep},
    restart{mdTimeStep > 0},
    mdTimeStep{mdTimeStep},
    classicalTimeStep{classicalTimeStep},
    h5file_name{h5file_name} {};  
  ~INAQSShared() { };

  int get_step(void) const {return mdTimeStep;}
  double get_dtc(void) const {return classicalTimeStep;}
  
  void inc_step(void){mdTimeStep++;} // FIXME: want this hidden

  // FIXME: better to have this be a variadic function
  template <typename T>
  void log(T * component, const std::string & message){
    char * name;
    {
      int status = -1;
      name = abi::__cxa_demangle(typeid(*component).name(), NULL, NULL, &status);
    }
    logstream << "[" << name << "] " << mdTimeStep << ": " << message << std::endl;
    free(name);
  }
  
  template <typename T>
  void saveh5(const T &tensor, std::string kind) const{
    if (! tensor.eval().save(arma::hdf5_name(h5file_name, kind + "/" +
                                             std::to_string(mdTimeStep),
                                             arma::hdf5_opts::replace)
                             )){
      throw std::runtime_error("Could not write " + kind + " at " + std::to_string(mdTimeStep) + "!");
    }
  };

  template <typename T>
  void loadh5(T &tensor, std::string kind, int idx) const{
    if (! tensor.load(arma::hdf5_name(h5file_name, kind + "/" + std::to_string(idx)))){
      throw std::runtime_error("Could not read " + kind + " at " + std::to_string(idx) + " from '" + h5file_name + "'!");
    }
  };

public:
  const int initial_step;
  const bool restart;
  
private:
  int mdTimeStep;
  const double classicalTimeStep;
  std::string h5file_name;
  std::ostream & logstream = std::cerr;
};

#endif
