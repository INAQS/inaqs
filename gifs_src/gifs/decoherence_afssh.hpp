#ifndef __GIFS_DECOHERENCE_AFSSH_HPP
#define __GIFS_DECOHERENCE_AFSSH_HPP
#include "decoherence.hpp"

/*
  Implements decoherence from the efficient "Augmented" FSSH scheme
  outlined in Jain et al. JCTC 2016.
*/

class AFSSH: public Decoherence{
public:
  explicit AFSSH(QMInterface *qm, double dtc): Decoherence{qm, dtc} {};
  virtual ~AFSSH() {};

  bool decohere(size_t active, Electronic &c) override;
  void hop(Electronic &c) override;

private:

};

#endif
