#ifndef GIFS_QM_NULL_HPP
#define GIFS_QM_NULL_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include "configreader.hpp"
#include <armadillo>


class QM_Null: public QMInterface{
public:
  QM_Null(std::shared_ptr<INAQSShared> shared,
          FileHandle& fh,
          const arma::uvec& in_qmids,
          arma::mat& in_qm_crd,
          arma::mat& in_mm_crd,
          arma::vec& in_mm_chg,
          const int charge,
          const int mult,
          const int excited_states,
          const int min_state):
    QMInterface(shared, in_qmids, in_qm_crd, in_mm_crd, in_mm_chg, charge, mult, excited_states, min_state)
  {(void) fh;};

  void get_properties(PropMap &props) override {(void) props;};
};

#endif
