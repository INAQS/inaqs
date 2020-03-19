#include "linkatoms.hpp"


LinkAtoms*
LinkAtoms::with_const_factors(arma::umat global_idx, arma::vec factors)
{
    return new LinkAtoms(global_idx, factors, arma::vec{});
};
//
LinkAtoms* 
LinkAtoms::with_const_length(arma::umat global_idx, arma::vec dist)
{
    arma::vec fac(dist.size());
    return new LinkAtoms(global_idx, fac, dist);
};
