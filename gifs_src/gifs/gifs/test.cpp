#include "stddef.h"
//
#include "constants.hpp"
//
#include "properties.hpp"
#include "qm.hpp"
#include "bomd.hpp"
#include "gifs.hpp"
//
#include "gifs_base.hpp"
//
#include <string>
#include <vector>
#include <iostream>


void
do_something()
{
    Gifs gifs_handle{};
}


int
main()
{
    std::vector<int> qmid = {6, 1, 1, 1, 1};
    Gifs gifs_handle(qmid.size(), qmid);

    std::vector<double> crd{};
    std::vector<double> grad{};

    double* null = nullptr; 

    crd.resize(3 * qmid.size());
    grad.resize(3 * qmid.size(), 10);

    std::cout << grad[0] << "\n";

    gifs_handle.get_gradient(crd.data(), 0, null, null, grad.data(), grad.data());

    std::cout << grad[0] << "\n";

    return 0;
}
