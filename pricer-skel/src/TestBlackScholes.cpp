#include "BlackScholesModel.hpp"
#include "pnl/pnl_vector.h"
#include <iostream>
int main() {
    std::cout << "Hello, World!" << std::endl;
    PnlVect* sigma = pnl_vect_create_from_scalar(2,0.2);
    PnlVect* spot = pnl_vect_create_from_scalar(2,100);
    PnlMat* path = pnl_mat_create(2,11);
    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel* blackScholesModel = new BlackScholesModel(2,0.04,0.0,sigma,spot);
    blackScholesModel->asset(path,5,10,rng);
}