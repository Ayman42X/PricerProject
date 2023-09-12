#include "BlackScholesModel.hpp"
#include "pnl/pnl_vector.h"
#include "Option.hpp"
#include "MonteCarlo.hpp"
#include "pnl/pnl_matrix.h"
#include <iostream>
#include "AsianOption.cpp"
#include "Option.cpp"
#include "MonteCarlo.cpp"

int main() {
    std::cout << "Prix des option : " << std::endl;
    PnlVect* sigma = pnl_vect_create_from_scalar(2,0.2);
    PnlVect* spot = pnl_vect_create_from_scalar(2,100);
    PnlMat* path = pnl_mat_create(11,2);
    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);

    PnlVect* weights = pnl_vect_create_from_scalar(2,0.5);

    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel* blackScholesModel = new BlackScholesModel(2,0.04,0.0,sigma,spot);
    //blackScholesModel->asset(path,5,10,rng);

    AsianOption* OptionTest = new AsianOption(5,10,2,100,weights);
    MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, OptionTest,rng, 1,1000000);

    BasketOption* OptionTest2 = new BasketOption(5,10,2,100,weights);
    MonteCarlo* monteCarlo2 = new MonteCarlo(blackScholesModel, OptionTest2,rng, 1,1000000);

    PerformanceOption* OptionTest3 = new PerformanceOption(5,10,2,weights);
    MonteCarlo* monteCarlo3 = new MonteCarlo(blackScholesModel, OptionTest3,rng, 1,1000000);

    double prix = 0.0;
    double std = 0.0;

    double prix2 = 0.0;
    double std2 = 0.0;

    double prix3 = 0.0;
    double std3 = 0.0;
    
    monteCarlo->price(prix,std);
    monteCarlo2->price(prix2,std2);
    monteCarlo3->price(prix3,std3);

    std::cout << "Le prix à l'intant 0 de l'option ASIATIQUE est de  : " << prix << std::endl;
    std::cout << "L'ecart type OPTION ASIATIQUE : " << std << std::endl;

    std::cout << "Le prix à l'intant 0 de l'option BASKET est de  : " << prix2 << std::endl;
    std::cout << "L'ecart type OPTION BASKET : " << std2 << std::endl;

    std::cout << "Le prix à l'intant 0 de l'option PERFORMANCE est de  : " << prix3 << std::endl;
    std::cout << "L'écart type OPTION PERFORMANCE : " << std3 << std::endl;
}

