#include "BlackScholesModel.hpp"
#include "pnl/pnl_vector.h"
#include "Option.hpp"
#include "MonteCarlo.hpp"
#include "pnl/pnl_matrix.h"
#include <iostream>
#include "AsianOption.cpp"
#include "BasketOption.cpp"
#include "PerformanceOption.cpp"
#include "MonteCarlo.cpp"

int main() {
    PnlVect* sigma = pnl_vect_create_from_scalar(2,0.2);
    PnlVect* spot = pnl_vect_create_from_scalar(2,100);
    
    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);

    PnlVect* weights = pnl_vect_create_from_scalar(2,0.5);

    pnl_rng_sseed(rng, time(NULL));
    BlackScholesModel* blackScholesModel = new BlackScholesModel(2,0.02,0.0,sigma,spot);
// Option Asiatique
    // AsianOption* OptionTest = new AsianOption(1.5,150,2,100,weights);
    // MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, OptionTest,rng, 0.0001,50000);
    // double prix = 0.0;
    // double std = 0.0;
    // monteCarlo->price(prix,std);

    // PnlVect* delta = pnl_vect_create(2);
    // PnlVect* stdVect = pnl_vect_create(2);
    // monteCarlo->deltaPrice(prix,std,delta,stdVect);
    // pnl_vect_print(delta);
    // PnlMat* past = pnl_mat_create(101,2);
    // blackScholesModel->asset(past,0.5,50,rng);
    // double prix1 = 0.0;
    // double std1 = 0.0;
    // monteCarlo->price(past,0.5009,prix1,std1);
//Option Basket
    PnlVect* sigmaB = pnl_vect_create_from_scalar(40,0.2);
    PnlVect* spotB = pnl_vect_create_from_scalar(40,100);
    BlackScholesModel* BasketblackScholesModel = new BlackScholesModel(40,0.04879,0.7,sigmaB,spotB);
    PnlVect* weightsBasket = pnl_vect_create_from_scalar(40,0.025);
    BasketOption* OptionTest2 = new BasketOption(1.0,1,40,100,weightsBasket);
    MonteCarlo* monteCarloBasket = new MonteCarlo(BasketblackScholesModel, OptionTest2,rng, 0.0001,50000);
    double prix2 = 0.0;
    double std2 = 0.0;
    PnlVect* delta = pnl_vect_create(40);
    PnlVect* stdVect = pnl_vect_create(40);
    monteCarloBasket->price(prix2,std2);
    PnlMat* past = pnl_mat_create(2,40);
    BasketblackScholesModel->asset(past,1,1,rng);
    double prix1 = 0.0;
    double std1 = 0.0;
    monteCarloBasket->price(past,0.5,prix1,std1);
    monteCarloBasket->delta(delta,stdVect);
    pnl_vect_print(delta);
//Option Performance  
//    PnlVect* sigmaP =   pnl_vect_create_from_scalar(5,0.2);
//    PnlVect* spotP = pnl_vect_create_from_scalar(5,100);
//    BlackScholesModel* PerformanceblackScholesModel = new BlackScholesModel(5,0.03,0.5,sigmaP,spotP);
//    PnlVect* weightsPerformance = pnl_vect_create_from_scalar(5,0.2);
//    PerformanceOption* OptionTest3 = new PerformanceOption(2.0,12,5,weightsPerformance);
//    MonteCarlo* monteCarloPerformance = new MonteCarlo(PerformanceblackScholesModel, OptionTest3,rng, 1,50000);
//     double prix3 = 0.0;
//     double std3 = 0.0;
//     monteCarloPerformance->price(prix3,std3);
}