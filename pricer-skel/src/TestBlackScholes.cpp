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
    //BlackScholesModel* blackScholesModel = new BlackScholesModel(2,0.02,0.0,sigma,spot);
// Option Asiatique
    // AsianOption* OptionTest = new AsianOption(1.5,150,2,100,weights);
    // MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, OptionTest,rng, 0.0001,50000);
    // double prix = 0.0;
    // double std = 0.0;

    // PnlVect* delta = pnl_vect_create(2);
    // PnlVect* stdVect = pnl_vect_create(2);
    // monteCarlo->deltaPrice(prix,std,delta,stdVect);
    // pnl_vect_print(delta);
    // PnlMat* past = pnl_mat_create(151,2);
    // pnl_mat_set_row(past,spot,0);
    // blackScholesModel->asset(past,0.5,50,rng);
    // double prix1 = 0.0;
    // double std1 = 0.0;
    // PnlVect* delta1 = pnl_vect_create(2);
    // PnlVect* stdDelta1 = pnl_vect_create(2);
    // monteCarlo->deltaPrice(past,0.100000001,prix1,std1,delta1,stdDelta1);
//Option Basket
    PnlVect* sigmaB = pnl_vect_create_from_scalar(40,0.2);
    PnlVect* spotB = pnl_vect_create_from_scalar(40,100);
    BlackScholesModel* BasketblackScholesModel = new BlackScholesModel(40,0.04879,0.7,sigmaB,spotB);
    PnlVect* weightsBasket = pnl_vect_create_from_scalar(40,0.025);
    BasketOption* OptionTest2 = new BasketOption(0.9973,1,40,100,weightsBasket);
    MonteCarlo* monteCarloBasket = new MonteCarlo(BasketblackScholesModel, OptionTest2,rng, 0.0001,50000);
    PnlVect* delta = pnl_vect_create(40);
    PnlVect* stdVect = pnl_vect_create(40);
    PnlMat* past = pnl_mat_create(2,40);
    pnl_mat_set_row(past,spotB,0);
    //pnl_mat_set_row(past,spotB,1);
    double prix1 = 0.0;
    double std1 = 0.0;
    std::vector<double> valeurs_reelles = {
        100.697076, 100.428997, 101.067177, 101.570036, 101.449229,
        101.159191, 100.136511, 101.842047, 101.293857, 101.486780,
        101.759902, 101.487399, 101.623716, 100.309059, 100.975381,
        101.061105, 102.049386, 101.352772, 101.056671, 100.728371,
        100.941014, 100.761932, 101.766708, 100.810520, 100.735441,
        101.559637, 101.770946, 101.005296, 101.698503, 101.456995,
        101.760860, 101.174944, 102.054293, 101.030278, 99.927319,
        99.661293, 101.436243, 100.881928, 101.199387, 101.247588
        };
    for (size_t i = 0; i < valeurs_reelles.size(); ++i) {
        MLET(past,0,i) = valeurs_reelles[i];
    }
monteCarloBasket->deltaPrice(prix1,std1,delta,stdVect);
    //monteCarloBasket->deltaPrice(past,0.0027,prix1,std1,delta,stdVect);
    pnl_vect_print(delta);
//Option Performance  
//    PnlVect* sigmaP =   pnl_vect_create_from_scalar(5,0.2);
//    PnlVect* spotP = pnl_vect_create_from_scalar(5,100);
//    BlackScholesModel* PerformanceblackScholesModel = new BlackScholesModel(5,0.03,0.5,sigmaP,spotP);
//    PnlVect* weightsPerformance = pnl_vect_create_from_scalar(5,0.2);
//    PerformanceOption* OptionTest3 = new PerformanceOption(2.0,12,5,weightsPerformance);
//    MonteCarlo* monteCarloPerformance = new MonteCarlo(PerformanceblackScholesModel, OptionTest3,rng, 0.0001,50000);
//     double prix3 = 0.0;
//     double std3 = 0.0;
//     PnlVect* delta3 = pnl_vect_create(5);
//     PnlVect* stdDelta3 = pnl_vect_create(5);
//     monteCarloPerformance->deltaPrice(prix3,std3,delta3,stdDelta3);
//     pnl_vect_print(delta3);
}