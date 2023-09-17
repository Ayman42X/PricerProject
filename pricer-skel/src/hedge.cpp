#include "json_helper.hpp"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include "AsianOption.cpp"
#include "BasketOption.cpp"
#include "PerformanceOption.cpp"
#include "MonteCarlo.hpp"
#include "HedgingResults.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Wrong number of arguments. Exactly two arguments are required" << std::endl;
        std::exit(0);
    }
    std::ifstream ifs(argv[2]);
    nlohmann::json j = nlohmann::json::parse(ifs);
    int size;
    std::string option_type;

    PnlVect* volatility;
    PnlVect* weights;
    PnlVect* spots;

    option_type = j.at("option type").get<std::string>();

    j.at("option size").get_to(size);

    j.at("volatility").get_to(volatility);
    if (volatility->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(volatility, size, GET(volatility, 0));
    }

    j.at("spot").get_to(spots);
    if (spots->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(spots, size, GET(spots, 0));
    }

    j.at("payoff coefficients").get_to(weights);
    if (weights->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(weights, size, GET(weights, 0));
    }

    double maturity;
    j.at("maturity").get_to(maturity);

    double strike;

    double r;
    j.at("interest rate").get_to(r);

    double correlation;
    j.at("correlation").get_to(correlation);

    int nbTimeSteps;
    j.at("timestep number").get_to(nbTimeSteps);

    int nbSamples;
    j.at("sample number").get_to(nbSamples);

    int nbHedgingdates;
    j.at("hedging dates number").get_to(nbHedgingdates);

    double fdStep;
    j.at("fd step").get_to(fdStep);

    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    BlackScholesModel* blackScholesModel = new BlackScholesModel(size, r, correlation, volatility, spots);

    Option* option = nullptr;

    if (option_type == "basket") {
        j.at("strike").get_to(strike);
        option = new BasketOption(maturity, nbTimeSteps, size, strike, weights);
    } else if (option_type == "performance") {
        option = new PerformanceOption(maturity, nbTimeSteps, size, weights);
    } else if (option_type == "asian") {
        j.at("strike").get_to(strike);
        option = new AsianOption(maturity, nbTimeSteps, size, strike, weights);
    }

    MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, option, rng, fdStep, nbSamples);

   
    double initialPrice = 0.0;
    double initialPriceStdDev = 0.0;
    double finalPnl = 0.0;
    PnlVect* initDelta = pnl_vect_create(size);
    PnlVect* initDeltaStdDev = pnl_vect_create(size);

    PnlVect* delta = pnl_vect_create(size);
    PnlVect* deltaStdDev = pnl_vect_create(size);
    
    PnlMat* path = pnl_mat_create_from_file(argv[1]);
    
    // prix et delta à t = 0
    monteCarlo->deltaPrice(initialPrice, initialPriceStdDev, initDelta, initDeltaStdDev);
    double prix = 0.0;
    double stdPrix = 0.0;


    PnlVect* initSpots = pnl_vect_create(size);
    pnl_mat_get_row(initSpots, path, 0);
    double CashValue = initialPrice - pnl_vect_scalar_prod(initDelta, initSpots);
    double currentPortfolioValue = 0.0;
    PnlMat* past = pnl_mat_create_from_scalar(nbTimeSteps+1,path->n,0.0);
    PnlVect* vect = pnl_vect_create(size);
    pnl_mat_set_row(past,initSpots,0);
    
    // PnlVect* vect = pnl_vect_create(path->n);
    // pnl_mat_get_row(vect,path,1);
    // pnl_mat_set_row(past,vect,1);
    double timeStepConsta = maturity / nbTimeSteps ;
    double timeStepHedging = maturity / nbHedgingdates;
    size_t nbTimeStepHedgDansTimeStepConsta = nbHedgingdates / nbTimeSteps;
    PnlVect* diffDeltas = pnl_vect_create(size);
    for (int i = 1; i <= nbHedgingdates; i++) {
        pnl_mat_get_row(vect, path, i);
        double dateHedging = i*timeStepHedging;
        //past = pnl_mat_create_from_scalar(nbTimeSteps+1,size,0.0);
        //pnl_mat_set_row(past,initSpots,0);
        size_t temps = static_cast<size_t>(std::ceil(dateHedging / timeStepConsta));
        pnl_mat_set_row(past, vect, temps);
        for (size_t j = 0; j < temps;j++){
            //vect = pnl_vect_create(size);
            pnl_mat_get_row(vect, path, j*nbTimeStepHedgDansTimeStepConsta);
            pnl_mat_set_row(past, vect, j);
            /*size_t tempsConstatH = j*nbHedgingdates/nbTimeSteps;
            if (i<tempsConstatH){
                pnl_mat_get_row(vect,path,i);
                pnl_mat_set_row(past,vect,j);
            }
            else{
                pnl_mat_get_row(vect,path,tempsConstatH);
                pnl_mat_set_row(past,vect,j);
            }
            pnl_vect_free(&vect);*/
        }
        // PnlVect* vectDebug0 = pnl_vect_create(size);
        // pnl_mat_get_row(vectDebug0,past,0);
        // pnl_vect_print(vectDebug0);
        // PnlVect* vectDebug1 = pnl_vect_create(size);
        // pnl_mat_get_row(vectDebug1,past,1);
        // pnl_vect_print(vectDebug1);
        monteCarlo->deltaPrice(past,dateHedging, prix, stdPrix, delta, deltaStdDev);
        currentPortfolioValue = CashValue*exp((r*maturity)/nbHedgingdates);
        for (int d = 0; d < size; d++) {
            LET(diffDeltas, d) = GET(delta, d) - GET(initDelta, d);
        }
        pnl_mat_get_row(initSpots, path, i);
        currentPortfolioValue-= pnl_vect_scalar_prod(diffDeltas, initSpots);
        CashValue = currentPortfolioValue;
        initDelta = delta;
        double pnl = CashValue + pnl_vect_scalar_prod(delta, initSpots) - option->payoff(path);
        std::cout << i << std::endl;
        std::cout << prix << std::endl;
    }

    finalPnl = CashValue + pnl_vect_scalar_prod(delta, initSpots) - option->payoff(path);

    HedgingResults hedge(initialPrice, initialPriceStdDev, finalPnl);;
    std::cout << hedge << std::endl;


    pnl_vect_free(&initDelta);
    pnl_vect_free(&initDeltaStdDev);
    pnl_vect_free(&delta);
    pnl_vect_free(&deltaStdDev);
    pnl_vect_free(&spots);
    pnl_vect_free(&weights);
    pnl_vect_free(&volatility);
    pnl_rng_free(&rng);
    exit(0);
}