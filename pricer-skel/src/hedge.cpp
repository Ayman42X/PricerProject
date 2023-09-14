#include "json_helper.hpp"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include "AsianOption.cpp"
#include "BasketOption.cpp"
#include "PerformanceOption.cpp"
#include "MonteCarlo.hpp"
#include "HedgingResults.hpp"
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
        pnl_vect_resize_from_scalar(weights, size, GET(volatility, 0));
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

    MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, option, rng, 0.0001, nbSamples);

   
    double initialPrice = 0.0;
    double initialPriceStdDev = 0.0;
    double finalPnl = 0.0;
    PnlVect* initDelta = pnl_vect_create(size);
    PnlVect* initDeltaStdDev = pnl_vect_create(size);

    PnlVect* delta = pnl_vect_create(size);
    PnlVect* deltaStdDev = pnl_vect_create(size);
    
    PnlMat* path = pnl_mat_create(nbTimeSteps + 1, size);
    double t = 0.0;
    
    // prix et delta à t = 0
    monteCarlo->price(initialPrice, initialPriceStdDev);
    monteCarlo->delta(initDelta, initDeltaStdDev);

    double lastPortfolioValue = initialPrice - GET(initDelta, 0) * GET(spots, 0);
    double currentPortfolioValue = 0.0;


    for (int i = 1; i <= nbTimeSteps; i++) {
        monteCarlo->delta(path, t, delta, deltaStdDev);
        currentPortfolioValue = lastPortfolioValue*exp(r*maturity/nbTimeSteps)-(GET(delta, i) - GET(delta, i-1))*GET(spots, i);
        lastPortfolioValue = currentPortfolioValue;
    }

    finalPnl = lastPortfolioValue + GET(delta, nbTimeSteps) * GET(spots, nbTimeSteps) - option->payoff(path);

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