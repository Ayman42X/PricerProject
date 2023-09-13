#include "Option.hpp"
// #pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <algorithm>
#include <iostream>

class EuropeanCallOption : public Option
{
public:
    double strike_; 
    EuropeanCallOption(double T, int nbTimeSteps, int size, double strike)
        : Option(T, nbTimeSteps, size), strike_(strike)
    {
    }
    double payoff(const PnlMat* path) override
    {
        double ST = pnl_mat_get(path, path->m - 1, 0); //sous jacent a la maturité
        return std::max(0.0, ST - strike_);
    }
};



// int main() {
//     double T = 1.0;
//     int nbTimeSteps = 100; 
//     int size = 3; 
//     double strike = 100.0; 
//     PnlVect* weights = pnl_vect_create_from_scalar(size, 1.0 / size); 
//     PnlMat* path = pnl_mat_create_from_double(nbTimeSteps + 1, size, 125.0); // Trajectoire constante à 100

//     BasketOption basketOption(T, nbTimeSteps, size, strike, weights);
//     double basketPayoff = basketOption.payoff(path);
//     std::cout << "Payoff de l'option panier : " << basketPayoff << std::endl;

//     PerformanceOption performanceOption(T, nbTimeSteps, size, weights);
//     double performancePayoff = performanceOption.payoff(path);
//     std::cout << "Payoff de l'option performance sur panier : " << performancePayoff << std::endl;
//     pnl_mat_free(&path);
//     pnl_vect_free(&weights);
//     return 0;
// }