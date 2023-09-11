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

class BasketOption : public Option
{
public:
    double strike_; 
    PnlVect* weights_; //vect poids
    
    BasketOption(double T, int nbTimeSteps, int size, double strike, const PnlVect* weights)
        : Option(T, nbTimeSteps, size), strike_(strike)
    {
        weights_ = pnl_vect_copy(weights); 
    }

    //destructeur    
    ~BasketOption()
    {
        pnl_vect_free(&weights_);
    }

    //payoff
    double payoff(const PnlMat* path) override
    {
        double weightedSum = 0.0;
        for (int d = 0; d < size_; ++d)
        {
            double assetPrice = pnl_mat_get(path, path->m - 1, d);
            weightedSum += pnl_vect_get(weights_, d) * assetPrice;
        }

        double payoff = weightedSum - strike_;
        return std::max(0.0, payoff); 
    }
};

//CLASSE OPTION PERFORMANCE
class PerformanceOption : public Option
{
    public:
    PnlVect* weights_; //vect poids
    PerformanceOption(double T, int nbTimeSteps, int size, const PnlVect* weights)
        : Option(T, nbTimeSteps, size), weights_(pnl_vect_copy(weights))
    {
    }
    ~PerformanceOption()
    {
        pnl_vect_free(&weights_);
    }
    double payoff(const PnlMat* path) override
    {

        double performance = 1.0;
        for (int i=1; i<= nbTimeSteps_; i++)
        {
            //somme des prix a l'intant ti
            double sum_ti = 0.0;
            double sum_ti_minus_1 = 0.0; //Somme  à l'instant ti-1
        
            for(int d= 0; d < size_; d++)
            {
                //somme pondérees des 2 instants
                sum_ti += pnl_mat_get(path, i, d) * pnl_vect_get(weights_, d);
                sum_ti_minus_1 += pnl_mat_get(path, i-1, d) * pnl_vect_get(weights_, d);
            }
            //terme de performance a l'intant ti, c'ets le rendement relatif entre les 2 instants
            double term = (sum_ti / sum_ti_minus_1) - 1.0;
            term = (term > 0.0) ? term : 0.0;
            //maj de la performance 
            performance *= (1.0 + term);   
        }
        return performance - 1.0;   
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