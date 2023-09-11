#include "MonteCarlo.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include <cmath>

MonteCarlo::MonteCarlo(BlackScholesModel* mod, Option* opt, PnlRng* rng, double fdStep, long nbSamples)
        : mod_(mod), opt_(opt), rng_(rng), fdStep_(fdStep), nbSamples_(nbSamples) 
{
}
/**
 * * Calcule le prix de l'option à la date 0
*
* @param[out] prix valeur de l'estimateur Monte Carlo
* @param[out] std_dev écart type de l'estimateur
*/
void price(double& prix, double& std_dev)
    {
    
    double sum = 0.0;
    double sum_squared = 0.0;

    // Initialization de la matrice path
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (long i = 0; i < nbSamples_; ++i)
        {
            mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
            double payoff = opt_->payoff(path);
            // Pour le payoff    
            sum += payoff;
            // Pour l'ecart type
            sum_squared += payoff * payoff;
        }

    prix = exp(-mod_->r_ * opt_->T_) * 1/nbSamples_ * (sum);
        
    double variance = exp(-2 * mod_->r_ * opt_->T_) * (sum_squared / nbSamples_ - (sum / nbSamples_) * (sum / nbSamples_));
        
    std_dev = sqrt(variance / nbSamples_);
        
    pnl_mat_free(&path);
    };

