#include "MonteCarlo.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
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
void MonteCarlo::price(double& prix, double& std_dev)
    {
    
    double sum = 0.0;
    double sum_squared = 0.0;

    // Initialization de la matrice path
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (long i = 0; i < nbSamples_; i++)
        {
            mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
            double payoff = opt_->payoff(path);
            // Pour le payoff    
            sum += payoff;
            // Pour l'ecart type
            sum_squared += payoff * payoff;
        }

    prix = exp(-1*mod_->r_ * opt_->T_)  * (sum/ nbSamples_);
        
    double variance = exp(-2 * mod_->r_ * opt_->T_) * (sum_squared / nbSamples_ - (sum / nbSamples_) * (sum / nbSamples_));
        
    std_dev = sqrt(variance / nbSamples_);
        
    pnl_mat_free(&path);
    };

     
    /**
     * Calcule le prix de l'option à la date t
     *
     * @param[in]  past contient la trajectoire du sous-jacent
     * jusqu'à l'instant t
     * @param[in] t date à laquelle le calcul est fait
     * @param[out] prix contient le prix
     * @param[out] std_dev contient l'écart type de l'estimateur
     */
   void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& std_dev)
{


     // Initialization de la matrice path
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    // Initialisation des variables pour stocker les prix simulés
    double sum = 0.0;
    double sum_squared = 0.0;

    for (long i = 0; i < nbSamples_; ++i)
    {
        // Génération d'une trajectoire du modèle à partir de t
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        double payoff = opt_->payoff(path);
        sum += payoff;

        sum_squared += payoff * payoff;
    }

    prix = exp(-mod_->r_ * (opt_->T_ - t)) * (sum / nbSamples_);
    // A revoir la formule de la variance
    double variance = exp(-2 * mod_->r_ * (opt_->T_ - t)) * (sum_squared / nbSamples_  - (sum / nbSamples_) * (sum / nbSamples_));
    std_dev =  sqrt(variance / nbSamples_);

    pnl_mat_free(&path);
}


    void MonteCarlo::delta(PnlVect* delta, PnlVect* std_dev) {
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    PnlMat* shift_path1 = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    PnlMat* shift_path2 = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    for (size_t d = 0; d < opt_->size_; d++)
    {
        double dfEstimator = 0;
        for (size_t i = 0; i < nbSamples_; i++)
        {
            mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
            pnl_mat_clone(shift_path1, path);
            pnl_mat_clone(shift_path2, path);
            mod_->shiftAsset(shift_path1, path, d, fdStep_, 0, opt_->nbTimeSteps_);
            mod_->shiftAsset(shift_path2, path, d, -fdStep_, 0, opt_->nbTimeSteps_);
            double dfPayOff = (opt_->payoff(shift_path1) - opt_->payoff(shift_path2)) / (2 * fdStep_ * MGET(path, 0, d));
            LET(delta, d) = GET(delta, d) +  dfPayOff;
            dfEstimator += dfPayOff * dfPayOff;
        }
        LET(delta, d) = exp(-mod_->r_ * opt_->T_) * (GET(delta, d) / nbSamples_);
        dfEstimator = exp(-2 * mod_->r_ * opt_->T_) * ((dfEstimator / nbSamples_) - GET(delta, d) * GET(delta, d));
        LET(std_dev, d) = sqrt(dfEstimator);
    }   
}



void MonteCarlo::deltaPrice(double& prix, double& std, PnlVect* delta, PnlVect* std_dev){
          PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    //PnlMat* shift_path1 = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    //PnlMat* shift_path2 = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    double var_price = 0;
    double dfPayOff = 0;
    for (size_t i = 0; i < nbSamples_; i++){
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        double payOff = opt_->payoff(path);
        prix += payOff;
        var_price += payOff * payOff;
        //pnl_mat_clone(shift_path1, path);
        //pnl_mat_clone(shift_path2, path);
        for (size_t d = 0; d < opt_->size_; d++){
            //pnl_mat_print(path);
            mod_->shiftAsset(path, path, d, fdStep_, 0, opt_->nbTimeSteps_);
            //pnl_mat_print(path);
            dfPayOff = opt_->payoff(path);
            mod_->shiftAsset(path, path, d, (-2 * fdStep_) / (1 + fdStep_), 0, opt_->nbTimeSteps_);
            //pnl_mat_print(path);
            dfPayOff -= opt_->payoff(path);
            mod_->shiftAsset(path, path, d, fdStep_ / (1 - fdStep_), 0, opt_->nbTimeSteps_);
            //pnl_mat_print(path);
            dfPayOff /= (2 * fdStep_ * MGET(path, 0, d));
            LET(delta, d) = GET(delta, d) +  dfPayOff;
            LET(std_dev, d) = GET(std_dev, d) +  dfPayOff * dfPayOff;
            dfPayOff = 0;
        }
    }
    prix = (prix / nbSamples_);
    var_price = exp(-2 * mod_->r_ * opt_->T_) * ((var_price / nbSamples_) - prix * prix);
    prix *= exp(-mod_->r_ * opt_->T_);
    std = sqrt(var_price/nbSamples_);
    for (size_t d = 0; d < opt_->size_; d++){
        LET(delta, d) = exp(-mod_->r_ * opt_->T_) * (GET(delta, d) / nbSamples_);
        LET(std_dev, d) = exp(-2 * mod_->r_ * opt_->T_) * ((GET(std_dev, d) / nbSamples_) - GET(delta, d) * GET(delta, d));
        LET(std_dev, d) = sqrt(GET(std_dev, d)/nbSamples_);
    }
    
}