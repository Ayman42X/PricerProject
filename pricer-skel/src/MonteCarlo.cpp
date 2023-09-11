#include "MonteCarlo.hpp"
MonteCarlo::MonteCarlo(BlackScholesModel* mod, Option* Option, PnlRng* rng, double& fdStep, long& nbSamples){
    mod_ = mod, opt_ = Option, rng_ = rng, fdStep_ = fdStep, nbSamples_ = nbSamples;
}

void MonteCarlo::price(double& prix, double& std_dev){
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    double var_M = 0;
    for (size_t i = 0; i < nbSamples_; i++)
    {
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        double payOff = opt_->payoff(path);
        prix += payOff;
        var_M += payOff * payOff;
    }
    prix /= nbSamples_;
    var_M = exp(-2 * mod_->r_ * opt_->T_) * ((var_M / nbSamples_) - prix * prix);
    std_dev = sqrt(var_M);
    
}

void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& std_dev){
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    double var_M = 0;
    for (size_t i = 0; i < nbSamples_; i++)
    {
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        double payOff = opt_->payoff(path);
        prix += payOff;
        var_M += payOff * payOff;
    }
    prix /= nbSamples_;
    var_M = exp(-2 * mod_->r_ * opt_->T_) * ((var_M / nbSamples_) - prix * prix);
    std_dev = sqrt(var_M);
}

