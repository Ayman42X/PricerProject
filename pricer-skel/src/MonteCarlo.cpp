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
    prix = exp(-mod_->r_ * opt_->T_) * (prix / nbSamples_);
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
    prix = exp(-mod_->r_ * (opt_->T_ - t)) * (prix / nbSamples_);
    var_M = exp(-2 * mod_->r_ * (opt_->T_ - t)) * ((var_M / nbSamples_) - prix * prix);
    std_dev = sqrt(var_M);
}

void MonteCarlo::delta(PnlVect* delta, PnlVect* std_dev){
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
        for (size_t d = 0; i < opt_->size_; i++){
            mod_->shiftAsset(path, path, d, fdStep_, 0, opt_->nbTimeSteps_);
            dfPayOff = opt_->payoff(path);
            mod_->shiftAsset(path, path, d, (-2 * fdStep_) / (1 + fdStep_), 0, opt_->nbTimeSteps_);
            dfPayOff -= opt_->payoff(path);
            mod_->shiftAsset(path, path, d, fdStep_ / (1 - fdStep_), 0, opt_->nbTimeSteps_);
            dfPayOff /= (2 * fdStep_ * MGET(path, 0, d));
            LET(delta, d) = GET(delta, d) +  dfPayOff;
            LET(std_dev, d) = GET(std_dev, d) +  dfPayOff * dfPayOff;
            dfPayOff = 0;
        }
    }
    prix = exp(-mod_->r_ * opt_->T_) * (prix / nbSamples_);
    var_price = exp(-2 * mod_->r_ * opt_->T_) * ((var_price / nbSamples_) - prix * prix);
    std = sqrt(var_price);
    for (size_t d = 0; d < opt_->size_; d++){
        LET(delta, d) = exp(-mod_->r_ * opt_->T_) * (GET(delta, d) / nbSamples_);
        LET(std_dev, d) = exp(-2 * mod_->r_ * opt_->T_) * ((GET(std_dev, d) / nbSamples_) - GET(delta, d) * GET(delta, d));
        LET(std_dev, d) = sqrt(GET(std_dev, d));
    }
    
}


