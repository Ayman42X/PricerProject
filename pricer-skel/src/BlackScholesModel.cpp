#include "BlackScholesModel.hpp"
#include <cmath>

BlackScholesModel::BlackScholesModel(int size, double r, double rho,PnlVect* sigma,PnlVect* spot){
    this->size_ = size;
    this->r_ = r;
    this->rho_ = rho;
    this->sigma_ = sigma;
    this->spot_ = spot;
}

void BlackScholesModel::asset(PnlMat* path, double T, int nbTimeSteps, PnlRng* rng){
    double timeStep = T/nbTimeSteps;
    int d = this->spot_->size;
    path = pnl_mat_create(nbTimeSteps+1,d);
    pnl_mat_set_col(path,this->spot_,0);
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(nbTimeSteps+1,d,this->rho_);
    for (int k = 0;k<nbTimeSteps+1;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    for (int t=1;t<nbTimeSteps+1;t++){
        for (int j=1;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*timeStep;
            pnl_vect_rng_normal(G,d,rng);
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(G,L);
            quantity = exp(quantity);
            double share = pnl_mat_get(path,t-1,j);
            share*=quantity;
            pnl_mat_set(path,t,j,share);
        } 
    }
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}

void BlackScholesModel::asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past){
    double timeStep = T/nbTimeSteps;
    int d = this->spot_->size;
    path = pnl_mat_create(nbTimeSteps+1,d);
    pnl_mat_set_subblock(path,past,0,0);
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(nbTimeSteps+1,d,this->rho_);
    for (int k = 0;k<nbTimeSteps+1;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    for (int t=1;t<nbTimeSteps+1;t++){
        for (int j=1;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*timeStep;
            pnl_vect_rng_normal(G,d,rng);
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(G,L);
            quantity = exp(quantity);
            double share = pnl_mat_get(path,t-1,j);
            share*=quantity;
            pnl_mat_set(path,t,j,share);
        } 
    }
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}

void BlackScholesModel::shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep){
    for (size_t i = (size_t) (t / timestep); i < path->m; i++){
        MLET(shift_path, i, d) = MGET(shift_path, i, d) * (1 + h);
    }
    
}


