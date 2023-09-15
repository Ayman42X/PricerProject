#include "BlackScholesModel.hpp"
#include <cmath>
#include <iostream>

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
    pnl_mat_set_row(path,this->spot_,0);
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    //pnl_mat_print(matriceCorrelation);
    pnl_mat_chol(matriceCorrelation);
    //pnl_mat_print(matriceCorrelation);
    for (int t=1;t<nbTimeSteps+1;t++){
        pnl_vect_rng_normal(G,d,rng);
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*timeStep;
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(L,G);
            quantity = exp(quantity);
            double share = MGET(path,t-1,j);
            share*=quantity;
            MLET(path,t,j) = share; 
        } 
    }
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}


void BlackScholesModel::asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past){
    // Initialisation de variables utilisées
    double timeStep = T/nbTimeSteps;
    int d = this->spot_->size;
    int iPlus1 = static_cast<int>(std::ceil(t/timeStep));
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlVect* s_t = pnl_vect_create(d);
    pnl_mat_set_subblock(path,past,0,0);
    // Calcul de la matrice de correlation
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    pnl_mat_get_row(s_t,past,iPlus1);
    double newStartingDate = (iPlus1)*timeStep - t;
    for (int temps=iPlus1;temps<nbTimeSteps+1;temps++){
        pnl_vect_rng_normal(G,d,rng);
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(timeStep);
            if (temps == iPlus1){
                quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(newStartingDate);
            }   
            pnl_mat_get_row(L,matriceCorrelation,j);
            if (temps == iPlus1){
                quantity += sigmaShare*sqrt(newStartingDate)*pnl_vect_scalar_prod(G,L);
            }
            else{
                quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(G,L);
            }
            quantity = exp(quantity);
            double share = MGET(path,temps-1,j);
            if (temps == iPlus1){
                share = MGET(path,temps,j);
            }
            share*=quantity;
            MLET(path,temps,j) = share;
        }
    }
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}


void BlackScholesModel::shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep){
    size_t timeTRowIndex = static_cast<size_t>(std::ceil(t / timestep));
    for (size_t i = timeTRowIndex; i < path->m; i++){
        MLET(shift_path, i, d) = MGET(shift_path, i, d) * (1 + h);
    }
}

void BlackScholesModel::simul_market(PnlMat* path, double T, int H, PnlRng* rng){
    double timeStep = T/H;
    int d = this->spot_->size;
    pnl_mat_set_row(path,this->spot_,0);
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    //pnl_mat_print(matriceCorrelation);
    pnl_mat_chol(matriceCorrelation);
    //pnl_mat_print(matriceCorrelation);
    for (int t=1;t<H+1;t++){
        pnl_vect_rng_normal(G,d,rng);
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double trendShare = GET(this->trend_,j);
            double quantity = (trendShare- (sigmaShare*sigmaShare)/2)*timeStep;
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(L,G);
            quantity = exp(quantity);
            double share = MGET(path,t-1,j);
            share*=quantity;
            MLET(path,t,j) = share; 
        } 
    }
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}