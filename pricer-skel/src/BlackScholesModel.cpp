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

void BlackScholesModel::price_handler(PnlVect* LastSpots,int d,double t,int myLastRaw,double timeStep,PnlMat* matriceCorrelation, PnlRng* rng){
    PnlVect* G = pnl_vect_create(d);
    PnlVect* L = pnl_vect_create(d);
    pnl_vect_rng_normal(G,d,rng);
    for(int j=0;j<d;j++){
        double sigmaShare = pnl_vect_get(this->sigma_,j);
        double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(t-myLastRaw*timeStep);
        pnl_mat_get_row(L,matriceCorrelation,j);
        quantity += sigmaShare*sqrt(t-myLastRaw*timeStep)*pnl_vect_scalar_prod(G,L);
        quantity = exp(quantity);
        double share = pnl_vect_get(LastSpots,j);
        share*=quantity;
        pnl_vect_set(LastSpots,j,share);
    }
}


void BlackScholesModel::asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past){
    double timeStep = T/nbTimeSteps;
    int d = this->spot_->size;
    pnl_mat_set_subblock(path,past,0,0);
    PnlVect* LastSpots = pnl_vect_create(d);
    int myLastRaw = static_cast<int>(floor(t*nbTimeSteps/T));
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    // Calcul des prix à S_t
    pnl_mat_get_row(LastSpots,past,myLastRaw);
    price_handler(LastSpots,d,t,myLastRaw,timeStep,matriceCorrelation,rng);
    double newStartingDate = (myLastRaw+1)*timeStep - t;
    for (int temps=myLastRaw+1;temps<nbTimeSteps+1;temps++){
        pnl_vect_rng_normal(G,d,rng);
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(newStartingDate+(temps-myLastRaw-1)*timeStep);   
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(newStartingDate+(temps-myLastRaw- 1)*timeStep)*pnl_vect_scalar_prod(G,L);
            quantity = exp(quantity);
            double share = pnl_vect_get(LastSpots,j);
            share*=quantity;
            pnl_mat_set(path,temps,j,share);
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