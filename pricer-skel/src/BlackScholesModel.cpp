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
    double timeStep = T/nbTimeSteps;
    int d = this->spot_->size;
    pnl_mat_set_subblock(path,past,0,0);
    PnlVect* LastSpots = pnl_vect_new();
    // Algo pour récupérer le dernier indice dans le past
    double myLastRaw = past->m-1;
    double myLastRawValue;
    for(int indice = 1; indice < past->m; indice++){
        myLastRawValue = MGET(past,indice,0);
        if(myLastRawValue == 0){
            myLastRaw = indice - 1;
            break;
        }
    }
    pnl_mat_get_row(LastSpots,past,myLastRaw);
    double x = MGET(past,myLastRaw,0);
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    //pnl_mat_print(path);
    double newStartingDate = (myLastRaw+1)*timeStep - t;
    for (int temps=myLastRaw+1;temps<nbTimeSteps+1;temps++){
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            double quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(newStartingDate+(temps-myLastRaw-1)*timeStep);
            pnl_vect_rng_normal(G,d,rng);
            pnl_mat_get_row(L,matriceCorrelation,j);
            quantity += sigmaShare*sqrt(newStartingDate+(temps-myLastRaw- 1)*timeStep)*pnl_vect_scalar_prod(G,L);
            quantity = exp(quantity);
            double share = pnl_vect_get(LastSpots,j);
            share*=quantity;
            pnl_mat_set(path,temps,j,share);
        }
        //std::cout << "Nouvelle Matrice" << std::endl;
        //pnl_mat_print(path);
    }
    double xxx = MGET(past,0,0);
    double xx = MGET(path, 0, 0);
    pnl_vect_free(&G);
    pnl_vect_free(&L);
    pnl_mat_free(&matriceCorrelation);
}


void BlackScholesModel::shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep){
    for (size_t i = (size_t) (t / timestep); i < path->m; i++){
        MLET(shift_path, i, d) = MGET(shift_path, i, d) * (1 + h);
    }
}