#include "BlackScholesModel.hpp"
#include <cmath>
#include <iostream>

BlackScholesModel::BlackScholesModel(int size, double r, double rho,PnlVect* sigma,PnlVect* spot){
    this->size_ = size;
    this->r_ = r;
    this->rho_ = rho;
    this->sigma_ = sigma;
    this->spot_ = spot;
    this->iPlus1 = 0;
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
    size_t iPlus1 = past->m -1;
    this->iPlus1 = past->m - 1;
    PnlVect* G = pnl_vect_new();
    PnlVect* L = pnl_vect_new();
    pnl_mat_set_subblock(path,past,0,0);
    // Calcul de la matrice de correlation
    PnlMat* matriceCorrelation = pnl_mat_create_from_scalar(d,d,this->rho_);
    for (int k = 0;k<d;k++){
        pnl_mat_set(matriceCorrelation,k,k,1.0);
    }
    pnl_mat_chol(matriceCorrelation);
    double newStartingDate = (iPlus1)*timeStep - t;
    for (int temps=iPlus1;temps<nbTimeSteps+1;temps++){
        pnl_vect_rng_normal(G,d,rng);
        for (int j=0;j<d;j++){
            double sigmaShare = pnl_vect_get(this->sigma_,j);
            pnl_mat_get_row(L,matriceCorrelation,j);
            double quantity = 0.0;
            double share = 0.0;
            if (temps == iPlus1){
                quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(newStartingDate);
                quantity += sigmaShare*sqrt(newStartingDate)*pnl_vect_scalar_prod(G,L);
                quantity = exp(quantity);
                share = MGET(path,temps,j);
            }
            else{
                quantity = (this->r_ - (sigmaShare*sigmaShare)/2)*(timeStep);
                quantity += sigmaShare*sqrt(timeStep)*pnl_vect_scalar_prod(G,L);
                quantity = exp(quantity);
                share = MGET(path,temps-1,j);
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
    size_t timeTRowIndex = this->iPlus1;
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

size_t BlackScholesModel::handler_time(double t, double timeStep){
    double a = t/timeStep;
    size_t time = static_cast<size_t>(std::ceil(a));
    double precision = 1e-8;
    if (time - (a) < 0){
        return time +1;
    }
    if (time - (t/timeStep) > 0){
        if (time - a > precision){
            if (time -a > (timeStep/3) - precision){
                return time;
            }
            else{
               return time -1;
            }

        }
    }
    if (time == t/timeStep){
        return time;
    }
    return time;
}