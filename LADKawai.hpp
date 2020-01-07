#ifndef _CLADKAWAILELEH_
#define _CLADKAWAILELEH_

#include "Macros.hpp"
#include "Options.hpp"
#include "Utils.hpp"
#include "AbstractLAD.hpp"
#include "AbstractFilter.hpp"
#include "AbstractDerivatives.hpp"
#include "WideGaussianFilter.hpp"
#include "Explicit4thOrder.hpp"

//Locally-Added Dissipation methods of Kawai, Shankar, & Lele (2010)

class LADKawai: public AbstractLAD{

  public:

    LADKawai(AbstractCSolver *cs){

	this->cs = cs;
	this->Nx = cs->Nx;
	this->Ny = cs->Ny;
	this->N  = Nx*Ny;

	//final outputs to the solver
	mu_star   = new double[N]{0.0};
	beta_star = new double[N]{0.0};
	k_star    = new double[N]{0.0};

	//scalar values pertaining to the velocity tensor
	S    = new double[N];
	dil  = new double[N];
	vort = new double[N];
   
	//Model coefficients 
	C_mu   = 0.002;
	C_k    = 0.01;
	C_beta = 1.75;

	//shock sensor...
	fsw = new double[N];

	dFbeta4dx04 = new double[N];
	dFbeta4dx14 = new double[N];
	dFmu4dx04 = new double[N];
	dFmu4dx14 = new double[N];

	viz = new double[N];
	
	//Our wide gaussian filter...
	filtX = new WideGaussianFilter(cs->dom, cs->bc, cs->bc->bcXType, AbstractDerivatives::DIRX);
	filtY = new WideGaussianFilter(cs->dom, cs->bc, cs->bc->bcYType, AbstractDerivatives::DIRY);

	//For now just use whatever the solver is using
	derivX = new Explicit4thOrder(cs->dom, cs->bc->bcXType, AbstractDerivatives::DIRX);
	derivY = new Explicit4thOrder(cs->dom, cs->bc->bcYType, AbstractDerivatives::DIRY);

    }

    void calcVelocityTensorStuff(double *gradU[2][2]);
    void calc4thOrderDerivative(double *phi, double *d4phi0, double *d4phi1, double *work1, double *work2);
    void calcLADViscosity(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE);
    void calcLADBeta(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE);
    void calcLADK(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE);

};

#endif
