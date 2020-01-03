#ifndef _CABSTRACTLADH_
#define _CABSTRACTLADH_

#include "Macros.hpp"
#include "Options.hpp"
#include "AbstractFilter.hpp"
#include "AbstractCSolver.hpp"

//Abstraction class for different Locally-Added Dissipation Methods

class AbstractLAD{

    public:

	int Nx, Ny, N;
	
	AbstractCSolver *cs;

	AbstractFilter *filtX, *filtY;
	//These are for computing non-high order derivatives
	AbstractDerivatives *derivX, *derivY;
	//These are for the high-order derivatives
	AbstractDerivatives *derivX4, *derivY4;

	double *mu_star;
	double *beta_star;
	double *k_star;

	double C_mu, C_k, C_beta;

	double *S, *dil, *vort;

	double *fsw;

	double *dFbeta4dx04, *dFbeta4dx14;
	double *dFmu4dx04,   *dFmu4dx14;

	double *viz;

	virtual void calcVelocityTensorStuff(double *gradU[2][2]) = 0;
	virtual void calcLADViscosity(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE) = 0;
	virtual void calcLADBeta(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE) = 0;
	virtual void calcLADK(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE) = 0;


};

#endif
