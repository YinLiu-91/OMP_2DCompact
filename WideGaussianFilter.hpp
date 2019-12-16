#ifndef _WIDEGAUSSFILTERH_
#define _WIDEGAUSSFILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractFilter.hpp"

//Wide Gaussian Filter taken from "Assessment of high-resolution methods for numerical simulations of compressible turbulence with shock waves", E. Johnsen et al. (2010)

class WideGaussianFilter: public AbstractFilter{

    public:

    //Interior scheme
    double a0, a1, a2, a3, a4;
    
    //NOTE: For non-periodic boundary conditions just going to do symmetry, probably should have like a compact boundary scheme here instead

    WideGaussianFilter(Domain *dom, BC *bc, Options::BCType bcType, AbstractDerivatives::Direct currentDir){

	this->Nx = dom->Nx;
	this->Ny = dom->Ny;

	this->bc = bc;
	this->currentDir = currentDir;
	this->bcType = bcType;

	if(currentDir == AbstractDerivatives::DIRX){
	    N = Nx;
	}else if(currentDir == AbstractDerivatives::DIRY){
	    N = Ny;
	}		

	a0 = 3565.0/10368.0;
	a1 = 3091.0/12960.0;
	a2 = 1997.0/25920.0;
	a3 =  149.0/12960.0;
	a4 =  107.0/103680.0;	

    }

    void filterPeriodic(double  *phi, double *phiF);
    void filterDirichlet(double *phi, double *phiF);

    void filter(double *phi, double *phiF);
    void filterField(double *dataIn, double *dataOut);

};

#endif
