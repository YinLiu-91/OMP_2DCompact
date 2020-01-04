#ifndef _E4ODH_
#define _E4ODH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"

class Explicit4thOrder: public AbstractDerivatives{


    public:

        double a, b;
        double a11, a21, a31, a41, a51, a61;
        double a12, a22, a32, a42, a52, a62;
        double a13, a23, a33, a43, a53;

	//Constructor
        Explicit4thOrder(Domain *dom, Options::BCType bcType, Direct currentDir){

   	    this->Nx = dom->Nx;
	    this->Ny = dom->Ny;
	    this->dx = dom->dx; 
	    this->dy = dom->dy;

	    this->currentDir = currentDir;

	    this->bcType = bcType;

	    if(currentDir == DIRX){
		N  = Nx;
		dd = dx;
	    }else if(currentDir == DIRY){
		N  = Ny;
		dd = dy;
	    }

	    //Coefficients for the interior 4th order, 4th derivative
	    a = 2.0, b = -1.0;

	    //2nd-order boundary scheme for fourth derivative
	    //boundary point 1
	    a11 =   3.0;
	    a21 = -14.0;
	    a31 =  26.0;
 	    a41 = -24.0;
	    a51 =  11.0;
	    a61 =  -2.0;

	    //boundary point 2
	    a12 =   2.0;
	    a22 =  -9.0;
	    a32 =  16.0;
	    a42 = -14.0;
	    a52 =   6.0;
	    a62 =  -1.0;

	    //2nd order central scheme for 4th derivative at boundary point 3
	    a13 =   1.0;
	    a23 =  -4.0;
	    a33 =   6.0;
	    a43 =  -4.0;
	    a53 =   1.0;

    }

    //Function's to call...
    void calc1stDerivField(double *dataIn, double *dataOut);
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm2, double *Nm1, double *Np1, double *Np2){};
    void calc1stDerivField_TPB(double *dataIn, double *dataOut, double *Nm3, double *Nm2, double *Nm1, double *Np1, double *Np2, double *Np3){};//empty call here since our bandwidth is 5

    void calc1stDeriv(double *phi, double *dphi);

    double calcNeumann(double *f){return 0.0;};

    void Calc4thPeriodic(double *phi, double *dphidx);
    void Calc4thDirichlet(double *phi, double *dphidx);
};

#endif
