#ifndef _CABSTRACTFILTERH_
#define _CABSTRACTFILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "AbstractDerivatives.hpp"
#include "Options.hpp"
#include "Macros.hpp"
#include "Domain.hpp"
#include "Utils.hpp"
#include "BC.hpp"

class AbstractFilter{

    public:

	int Nx, Ny, N;

	double *diagF, *offlowerF, *offupperF;

	double *offlowerF2, *offupperF2;

	AbstractDerivatives::Direct currentDir;
	Options::BCType bcType;

	BC *bc;

	//Functions that a filter object have to define
	virtual void filterField(double *dataIn, double *dataOut) = 0;

};

#endif
