#ifndef _SPONGEBCH_
#define _SPONGEBCH_

#include <iostream>
#include "Macros.hpp"
#include "Options.hpp"
#include "Domain.hpp"
#include "IdealGas.hpp"
#include "BC.hpp"
#include "Utils.hpp"
#include "AbstractSingleBlockMesh.hpp"

class SpongeBC{

    public:

	AbstractSingleBlockMesh *msh;
	Domain *domain;
	IdealGas *idealGas;
	BC *bc;
	Options *opt;

	int Nx, Ny, N;

	double avgT;
	double epsP;
	double spongeP;
	double spongeStrength;
	double spongeLX0, spongeLX1;
	double spongeLY0, spongeLY1;

	double *sigma;
	double *spongeRhoAvg;
	double *spongeRhoUAvg;
	double *spongeRhoVAvg;
	double *spongeRhoEAvg;

	Options::SpongeKind spongeKind;
    
	SpongeBC(AbstractSingleBlockMesh *msh, Domain *domain, IdealGas *idealGas, BC *bc, Options *opt);
	void initRectSpongeBC();
	void initCylSpongeBC();
	void dumpSpongeAvg(int timeStep);
	void readSpongeAvgFromRestart();
};

#endif
