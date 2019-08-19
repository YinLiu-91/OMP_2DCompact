#ifndef _CABSTRACTSINGLEBLOCKMESHH_
#define _CABSTRACTSINGLEBLOCKMESHH_

#include <iostream>
#include <fstream>

#include "Options.hpp"
#include "AbstractCSolver.hpp"
#include "Domain.hpp"
#include "AbstractDerivatives.hpp"
#include "Adt.hpp"

//There's some cyclic dependency going on here...
class AbstractCSolver;

class AbstractSingleBlockMesh{

    public:

	AbstractCSolver *cs;
	Domain *d;

	int mpiRank;
	double *x, *y, *z;

	double *J;
	double *J11, *J12;
	double *J21, *J22;

	double x_min, x_max;
	double y_min, y_max;
	double max_xi, max_eta;
	int Nx, Ny;

        double periodicXTranslation[2];
        double periodicYTranslation[2];

	bool periodicBCX, periodicBCY;
	bool transPeriodicX, transPeriodicY;
	bool interPeriodicX, interPeriodicY;

	AbstractDerivatives *derivX, *derivY;

	Adt<double> *adt;

	void allocateForMesh();
	void solveForJacobians();
	void dumpGrid();
	void getOrderedBlockCoordinates(int ip, int jp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]);
	void getOrderedBlockXiCoordinates(int ip, int jp, double box_pxi[8][3]);
	void generateCoordinateHaloArrays(double *&x_halo, double *&y_halo);
	int findCVForPoint(double p[3], double *x_halo, double *y_halo);
	void initMeshADT();
	void handlePeriodicStuff();
	

	//This is the function that makes this class abstract, mesh is never read or generated within this class...
	virtual void getMesh() = 0;

};


#endif
