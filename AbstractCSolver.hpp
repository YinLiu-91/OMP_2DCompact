#ifndef _CABSTRACTSOLVERH_
#define _CABSTRACTSOLVERH_

#include "Options.hpp"
#include "AbstractDerivatives.hpp"
#include "AbstractFilter.hpp"
#include "IdealGas.hpp"
#include "TimeStepping.hpp"
#include "AbstractSingleBlockMesh.hpp"
#include <vector>

//cyclic dependency here...
class AbstractSingleBlockMesh;

class AbstractCSolver{

    public:

	int rkStep;
	bool rkLast;
	bool endFlag;

	int Nx, Ny, N;

	int timeStep;
	double time;

	int baseDirection;

	//Stuff for timing functions
	bool useTiming;
	double ft1, ft2;
	 
	//Initial Conditions  
	double *rho0, *p0, *U0, *V0;

	//Common Solver Data
        double *rho1,  *rhok,  *rhok2,  *rho2;
        double *rhoU1, *rhoUk, *rhoUk2, *rhoU2;
        double *rhoV1, *rhoVk, *rhoVk2, *rhoV2;
        double *rhoE1, *rhoEk, *rhoEk2, *rhoE2;       

	double *p, *T, *U, *V;
 
	Options *opt;
	Domain *dom;
	BC *bc;
	TimeStepping *ts;
	IdealGas *ig;
        AbstractDerivatives *derivX, *derivY;
        AbstractFilter *filtX, *filtY;


	//Temporary data storage...
	double *tempX1,  *tempX2,  *tempX3,  *tempX4, *tempX5;
	double *tempX6,  *tempX7,  *tempX8,  *tempX9, *tempX10;

	double *tempY1,  *tempY2,  *tempY3,  *tempY4,  *tempY5;
	double *tempY6,  *tempY7,  *tempY8,  *tempY9,  *tempY10;
	double *tempY11, *tempY12, *tempY13, *tempY14, *tempY15;
	double *tempY16, *tempY17, *tempY18, *tempY19, *tempY20;
	double *tempY21, *tempY22, *tempY23, *tempY24;

	//Common Msh Stuff if needed
	AbstractSingleBlockMesh *msh;
	double *J, *J11, *J12, *J21, *J22;

	//double pointer list that we can attach things to
	vector<double*> varData;

	//Each Solver Class needs to have these functions to overwrite the pure virtual ones
	virtual void initializeSolverData() = 0;
	virtual void setInitialConditions() = 0;
	virtual void preStep() = 0;
	virtual void preSubStep() = 0;
	virtual void solveEqnSet() = 0;
	virtual void postSubStep() = 0;
	virtual void updateData() = 0;
	virtual void postStep() = 0;
	virtual void reportAll() = 0;

	//Hook functions
	virtual void initialHook(){};
	virtual void subStepTemporalHook(){};
	virtual void fullStepTemporalHook(){};
	virtual void preStepBoundaryHook(){};
	virtual void postStepBoundaryHook(){};

	//RHS source terms
	virtual double contRHSSource(int ip) = 0;
	virtual double xmomRHSSource(int ip) = 0;
	virtual double ymomRHSSource(int ip) = 0;
	virtual double engyRHSSource(int ip) = 0;
};

#endif
