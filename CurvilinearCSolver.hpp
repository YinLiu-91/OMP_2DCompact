#ifndef _CURVILINEARCSOLVERH_
#define _CURVILINEARCSOLVERH_

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "Macros.hpp"
#include "Options.hpp"
#include "Utils.hpp"
#include "SpongeBC.hpp"
#include "AbstractCSolver.hpp"
#include "PngWriter.hpp"
#include "CurvilinearInterpolator.hpp"
#include "Pade6.hpp"
//#include "Penta10.hpp"
//#include "CD2.hpp"
//#include "Compact8Filter.hpp"
#include "Compact10Filter.hpp"
//#include "Stats.hpp"
//#include "AbstractSGS.hpp"
//#include "VremanSGS.hpp"
//#include "DSMSGS.hpp"

class CurvilinearCSolver: public AbstractCSolver{

    public:

	double alphaF;
	double mu_ref;

	//Make local copies for macros...

	double t1, t2;

	
	
	//Track the current time and timestep
	int filterTimeStep;

	//Kill solver condition
        bool done;

	//non-conserved data
	double  *U, *Ucurv,
		*V, *Vcurv,
	        *T,  
		*p,
	       *mu,
	      *sos;

	//derivatives of data
	double *dU1, *dU2;
	double *dV1, *dV2;
	double *dW1, *dW2;
	double *dT1, *dT2;


	double *Tau11,  *Tau12;
	double *Tau21,  *Tau22;

	double *cont_1, *cont_2;
	double *mom1_1, *mom1_2;
	double *mom2_1, *mom2_2;
	double *engy_1, *engy_2;

	vector<double*> tempVec;

	//Stuff for different styles of interprocessor computation...
	enum CompStyle {VANILLA};
        CompStyle compStyle;

	bool spongeFlag;
	SpongeBC *spg; 

	//Moving Wall BC Velocities
	double X0WallV, X1WallV;
	double Y0WallU, Y1WallU;

	//For drawing images
	list<PngWriter*> imageList;

	//Statstics object
	//Stats *stats;
	//bool statsFlag;

	//LES Object
	//AbstractSGS *les;
	//bool LESFlag;
	//double Pr_t;

	//Alias'd derivative objects
	AbstractDerivatives *derivXi1, *derivXi2;
	AbstractFilter *filtXi1, *filtXi2;

	//Constructor to use for this class...
	CurvilinearCSolver(Domain *dom, BC *bc, TimeStepping *ts, Options *opt){

	    //Take in input information and initialize data structures...
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->opt = opt;
	    this->alphaF = opt->alphaF;
	    this->mu_ref = opt->mu_ref;
	    this->useTiming = opt->useTiming;

	    baseDirection = 0;

	    ig = new IdealGas(dom, mu_ref);

	    //give ourselves the local copies of the domain sizes
	    Nx = dom->Nx;
	    Ny = dom->Ny;
	    N  = dom->N;

	    //initialize time and timestep
   	    time = opt->time;
	    timeStep = opt->timeStep;
	    filterTimeStep = 0;
	    endFlag = false;
	    done = false;
	    rkLast = false;

	    //Computational style, should be in inputs
	    compStyle = VANILLA;
	
	    //Allocate our arrays for the solver data
	    initializeSolverData();		    	    

	    //Initialize our derivative calculations for each direction...
            derivX = new Pade6(dom, bc->bcXType, AbstractDerivatives::DIRX);
  	    derivY = new Pade6(dom, bc->bcYType, AbstractDerivatives::DIRY);

	    //if(opt->xFDType == Options::CD2){
	    //    derivX = new CD2(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    //}else if(opt->xFDType == Options::PADE6){
	    //    derivX = new Pade6(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    //}else if(opt->xFDType == Options::PENTA10){
	    //    derivX = new Penta10(dom, bc->bcXType, AbstractDerivatives::DIRX);
	    //}else{
	    //	cout << "Should never get here? unknown x-derivative" << endl;
 	    //	MPI_Abort(MPI_COMM_WORLD, -10);
	    //}
/*
	    if(opt->yFDType == Options::CD2){
	        derivY = new CD2(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else if(opt->yFDType == Options::PADE6){
	        derivY = new Pade6(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else if(opt->yFDType == Options::PENTA10){
	        derivY = new Penta10(dom, bc->bcYType, AbstractDerivatives::DIRY);
	    }else{
		cout << "Should never get here? unknown y-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);
	    }
*/

	    derivXi1 = derivX;
	    derivXi2 = derivY;

	    //Initialize the filters we're going to use for each direction
	    filtX  = new Compact10Filter(alphaF, dom, bc, bc->bcXType, AbstractDerivatives::DIRX);
	    filtY  = new Compact10Filter(alphaF, dom, bc, bc->bcYType, AbstractDerivatives::DIRY);

/*
	    if(opt->filterType == Options::COMPACT8){
	        filtX  = new Compact8Filter(alphaF, dom, bc, bc->bcXType, AbstractDerivatives::DIRX);
	        filtY  = new Compact8Filter(alphaF, dom, bc, bc->bcYType, AbstractDerivatives::DIRY);
	        filtZ  = new Compact8Filter(alphaF, dom, bc, bc->bcZType, AbstractDerivatives::DIRZ);
	    }else if(opt->filterType == Options::COMPACT10){
	        filtX  = new Compact10Filter(alphaF, dom, bc, bc->bcXType, AbstractDerivatives::DIRX);
	        filtY  = new Compact10Filter(alphaF, dom, bc, bc->bcYType, AbstractDerivatives::DIRY);
	        filtZ  = new Compact10Filter(alphaF, dom, bc, bc->bcZType, AbstractDerivatives::DIRZ);

	    }else{
		cout << "Should never get here? unknown z-derivative" << endl;
		MPI_Abort(MPI_COMM_WORLD, -10);

	    }
*/
	    filtXi1 = filtX;
	    filtXi2 = filtY;

 	    X0WallV = 0.0; X1WallV = 0.0; 
	    Y0WallU = 0.0; Y1WallU = 0.0; 

/*
	    //Initialize the statistics object
	    statsFlag = opt->velocityStats || opt->thermoStats;
	    if(statsFlag){
	        stats = new Stats(this, opt);
	    }else{
		stats = NULL;
	    }

	    //Initialize the LES stuff if we need to...
	    Pr_t = 0.9;
	    if(opt->lesModel == Options::NOMODEL){
		LESFlag = false;
		les = NULL;
	    }else if(opt->lesModel == Options::VREMAN){
		LESFlag = true;
		les = new VremanSGS(this);
	    }else if(opt->lesModel == Options::DSM){
		LESFlag = true;
		les = new DSMSGS(this, opt->testFilterType, opt->useTaukk, opt->dumpCoeffRange);
	    } 

	    t1 = MPI_Wtime();
*/
	}

	//Pre solver utility functions
	void addImageOutput(PngWriter *pw);
	void generateImagePlane(PngWriter *pw);

	//Pre Step Functions...
	void calcDtFromCFL();

	//Pre Sub Step Functions...
	void preStepBCHandling();
	void preStepDerivatives();

	//Solve Eqn Set Functions...
	void solveContinuity();
	void solveXMomentum();
	void solveYMomentum();
	void solveEnergy();

	//Post Sub Step Functions...
	void postStepBCHandling();
	void filterConservedData();
	void updateNonConservedData();

	//Post Step Functions
	void updateSponge();
	void writeImages();
	void writePlaneImageForVariable(PngWriter *pw);
	void checkSolution();
	void dumpSolution();
	void checkEnd();
	void reportAll();

	//Hook functions
	void initialHook();
	void fullStepTemporalHook();
	void subStepTemporalHook();
	void preStepBoundaryHook();
	void postStepBoundaryHook();

	double contRHSSource(int ip);
	double xmomRHSSource(int ip);
	double ymomRHSSource(int ip);
	double engyRHSSource(int ip);

	//Inline, Or general solver functions
	inline double calcSpongeSource(double phi, double phiSpongeAvg, double sigma){
        	return sigma*(phiSpongeAvg - phi);
	};

	//Derivative calculation functions
	void computeGradient(vector<double*> vecIn, vector<double*>vecOut);
	void computeGradDotComponents(vector<double*> vecIn, vector<double*>vecOut);

	bool checkForAndDeleteKillFile(string killFileName);

	/////////////////////////////////////
	//Our Generalized Solver Functions //
	/////////////////////////////////////

	void setInitialConditions();
	void initializeSolverData();

	void preStep(){

   	    if(timeStep == opt->timeStep){
//        	dumpSolution();
        	writeImages();
    	    }
    	    calcDtFromCFL();

	}

	void preSubStep(){

    	    preStepBCHandling();
	    preStepBoundaryHook();

    	    preStepDerivatives();

	    //reportAll();
	}

	void solveEqnSet(){

    	    solveContinuity();
            solveXMomentum();
    	    solveYMomentum();
    	    solveEnergy();
	    
	}

	void postSubStep(){

    	    postStepBCHandling();
	    postStepBoundaryHook();
	    subStepTemporalHook();
	}

	void updateData(){
 
     	    if(rkLast || opt->subStepFiltering){
	       	filterConservedData();
   	    }
 
	    updateNonConservedData();
	}

	void postStep(){


	    fullStepTemporalHook();

    	    updateSponge();
    	    checkSolution();

	    //if(timeStep%opt->stats_interval == 0 && statsFlag){
		//stats->updateStatsFields();
	    //} 


    	    //if(timeStep%ts->dumpStep == 0){
       // 	dumpSolution();
//	    }

       	    writeImages();

    	    checkEnd();

	}



};



#endif


