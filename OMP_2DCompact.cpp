#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

#include "Macros.hpp"
#include "Options.hpp"
#include "TimeStepping.hpp"
#include "Domain.hpp"
#include "BC.hpp" 

#include "AbstractCSolver.hpp"
#include "CurvilinearCSolver.hpp"

#include "AbstractSingleBlockMesh.hpp"
#include "AlgebraicSingleBlockMesh.hpp"

#include "AbstractRK.hpp"
#include "TVDRK3.hpp"
//#include "RK4.hpp"
//#include "KenRK4.hpp"
#include "LSLDDRK4.hpp"


int main(int argc, char *argv[]){

    int numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }


    cout << endl;
    cout << " --------------------" << endl;
    cout << "  Starting up Solver " << endl;
    cout << " --------------------" << endl;
    cout << endl;
    cout << " > Running on " << numThreads << " threads!" << endl;
    cout << endl;

    //Have the root rank parse the input file and broadcast it out...
    Options *opt = new Options;

    
    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping *ts = new TimeStepping(opt);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////

    bool periodicBC[2];
    BC *bc = new BC(opt, periodicBC);


    int TS, Nx, Ny;
    double time;
    if(opt->fromRestart || opt->onlyGridFromRestart){

	//If we need to read from a file, pull the dimensions from
	//the leading three doubles from the file...

	double cN[4];
	FILE *ptr;
	ptr = fopen(opt->filename.c_str(), "rb");
	if(ptr == NULL){
	    cout << "ERROR: Couldn't open file " << opt->filename << endl;
	    abort();
	}else{
	    fread(cN, sizeof(double), 4, ptr);
	}
	fclose(ptr);

	TS = (int)cN[0];
	time = cN[1];
	Nx = (int)cN[2];
	Ny = (int)cN[3];

	if(opt->fromRestart){
	    opt->timeStep = TS;
	    opt->time = time;
	}else{
	    opt->timeStep = 0;
	    opt->time = 0.0;
	}
	opt->Nx = Nx;
	opt->Ny = Ny;

    }else{
	TS = 0;
	time = 0.0;
	opt->timeStep = TS;
	opt->time = time;
	Nx = opt->Nx;
	Ny = opt->Ny;
    }

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    //For curvilinear coordinates these should all correspond to the max xi, eta, and zeta values
    double Lx = 1.0, Ly = 1.0;
    Domain *d = new Domain(opt, Lx, Ly);

    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    AbstractCSolver *cs;
    cs = new CurvilinearCSolver(d, bc, ts, opt);

    //Attach the mesh object to the solver...
    cs->msh = new AlgebraicSingleBlockMesh(cs, d);

    //ADT testing
    double p[2] = {2.0, 2.0};
    int icv = cs->msh->findCVForPoint(p);
    cout << icv << ", " << cs->msh->x[icv] << "," << cs->msh->y[icv] << endl;

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk;

    if(opt->rkType == Options::TVDRK3){
        rk = new TVDRK3(cs);
    }else if(opt->rkType == Options::RK4){
//	rk = new RK4(cs);
    }else if(opt->rkType == Options::KENRK4){
//        rk = new KenRK4(cs);
    }else if(opt->rkType == Options::LSLDDRK4){
        rk = new LSLDDRK4(cs);
    }else{
	cout << "SHOULD NEVER GET HERE!" << endl;
    }

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////

    bool fromRestart = opt->fromRestart;
    if(!fromRestart){

	#pragma omp parallel for collapse(2)
        FOR_Y{
            FOR_X{

                int ip = GET2DINDEX_XY;

		double r2, pt, rt;
		
//		r2 = (cs->msh->x[ip] - M_PI)*(cs->msh->x[ip] - M_PI) + (cs->msh->y[ip] - M_PI)*(cs->msh->y[ip] - M_PI); 
//		r2 = (cs->msh->x[ip] - M_PI)*(cs->msh->x[ip] - M_PI); 

		if(cs->msh->x[ip] > 0.5*M_PI && cs->msh->x[ip] < 1.5*M_PI && cs->msh->y[ip] > 0.5*M_PI && cs->msh->y[ip] < 1.5*M_PI){
		    pt = 1.0;
		    rt = 1.0;
		}else{
		    pt = 0.125;
		    rt = 0.1;
		}

                cs->rho0[ip] = rt;
                cs->p0[ip]   = pt/cs->ig->gamma;
                cs->U0[ip]   = 0.0;
                cs->V0[ip]   = 0.0;
            }
        }
    }else{

	string filename = opt->filename;

	/* Need to set up just a serial read of the data...
	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	//Need to displace by the first 3 doubles and then three double fields
	disp = 5*sizeof(double)+ sizeof(double)*3.0*opt->Nx*opt->Ny*opt->Nz;

	IF_RANK0{
	    cout << " " << endl;
	}

	double *rho_in,
	       *rhoU_in,
	       *rhoV_in,
	       *rhoW_in,
	       *rhoE_in;

	cs->c2d->allocY(rho_in);
	cs->c2d->allocY(rhoU_in);
	cs->c2d->allocY(rhoV_in);
	cs->c2d->allocY(rhoW_in);
	cs->c2d->allocY(rhoE_in);

	cs->c2d->readVar(fh, disp, 1, rho_in);
	cs->c2d->readVar(fh, disp, 1, rhoU_in);
	cs->c2d->readVar(fh, disp, 1, rhoV_in);
	cs->c2d->readVar(fh, disp, 1, rhoW_in);
	cs->c2d->readVar(fh, disp, 1, rhoE_in);

	MPI_File_close(&fh);

        FOR_Z_YPEN{
            FOR_Y_YPEN{
                FOR_X_YPEN{

                    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

	            cs->rho0[ip] = rho_in[jp];
	            cs->U0[ip] = cs->ig->solveU(rho_in[jp], rhoU_in[jp]);
	            cs->V0[ip] = cs->ig->solveU(rho_in[jp], rhoV_in[jp]);
	            cs->W0[ip] = cs->ig->solveU(rho_in[jp], rhoW_in[jp]);
	            cs->p0[ip] = cs->ig->solvep(rho_in[jp], rhoE_in[jp], cs->U0[ip], cs->V0[ip], cs->W0[ip]);

		}
	    }
	}

	cs->c2d->deallocXYZ(rho_in);
	cs->c2d->deallocXYZ(rhoU_in);
	cs->c2d->deallocXYZ(rhoV_in);
	cs->c2d->deallocXYZ(rhoW_in);
	cs->c2d->deallocXYZ(rhoE_in);
*/
    }

    cs->setInitialConditions();



    ///////////////////////////
    //Add images to be output//
    ///////////////////////////


    //This is probably bad programming, but we'll downcast the abstract solver class pointer to the
    //solver pointer so we can access the add image function and the solver member we want to print out

    CurvilinearCSolver *cs_downcast = static_cast<CurvilinearCSolver*>(cs);

    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->p, "P", PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->rho2, "RHO", PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->rho2, "RHOBOW", PngWriter::RAINBOW));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->U, "U", PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->T, "T", PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->mu_star,   "MUS",   PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->beta_star, "BETAS", PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->k_star,    "KS",    PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->fsw,    "FSW",    PngWriter::GREYSCALE));
    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->vort,    "VORT",    PngWriter::GREYSCALE));
//    cs_downcast->addImageOutput(new PngWriter(1, 1024, 1024, cs_downcast->lad->dFbeta4dx04,    "dFbeta0",    PngWriter::GREYSCALE));
//    cs_downcast->addImageOutput(new PngWriter(1, 1024, 1024, cs_downcast->lad->dFbeta4dx14,    "dFbeta1",    PngWriter::GREYSCALE));
//    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->dFmu4dx04,    "dFmu0",    PngWriter::GREYSCALE));
//    cs_downcast->addImageOutput(new PngWriter(25, 1024, 1024, cs_downcast->lad->dFmu4dx14,    "dFmu1",    PngWriter::GREYSCALE));
//    cs_downcast->addImageOutput(new PngWriter(1, 1024, 1024, cs_downcast->preMom2_2, "premom22_O0",    PngWriter::GREYSCALE));


    ////////////////////////////////////////
    //Execute the solver timestepping loop//
    ////////////////////////////////////////

    rk->executeSolverLoop();  

    return 0;
}

void CurvilinearCSolver::initialHook(){};
void CurvilinearCSolver::fullStepTemporalHook(){};
void CurvilinearCSolver::subStepTemporalHook(){};
void CurvilinearCSolver::preStepBoundaryHook(){};
void CurvilinearCSolver::postStepBoundaryHook(){};
double CurvilinearCSolver::contRHSSource(int ip){

	return 0.0;
};

double CurvilinearCSolver::xmomRHSSource(int ip){return 0.0;};
double CurvilinearCSolver::ymomRHSSource(int ip){

/*
	double x = msh->x[ip];
	double y = msh->y[ip];

	double r2 = (x - 3.0*M_PI)*(x - 3.0*M_PI) + (y - 2.0*M_PI)*(y - 2.0*M_PI);
	double pt = 250.*exp(-r2/0.005);

	double p_source;
        p_source = sin(4.0*time)*pt; 	
*/
	double p_source = 0.0;

	return p_source;
};

double CurvilinearCSolver::engyRHSSource(int ip){

	return 0.0;
};







