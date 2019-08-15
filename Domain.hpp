#ifndef _DOMAINH_
#define _DOMAINH_

#include <iostream>
#include "Options.hpp"
#include "BC.hpp"
#include "Macros.hpp"

class Domain{

    public:

	//Global variables are denoted with a g
        int Nx, Ny, N;
        double Lx, Ly; 

	//These are going to be based on the global indices
        double *x, *y;

	double dx, dy;

    Domain(Options *opt, double Lx, double Ly){
	
	Nx = opt->Nx;
	Ny = opt->Ny;
	N = Nx*Ny;

	Lx = Lx;
	Ly = Ly;

	x = new double[Nx];	
	y = new double[Ny];	

	if(opt->bcXType == Options::PERIODIC_SOLVE){
	    for(int ip = 0; ip < Nx; ip++){
		x[ip] = (double)ip*Lx/(double)Nx;
	    }
	}else{
            for(int ip = 0; ip < Nx; ip++){
	        x[ip] = (((double)ip)/((double)Nx - 1.0))*Lx;
	    }	
	}

	if(opt->bcYType == Options::PERIODIC_SOLVE){
	    for(int jp = 0; jp < Ny; jp++){
	        y[jp] = (double)jp*Ly/(double)Ny;
	    }
	}else{
            for(int jp = 0; jp < Ny; jp++){
	        y[jp] = (((double)jp)/((double)Ny - 1.0))*Ly;
	    }	
	}

	dx = x[1]-x[0];
	dy = y[1]-y[0];

	std::cout << " > Global Domain initialization..." << std::endl;
	std::cout << " > Domain: " << Lx << "x" << Ly << std::endl;
	std::cout << " > Mesh: " << Nx << "x" << Ny << std::endl;
	std::cout << " > Total Points: " << Nx*Ny << std::endl; 
    }
  
};



#endif
