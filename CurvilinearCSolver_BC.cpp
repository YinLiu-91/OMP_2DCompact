#include "CurvilinearCSolver.hpp"

void CurvilinearCSolver::setInitialConditions(){


//    if(useTiming) ft1 = MPI_Wtime();

    cout << endl;
    cout << " > Setting initial conditions..." << endl; 


    //Initialize sponge boundary conditions if necessary
    spg = NULL;

    //TODO NEED TO ADD EXCEPTIONS FOR WHEN WRONG KIND OF SPONGE OR UNIMPLEMENTED SPONGE BC IS TRYING TO BE USED
    if( bc->bcX0 == Options::SPONGE || \
        bc->bcX1 == Options::SPONGE || \
        bc->bcY0 == Options::SPONGE || \
        bc->bcY1 == Options::SPONGE){
        spongeFlag = true;
        spg = new SpongeBC(msh, dom, ig, bc, opt);
    }else{
        spongeFlag = false;
    }

    //just do the simple stuff in a loop...
    FOR_XY{
	U[ip]	 = U0[ip];
	V[ip] 	 = V0[ip];
	p[ip]	 = p0[ip];

	rho1[ip]  = rho0[ip];
	rhoU1[ip] = rho1[ip]*U[ip];	
	rhoV1[ip] = rho1[ip]*V[ip];	

	Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
    }

    //Can now release initial condition data...
    delete[] rho0;
    delete[] U0;
    delete[] V0;
    delete[] p0;

    //Call the ideal gas relations for the slightly more involved stuff..
    FOR_XY{
	rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip]);
    	T[ip]     = ig->solveT(rho1[ip], p[ip]);
        mu[ip]    = ig->solveMu(T[ip]);
        sos[ip]   = ig->solveSOS(rho1[ip], p[ip]);
    }
    //This is where we'll do the boundary condition specific stuff...
    bool wallBCFlag = false;

    //--------------------------------
    //No-Slip Wall Boundary Conditions
    //--------------------------------

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
	    if(bc->bcX0 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GET_Xp;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);

	    }
        }END_FORX0
    }


    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
	    if(bc->bcX1 == Options::CONST_T_WALL){

		int index[10] = FILL_GET_Xm;
	 	double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);

	    }
        }END_FORX1
    }

    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Y0{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
	    if(bc->bcY0 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GET_Yp;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
	        T[ip] = derivY->calcNeumann(T_out);

	    }
        }END_FORY0
    }

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){
	wallBCFlag = true;
        FOR_Y1{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
	    if(bc->bcY1 == Options::ADIABATIC_WALL){

		int index[10] = FILL_GET_Ym;
		double T_out[10];
		getDataFromIndex(T, index, 10, T_out);
	        T[ip] = derivY->calcNeumann(T_out);

	    }
        }END_FORY1
    }

    //-------------------------------
    //Moving Wall Boundary Conditions
    //-------------------------------

    if(bc->bcX0 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X0WallV;

	    int index[10] = FILL_GET_Xp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

        }END_FORX0
    }

    if(bc->bcX1 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X1WallV;

	    int index[10] = FILL_GET_Xm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

        }END_FORX1
    }

    if(bc->bcY0 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y0{
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Y0WallU;
            rhoV1[ip] = 0.0;

	    int index[10] = FILL_GET_Yp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

        }END_FORY0

   }

    if(bc->bcY1 == Options::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Y1{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Y1WallU;
            rhoV1[ip] = 0.0;

	    int index[10] = FILL_GET_Ym;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

        }END_FORY1

   }

    if(wallBCFlag == true){
	//Need to update the pressure, sos, and rhoE fields at the boundaries with walls...
	FOR_XY{
	    p[ip]     = ig->solvep_idealgas(rho1[ip], T[ip]);
	    sos[ip]   = ig->solveSOS(rho1[ip],p[ip]);
	    rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip]);

	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip]; 
	}
    }

    if(spongeFlag == true){

        if(opt->spongeFromRestart){

            spg->readSpongeAvgFromRestart();   

    	    getRange(spg->spongeRhoAvg,   "spongeRhoAvg",  Nx, Ny); 
    	    getRange(spg->spongeRhoUAvg,  "spongeRhoUAvg", Nx, Ny); 
    	    getRange(spg->spongeRhoVAvg,  "spongeRhoVAvg", Nx, Ny); 
    	    getRange(spg->spongeRhoEAvg,  "spongeRhoEAvg", Nx, Ny); 

	}else{
	    FOR_XY{
	        spg->spongeRhoAvg[ip]  = rho1[ip];
	        spg->spongeRhoUAvg[ip] = rhoU1[ip];
	        spg->spongeRhoVAvg[ip] = rhoV1[ip];
	        spg->spongeRhoEAvg[ip] = rhoE1[ip];
 	    }
	}
    }



    std::cout << " > Finished initialization of flow field " << std::endl;

    getRange(rho1,  "rho0",  Nx, Ny); 
    getRange(rhoU1, "rhoU0", Nx, Ny); 
    getRange(rhoV1, "rhoV0", Nx, Ny); 
    getRange(rhoE1, "rhoE0", Nx, Ny); 
    getRange(T, "T0", Nx, Ny); 
    getRange(p, "p0", Nx, Ny); 

    //Run the initial hook
    initialHook();

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > setInitCond Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/

}



void CurvilinearCSolver::preStepBCHandling(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP, *rhoUP, *rhoVP, *rhoEP;
    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoEP = rhoE1;
    }else{
	rhoP  = rhok;
	rhoUP = rhoUk; 
	rhoVP = rhoVk;
	rhoEP = rhoEk;
    }


    //--------------------------------
    //No-slip wall boundary conditions
    //--------------------------------

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	FOR_X0{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    //rhoP[ip] = rhoP[GETMAJIND_YPEN_Xp1];
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    if(bc->bcX0 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GET_Xp;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);
	    }
	   p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	   rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	   Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	   Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip]; 

	}END_FORX0

    }
    

    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	FOR_X1{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    //rhoP[ip] = rhoP[GETMAJIND_YPEN_Xm1];
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    if(bc->bcX1 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GET_Xm;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivX->calcNeumann(T_out);
	    }
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];

	}END_FORX1
    }   


    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){

	FOR_Y0{
	    if(bc->bcY0 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GET_Yp;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivY->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    //rhoP[ip] = rhoP[GETMAJIND_YPEN_Yp1];
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]);
   	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip];
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip]; 
	}END_FORY0
    }
    

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){

	FOR_Y1{

	    if(bc->bcY1 == Options::ADIABATIC_WALL){
	        int index[10] = FILL_GET_Ym;
	        double T_out[10];
	        getDataFromIndex(T, index, 10, T_out);
                T[ip] = derivY->calcNeumann(T_out);
	    }
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    //rhoP[ip] = rhoP[GETMAJIND_YPEN_Ym1];
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip];
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];

	}END_FORY1
    }

    //-------------------------------
    //Moving wall boundary conditions
    //-------------------------------

    if(bc->bcX0 == Options::MOVING_ADIABATIC_WALL){
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X0WallV;

	    int index[10] = FILL_GET_Xp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
        }END_FORX0
    }

    if(bc->bcX1 == Options::MOVING_ADIABATIC_WALL){
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X1WallV;

	    int index[10] = FILL_GET_Xm;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivX->calcNeumann(T_out);

	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
        }END_FORX1
    }

    if(bc->bcY0 == Options::MOVING_ADIABATIC_WALL){

        FOR_Y0{

	    int index[10] = FILL_GET_Yp;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Y0WallU;
            rhoVP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip];
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip]; 
        }END_FORY0
    }

    if(bc->bcY1 == Options::MOVING_ADIABATIC_WALL){

        FOR_Y1{

	    int index[10] = FILL_GET_Ym;
	    double T_out[10];
	    getDataFromIndex(T, index, 10, T_out);
            T[ip] = derivY->calcNeumann(T_out);

            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Y1WallU;
            rhoVP[ip] = 0.0;
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip]); //Not 100% about this?
  	    Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip]; 
	    Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
        }END_FORY1
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > preBCHandl Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/


}



void CurvilinearCSolver::postStepBCHandling(){

//    if(useTiming) ft1 = MPI_Wtime();

    ////////////////////////////////
    //ADIABATIC AND MOVING WALL BC// 
    ////////////////////////////////

    if(bc->bcX0 == Options::ADIABATIC_WALL || bc->bcX0 == Options::MOVING_ADIABATIC_WALL || bc->bcX0 == Options::CONST_T_WALL){
	FOR_X0{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_1[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX0
    }

    if(bc->bcX1 == Options::ADIABATIC_WALL || bc->bcX1 == Options::MOVING_ADIABATIC_WALL || bc->bcX1 == Options::CONST_T_WALL){
	FOR_X1{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_1[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX1
    }   

    if(bc->bcY0 == Options::ADIABATIC_WALL || bc->bcY0 == Options::MOVING_ADIABATIC_WALL || bc->bcY0 == Options::CONST_T_WALL){
	FOR_Y0{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_2[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY0
    }

    if(bc->bcY1 == Options::ADIABATIC_WALL || bc->bcY1 == Options::MOVING_ADIABATIC_WALL || bc->bcY1 == Options::CONST_T_WALL){
	FOR_Y1{
	    rhok2[ip]  = -ts->dt*J[ip]*cont_2[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY1
    }

    /////////////
    //SPONGE BC//
    /////////////

    if(bc->bcX0 == Options::SPONGE){
	FOR_X0{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX0
    }

    if(bc->bcX1 == Options::SPONGE){
	FOR_X1{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX1
    }   

    if(bc->bcY0 == Options::SPONGE){
	FOR_Y0{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY0
    }

    if(bc->bcY1 == Options::SPONGE){
	FOR_Y1{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY1
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > postBCHand Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/
}


void CurvilinearCSolver::updateSponge(){

//    if(useTiming) ft1 = MPI_Wtime();

    if(spongeFlag){
	double eps = 1.0/(spg->avgT/ts->dt + 1.0);
	FOR_XY{
	    spg->spongeRhoAvg[ip]  += eps*(rho1[ip]  - spg->spongeRhoAvg[ip]);	
	    spg->spongeRhoUAvg[ip] += eps*(rhoU1[ip] - spg->spongeRhoUAvg[ip]);	
	    spg->spongeRhoVAvg[ip] += eps*(rhoV1[ip] - spg->spongeRhoVAvg[ip]);	
	    spg->spongeRhoEAvg[ip] += eps*(rhoE1[ip] - spg->spongeRhoEAvg[ip]);	
	    spg->spongeRhoEAvg[ip] = spg->epsP*spg->spongeRhoEAvg[ip] + (1.0 -  spg->epsP)*(spg->spongeP/(ig->gamma-1.0) \
					 + 0.5*(spg->spongeRhoUAvg[ip]*spg->spongeRhoUAvg[ip] + spg->spongeRhoVAvg[ip]*spg->spongeRhoVAvg[ip] \
					 )/spg->spongeRhoAvg[ip]);
	}
	
        if(bc->bcX0 == Options::SPONGE){
	    FOR_X0{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX0
        }

        if(bc->bcX1 == Options::SPONGE){
	    FOR_X1{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX1
        }   

        if(bc->bcY0 == Options::SPONGE){
	    FOR_Y0{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY0
        }

        if(bc->bcY1 == Options::SPONGE){
	    FOR_Y1{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY1
        }
    }
/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > updateSpge Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/
};


