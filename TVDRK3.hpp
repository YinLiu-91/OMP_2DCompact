#ifndef _CTVDRK3H_
#define _CTVDRK3H_

#include "Utils.hpp"
#include "AbstractRK.hpp"

class TVDRK3:public AbstractRK{

    public:

	double a1, a2, a3;
	double b2, b3;
	double c1, c2, c3;

	int baseDirection;
	
	TVDRK3(AbstractCSolver *cs){

	    this->cs = cs;
	    baseDirection = cs->baseDirection;

	    Nx = cs->Nx;
	    Ny = cs->Ny;
	    N  = cs->N;

	    a1 = 1.0;
	    a2 = 3.0/4.0;
	    a3 = 1.0/3.0;

	    b2 = 1.0/4.0;
	    b3 = 2.0/3.0;

	    c1 = 1.0;
	    c2 = 1.0/4.0;
	    c3 = 2.0/3.0;
	   
	}


        void executeSolverLoop();
        void updateConservedData();

};


void TVDRK3::updateConservedData(){
 
//    if(cs->useTiming) cs->ft1 = MPI_Wtime();

    //Method outline (Gottlieb and Shu, 1998)
    // u(1)   =       u(n) + dt*L(u(n))
    // u(2)   = (3/4)*u(n) + (1/4)*u(1) + (1/4)*dt*L(u(1))
    // u(n+1) = (1/3)*u(n) + (2/3)*u(2) + (2/3)*dt*L(u(2))
    // TVD for CFL = 1

    if(cs->rkStep == 1){

        #pragma omp parallel for
        FOR_XY{
            //Calculate intermediate solution (u(1))
            cs->rhok[ip]  = cs->rho1[ip]  + c1*cs->rhok2[ip]; 
            cs->rhoUk[ip] = cs->rhoU1[ip] + c1*cs->rhoUk2[ip];
            cs->rhoVk[ip] = cs->rhoV1[ip] + c1*cs->rhoVk2[ip];
            cs->rhoEk[ip] = cs->rhoE1[ip] + c1*cs->rhoEk2[ip];
        }   

    }else if(cs->rkStep == 2){

        #pragma omp parallel for
        FOR_XY{
            //Calculate intermediate solution (u(2))
            cs->rhok[ip]  = a2*cs->rho1[ip]  + b2*cs->rhok[ip]  + c2*cs->rhok2[ip]; 
            cs->rhoUk[ip] = a2*cs->rhoU1[ip] + b2*cs->rhoUk[ip] + c2*cs->rhoUk2[ip];
            cs->rhoVk[ip] = a2*cs->rhoV1[ip] + b2*cs->rhoVk[ip] + c2*cs->rhoVk2[ip];
            cs->rhoEk[ip] = a2*cs->rhoE1[ip] + b2*cs->rhoEk[ip] + c2*cs->rhoEk2[ip];
        }

    }else if(cs->rkStep == 3){

        #pragma omp parallel for
        FOR_XY{
            //Get the final solution (u(n+1))
            cs->rho2[ip]  = a3*cs->rho1[ip]  + b3*cs->rhok[ip]  + c3*cs->rhok2[ip];
            cs->rhoU2[ip] = a3*cs->rhoU1[ip] + b3*cs->rhoUk[ip] + c3*cs->rhoUk2[ip];
            cs->rhoV2[ip] = a3*cs->rhoV1[ip] + b3*cs->rhoVk[ip] + c3*cs->rhoVk2[ip];
            cs->rhoE2[ip] = a3*cs->rhoE1[ip] + b3*cs->rhoEk[ip] + c3*cs->rhoEk2[ip];
        }

    }

}

void TVDRK3::executeSolverLoop(){

    cout << "HERE" << endl;

    while(cs->endFlag == false){

	cs->rkLast = false;
	
	cs->preStep();


	//Step 1
	cs->rkStep = 1;
	
	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();

	updateConservedData();

	cs->updateData();

	//Step 2
	cs->rkStep = 2;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 3
	cs->rkStep = 3;
	cs->rkLast = true;	

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	cs->postStep();	

    }
}
        
#endif
