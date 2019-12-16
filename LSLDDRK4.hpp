#ifndef _CLSLLDRK4H_
#define _CLSLLDRK4H_

#include "Utils.hpp"
#include "AbstractRK.hpp"

//LSLDDRK4 = Low Storage, Low Dissipation, Dispersion RK4
//Low Dissipation, Dispersion 2-Step 11-Stage RK4
//This is the two-register low storage scheme of Stanescu & Habashi


class LSLDDRK4:public AbstractRK{

    public:

	double a1_1, a2_1, a3_1, a4_1, a5_1;
	double b1_1, b2_1, b3_1, b4_1, b5_1;
	double c1_1, c2_1, c3_1, c4_1, c5_1;

	double a1_2, a2_2, a3_2, a4_2, a5_2, a6_2;
	double b1_2, b2_2, b3_2, b4_2, b5_2, b6_2;
	double c1_2, c2_2, c3_2, c4_2, c5_2, c6_2;

	int rkDualStepStep;

	LSLDDRK4(AbstractCSolver *cs){

	    this->cs = cs;

	    Nx = cs->Nx;
	    Ny = cs->Ny;
	    N  = cs->N;
	

	    //Five Stage RK4 parameters
	    a1_1 = 0.0;
	    a2_1 = -0.60512264332862261228;
	    a3_1 = -2.04375640234761394333;
	    a4_1 = -0.74069990637544192841;
	    a5_1 = -4.42317651302968168941;

	    b1_1 = 0.26874543888713438496;
	    b2_1 = 0.8014706973220802933;
	    b3_1 = 0.50515704269422722538;
	    b4_1 = 0.56235680379000296407;
	    b5_1 = 0.05900655127758823335;

	    c1_1 = 0.0;
	    c2_1 = 0.26874543888713438496;
	    c3_1 = 0.58522806929524303469;
	    c4_1 = 0.68270664478424678821;
	    c5_1 = 1.1646854837729261436;

	    //Six-Stage RK4 parameters
	    a1_2 = 0.0;
	    a2_2 = -0.44127377153877382565;
	    a3_2 = -1.073982008079781868;
	    a4_2 = -1.7063570791256758809;
	    a5_2 = -2.7979293162682443056;
	    a6_2 = -4.0913537120919160454;

	    b1_2 = 0.11584888181285561688;
	    b2_2 = 0.37287699051652864918;
	    b3_2 = 0.73795368921435295698;
	    b4_2 = 0.57981109366311039538;
	    b5_2 = 1.031284991300145194;
	    b6_2 = 0.15;

	    c1_2 = 0.0;
	    c2_2 = 0.11584888181285561688;
	    c3_2 = 0.32418503640412806853;
	    c4_2 = 0.61932082035177792368;
	    c5_2 = 0.80344726663359079059;
	    c6_2 = 0.91841664452065965078;


	}


        void executeSolverLoop();
        void updateConservedData();
	void updateConservedDataStep1();
	void updateConservedDataStep2();

};

void LSLDDRK4::updateConservedDataStep1(){



    if(cs->rkStep == 1){

	    #pragma omp parallel for
            FOR_XY{
		cs->rhok[ip]  = cs->rho1[ip]  + b1_1*cs->rhok2[ip]; 
                cs->rhoUk[ip] = cs->rhoU1[ip] + b1_1*cs->rhoUk2[ip];
                cs->rhoVk[ip] = cs->rhoV1[ip] + b1_1*cs->rhoVk2[ip];
                cs->rhoEk[ip] = cs->rhoE1[ip] + b1_1*cs->rhoEk2[ip];
	
		cs->rho1[ip]   = cs->rhok2[ip]; 	
		cs->rhoU1[ip]  = cs->rhoUk2[ip]; 	
		cs->rhoV1[ip]  = cs->rhoVk2[ip]; 	
		cs->rhoE1[ip]  = cs->rhoEk2[ip]; 	
            }

    }else if(cs->rkStep == 2){
	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a2_1*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a2_1*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a2_1*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a2_1*cs->rhoE1[ip];

		cs->rhok[ip]  += b2_1*cs->rho1[ip]; 
                cs->rhoUk[ip] += b2_1*cs->rhoU1[ip];
                cs->rhoVk[ip] += b2_1*cs->rhoV1[ip];
                cs->rhoEk[ip] += b2_1*cs->rhoE1[ip];
            }

    }else if(cs->rkStep == 3){

	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a3_1*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a3_1*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a3_1*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a3_1*cs->rhoE1[ip];

		cs->rhok[ip]  += b3_1*cs->rho1[ip]; 
                cs->rhoUk[ip] += b3_1*cs->rhoU1[ip];
                cs->rhoVk[ip] += b3_1*cs->rhoV1[ip];
                cs->rhoEk[ip] += b3_1*cs->rhoE1[ip];
	    }

    }else if(cs->rkStep == 4){
	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a4_1*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a4_1*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a4_1*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a4_1*cs->rhoE1[ip];

		cs->rhok[ip]  += b4_1*cs->rho1[ip]; 
                cs->rhoUk[ip] += b4_1*cs->rhoU1[ip];
                cs->rhoVk[ip] += b4_1*cs->rhoV1[ip];
                cs->rhoEk[ip] += b4_1*cs->rhoE1[ip];

	}
    }else if(cs->rkStep == 5){
	    #pragma omp parallel for
	    FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a5_1*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a5_1*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a5_1*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a5_1*cs->rhoE1[ip];

		cs->rho2[ip]  = cs->rhok[ip]  + b5_1*cs->rho1[ip]; 
                cs->rhoU2[ip] = cs->rhoUk[ip] + b5_1*cs->rhoU1[ip];
                cs->rhoV2[ip] = cs->rhoVk[ip] + b5_1*cs->rhoV1[ip];
                cs->rhoE2[ip] = cs->rhoEk[ip] + b5_1*cs->rhoE1[ip];

	    }
    }


};

void LSLDDRK4::updateConservedDataStep2(){


    if(cs->rkStep == 1){

	    #pragma omp parallel for
            FOR_XY{
		cs->rhok[ip]  = cs->rho1[ip]  + b1_2*cs->rhok2[ip]; 
                cs->rhoUk[ip] = cs->rhoU1[ip] + b1_2*cs->rhoUk2[ip];
                cs->rhoVk[ip] = cs->rhoV1[ip] + b1_2*cs->rhoVk2[ip];
                cs->rhoEk[ip] = cs->rhoE1[ip] + b1_2*cs->rhoEk2[ip];
	
		cs->rho1[ip]   = cs->rhok2[ip]; 	
		cs->rhoU1[ip]  = cs->rhoUk2[ip]; 	
		cs->rhoV1[ip]  = cs->rhoVk2[ip]; 	
		cs->rhoE1[ip]  = cs->rhoEk2[ip]; 	
            }

    }else if(cs->rkStep == 2){
	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a2_2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a2_2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a2_2*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a2_2*cs->rhoE1[ip];

		cs->rhok[ip]  += b2_2*cs->rho1[ip]; 
                cs->rhoUk[ip] += b2_2*cs->rhoU1[ip];
                cs->rhoVk[ip] += b2_2*cs->rhoV1[ip];
                cs->rhoEk[ip] += b2_2*cs->rhoE1[ip];
            }

    }else if(cs->rkStep == 3){

	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a3_2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a3_2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a3_2*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a3_2*cs->rhoE1[ip];

		cs->rhok[ip]  += b3_2*cs->rho1[ip]; 
                cs->rhoUk[ip] += b3_2*cs->rhoU1[ip];
                cs->rhoVk[ip] += b3_2*cs->rhoV1[ip];
                cs->rhoEk[ip] += b3_2*cs->rhoE1[ip];
	    }

    }else if(cs->rkStep == 4){
	    #pragma omp parallel for
            FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a4_2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a4_2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a4_2*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a4_2*cs->rhoE1[ip];

		cs->rhok[ip]  += b4_2*cs->rho1[ip]; 
                cs->rhoUk[ip] += b4_2*cs->rhoU1[ip];
                cs->rhoVk[ip] += b4_2*cs->rhoV1[ip];
                cs->rhoEk[ip] += b4_2*cs->rhoE1[ip];

	}
    }else if(cs->rkStep == 5){
	    #pragma omp parallel for
	    FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a5_2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a5_2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a5_2*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a5_2*cs->rhoE1[ip];

		cs->rhok[ip]  += b5_2*cs->rho1[ip]; 
                cs->rhoUk[ip] += b5_2*cs->rhoU1[ip];
                cs->rhoVk[ip] += b5_2*cs->rhoV1[ip];
                cs->rhoEk[ip] += b5_2*cs->rhoE1[ip];

	    }
    }else if(cs->rkStep == 6){
	    #pragma omp parallel for
	    FOR_XY{
                cs->rho1[ip]  = cs->rhok2[ip]  + a6_2*cs->rho1[ip]; 
                cs->rhoU1[ip] = cs->rhoUk2[ip] + a6_2*cs->rhoU1[ip];
                cs->rhoV1[ip] = cs->rhoVk2[ip] + a6_2*cs->rhoV1[ip];
                cs->rhoE1[ip] = cs->rhoEk2[ip] + a6_2*cs->rhoE1[ip];

		cs->rho2[ip]  = cs->rhok[ip]  + b6_2*cs->rho1[ip]; 
                cs->rhoU2[ip] = cs->rhoUk[ip] + b6_2*cs->rhoU1[ip];
                cs->rhoV2[ip] = cs->rhoVk[ip] + b6_2*cs->rhoV1[ip];
                cs->rhoE2[ip] = cs->rhoEk[ip] + b6_2*cs->rhoE1[ip];

	    }
    }


}


void LSLDDRK4::updateConservedData(){

    if(rkDualStepStep == 1){
	updateConservedDataStep1();
    }else if(rkDualStepStep == 2){
	updateConservedDataStep2();
    }

}

void LSLDDRK4::executeSolverLoop(){

    while(cs->endFlag == false){

	rkDualStepStep = 1;
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

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();


	//Step 4
	cs->rkStep = 4;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();


	//Step 5
	cs->rkStep = 5;
	cs->rkLast = true;	

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	cs->postStep();	

	//bust out early if we hit the end in the middle...
	if(cs->endFlag){
	    break;		
	}

	//Start the second RK step...

	rkDualStepStep = 2;
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

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();


	//Step 4
	cs->rkStep = 4;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 5
	cs->rkStep = 5;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 6
	cs->rkStep = 6;
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
