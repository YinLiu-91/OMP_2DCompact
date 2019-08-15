#ifndef _BCH_
#define _BCH_

#include <iostream>
#include "Macros.hpp"
#include "Options.hpp"

class BC{

    public:

	Options::BCType bcXType, bcYType;
	Options::BCKind bcX0, bcX1, bcY0, bcY1;

	double periodicDisp[2][2];

	BC(Options *opt, bool bcPeriodic[2]){

	    this->bcXType = opt->bcXType;
	    this->bcYType = opt->bcYType;

	    this->bcX0 = opt->bcX0;
	    this->bcX1 = opt->bcX1;
	    this->bcY0 = opt->bcY0;
	    this->bcY1 = opt->bcY1;

	    for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
	            this->periodicDisp[i][j] = opt->periodicDisp[i][j];
	        }
	    }

	    std::cout << std::endl;
	    std::cout << " > Initializing boundary conditions..." << std::endl;
	
	    std::cout << " >     X BOUNDARY CONDITIONS    " << std::endl;
	    if(bcXType == Options::PERIODIC_SOLVE){
		std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[0] = true;
	    }else{
		if(bcX0 == Options::SPONGE){
		   std::cout << " > X0=SPG";
		}else if(bcX0 == Options::ADIABATIC_WALL){
		   std::cout << " > X0=ABW";
		}else if(bcX0 == Options::MOVING_ADIABATIC_WALL){
		   std::cout << " > X0=MOW";
		}

		if(bcX1 == Options::SPONGE){
		   std::cout << "----------------SPG=X1" << std::endl;
		}else if(bcX1 == Options::ADIABATIC_WALL){
		   std::cout << "----------------ABW=X1" << std::endl;
		}else if(bcX1 == Options::MOVING_ADIABATIC_WALL){
		   std::cout << "----------------MOW=X1" << std::endl;
		}

		bcPeriodic[0] = false;
	    }

	    std::cout << " >     Y BOUNDARY CONDITIONS    " << std::endl;
	    if(bcYType == Options::PERIODIC_SOLVE){
	        std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[1] = true;
	    }else{
		if(bcY0 == Options::SPONGE){
		   std::cout << " > Y0=SPG";
		}else if(bcY0 == Options::ADIABATIC_WALL){
		   std::cout << " > Y0=ABW";
		}else if(bcY0 == Options::MOVING_ADIABATIC_WALL){
		   std::cout << " > Y0=MOW";
		}

		if(bcY1 == Options::SPONGE){
		   std::cout << "----------------SPG=Y1" << std::endl;
		}else if(bcY1 == Options::ADIABATIC_WALL){
		   std::cout << "----------------ABW=Y1" << std::endl;
		}else if(bcY1 == Options::MOVING_ADIABATIC_WALL){
		    std::cout << "----------------MOW=Y1" << std::endl;
		}

		bcPeriodic[1] = false;
	    }

		std::cout << std::endl;

	
	     //Clear out the displacement values if not periodic
	     if(bcPeriodic[0] == false){
		periodicDisp[0][0] = 0.0;
		periodicDisp[0][1] = 0.0;
	     }

	     if(bcPeriodic[1] == false){
		periodicDisp[1][0] = 0.0;
		periodicDisp[1][1] = 0.0;
	     }

	}

};

#endif
