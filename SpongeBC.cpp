#include "SpongeBC.hpp"

SpongeBC::SpongeBC(AbstractSingleBlockMesh *msh, Domain *domain, IdealGas *idealGas, BC *bc, Options *opt){

	    
	    std::cout << endl;
	    std::cout << " > Sponge BC found, initializing Sponge average fields and strength fields..." << std::endl;
	

	    this->msh = msh;
	    this->domain = domain;
	    this->idealGas = idealGas;
	    this->bc = bc;
	    this->opt = opt;

	    this->Nx = domain->Nx;
	    this->Ny = domain->Ny;
	    N = Nx*Ny;

	    this->mpiRank = mpiRank;

	    sigma = new double[N];
	    spongeRhoAvg  = new double[N];
	    spongeRhoUAvg = new double[N];
	    spongeRhoVAvg = new double[N];
	    spongeRhoEAvg = new double[N];

	    avgT = opt->spongeAvgT;
	    epsP = 0.005;
	    spongeP = opt->spongeP;
	    spongeStrength = opt->spongeStrength;

	    //Need to initialize the sponge sigma to zero
	    FOR_XY sigma[ip] = 0.0;

	    //If rectangular sponge BC
	    if(opt->spongeKind == Options::RECTILINEAR){
	    	initRectSpongeBC();
	    //If cylindrical sponge BC
 	    }else if(opt->spongeKind == Options::CYLINDRICAL){
	    	initCylSpongeBC();	   
	    }else{
		//shouldn't get here for now...
	    }
	    //If spherical sponge BC


            getRange(sigma, "SPONGE SIGMA", Nx, Ny); 
	    std::cout << " > Done initializing sponge!" << std::endl;

	}


void SpongeBC::initRectSpongeBC(){


	    //Default the maximum ends of the sponge to the domain max, can and may be changed for curvilinear domains
	    double spongeXMin = msh->x_min;
	    double spongeXMax = msh->x_max;

	    double spongeYMin = msh->y_min;
	    double spongeYMax = msh->y_max;

	    //Default to an 1/8th of the domain size in that direction 
	    spongeLX0 = opt->spongeRectX0Perc*(spongeXMax-spongeXMin);
	    spongeLY0 = opt->spongeRectY0Perc*(spongeYMax-spongeYMin);

	    spongeLX1 = opt->spongeRectX1Perc*(spongeXMax-spongeXMin);
	    spongeLY1 = opt->spongeRectY1Perc*(spongeYMax-spongeYMin);



	    cout << " > spongeLX0 = " << spongeLX0 << endl;	    
	    cout << " > spongeLY0 = " << spongeLY0 << endl;	    
	    cout << " > spongeLX1 = " << spongeLX1 << endl;	    
	    cout << " > spongeLY1 = " << spongeLY1 << endl;	    
	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == Options::SPONGE){
		FOR_Y{
		    FOR_X{
		    	    int ip = GET2DINDEX_XY;
			    double dx = msh->x[ip]-spongeXMin;
		            if(dx < spongeLX0){
		        	double spongeX = (spongeLX0 - dx)/spongeLX0;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		    }
	  	}
	    }

	    if(bc->bcX1 == Options::SPONGE){
	        FOR_Y{
		    FOR_X{
			    int ip = GET2DINDEX_XY;
			    double dx = spongeXMax - msh->x[ip];
		    	    if(dx < spongeLX1){
		        	double spongeX = (msh->x[ip] - (spongeXMax - spongeLX1))/spongeLX1;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		    }
	  	}
	    }     

	    if(bc->bcY0 == Options::SPONGE){
	        FOR_Y{
	            FOR_X{
		            int ip = GET2DINDEX_XY;
			    double dy = msh->y[ip]-spongeYMin;	
		            if(dy < spongeLY0){
		                double spongeY = (spongeLY0 - dy)/spongeLY0;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		    }
		}
	    }
	
	    if(bc->bcY1 == Options::SPONGE){
		FOR_X{
		    FOR_Y{
		            int ip = GET2DINDEX_XY;
			    double dy = spongeYMax-msh->y[ip];
		            if(dy < spongeLY1){
		                double spongeY = (msh->y[ip] - (spongeYMax - spongeLY1))/spongeLY1;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		    }
		}
	    }    

}

void SpongeBC::initCylSpongeBC(){



        double spongeXMin = msh->x_min[0];
        double spongeXMax = msh->x_max[0];
	double spongeYMin = msh->x_min[1];
	double spongeYMax = msh->x_max[1];

	double spongeCylAxis[3];
	spongeCylAxis[2] = opt->spongeCylAxisZ;


	double spongeRmin = opt->rMin;

	double spongeRmax;

	if(opt->spongeCylAxisOrient == 2){

	     //This maybe isn't the most general formulation, just takes the most distant boundary point  
	     //from the cylindrical axis as the the maximum sponge radius
	     double spongeRXmax1 = fabs(spongeXMin - spongeCylAxis[0]);
	     double spongeRXmax2 = fabs(spongeXMax - spongeCylAxis[0]);
	     double spongeRXmax  = fmax(spongeRXmax1, spongeRXmax2);

	     double spongeRYmax1 = fabs(spongeYMin - spongeCylAxis[1]);
	     double spongeRYmax2 = fabs(spongeYMax - spongeCylAxis[1]);
	     double spongeRYmax  = fmax(spongeRYmax1, spongeRYmax2);

	     spongeRmax = fmax(spongeRXmax, spongeRYmax);

	     //Is it possible for other faces to be sponges in this kind of configuration? I think so?! Need to implement
	     if(bc->bcY1 == Options::SPONGE){
		 FOR_X{
		     FOR_Y{
			     int ip = GET2DINDEX_XY;
			     double r = sqrt(pow(msh->x[ip]-spongeCylAxis[0],2.0) + pow(msh->y[ip]-spongeCylAxis[1],2.0));
			     if(r > spongeRmin){
				 double spongeR = (r-spongeRmin)/(spongeRmax-spongeRmin);
				 sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeR,2.0) + 0.845*pow(spongeR, 8.0)), sigma[ip]);
			     }
		     }
	 	 }
	     }


	
	}else{
	
	    cout << "SPONGECYLINDRICAL SHAPE NEEDS TO BE IMPLEMENTED FOR ORIENT = " << opt->spongeCylAxisOrient << endl;
	    abort();
	}

	    //Default the maximum ends of the sponge to the domain max, can and may be changed for curvilinear domains

/*		
	    if(bc->bcY1 == Options::SPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
			FOR_Z_YPEN{
			    int ip = GETMAJIND_YPEN;
			    double r = sqrt(pow(msh->x[ip]-spongeXAvg,2.0) + pow(msh->y[ip]-spongeYAvg,2.0));
			    if(r > spongeYRmin){
				double spongeR = (r-spongeYRmin)/(spongeYRmax-spongeYRmin);
				sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeR,2.0) + 0.845*pow(spongeR, 8.0)), sigma[ip]);
			    }
			}
		    }
		}
	    }

	    //Default to an 1/8th of the domain size in that direction 
	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
			FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dx = msh->x[ip]-spongeXMin;
		            if(dx < spongeLX){
		        	double spongeX = (spongeLX - dx)/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }

	    if(bc->bcX1 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
			    int ip = GETMAJIND_YPEN;
			    double dx = spongeXMax - msh->x[ip];
		    	    if(dx < spongeLX){
		        	double spongeX = (msh->x[ip] - (spongeXMax - spongeLX))/spongeLX;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ip]);
			    }
		        }
		    }
	  	}
	    }     

	    if(bc->bcY0 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
			FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = msh->y[ip]-spongeYMin;	
		            if(dy < spongeLY){
		                double spongeY = (spongeLY - dy)/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcY1 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dy = spongeYMax-msh->y[ip];
		            if(dy < spongeLY){
		                double spongeY = (msh->y[ip] - (spongeYMax - spongeLY))/spongeLY;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcZ0 == Options::RECT_CURVILINEARSPONGE){
		FOR_X_YPEN{
		    FOR_Y_YPEN{
		        FOR_Z_YPEN{
		            int ip = GETMAJIND_YPEN;
			    double dz = msh->z[ip] - spongeZMin;
		            if(dz < spongeLZ){
		                double spongeZ = (spongeLZ - dz)/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
  	    }
	
	    if(bc->bcZ1 == Options::RECT_CURVILINEARSPONGE){
	        FOR_X_YPEN{
	            FOR_Y_YPEN{
		        FOR_Z_YPEN{
		    	    int ip = GETMAJIND_YPEN;
			    double dz = spongeZMax-msh->z[ip];
		    	    if(dz < spongeLZ){
		        	double spongeZ = (msh->z[ip] - (spongeZMax-spongeLZ))/spongeLZ;
			        sigma[ip] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ip]);
			    }
		        }
		    }
		}
	    }    
*/

}

void SpongeBC::dumpSpongeAvg(int timeStep){
/*
    double *tempY1, *tempY2, *tempY3, *tempY4, *tempY5;

    tempY1 = new double[N];  
    tempY2 = new double[N];  
    tempY3 = new double[N];  
    tempY4 = new double[N];  


    //Move stuff over to x-major y-pencils for writing out
    
	FOR_Y_YPEN{
	    FOR_X_YPEN{
		int ip = GETMAJIND_YPEN;
		int jp = GETIND_YPEN;

		tempY1[jp] = spongeRhoAvg[ip]; 		
		tempY2[jp] = spongeRhoUAvg[ip]; 		
		tempY3[jp] = spongeRhoVAvg[ip]; 		
		tempY4[jp] = spongeRhoEAvg[ip]; 		
	    }
	}

    ofstream outfile;
    outfile.precision(17);
    string outputFileName;
    outputFileName = "SpongeAvgDump.";

    ostringstream timeStepString;
    timeStepString << timeStep;  

    outputFileName.append(timeStepString.str());

    MPI_File fh;
    MPI_Offset disp, filesize;

    MPI_File_open(MPI_COMM_WORLD, outputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    filesize = 0;
    MPI_File_set_size(fh, filesize);
  
    disp = 0;

    c2d->writeVar(fh, disp, 1, tempY1); 
    c2d->writeVar(fh, disp, 1, tempY2); 
    c2d->writeVar(fh, disp, 1, tempY3); 
    c2d->writeVar(fh, disp, 1, tempY4); 
    c2d->writeVar(fh, disp, 1, tempY5); 

    MPI_File_close(&fh);

    c2d->deallocXYZ(tempY1);
    c2d->deallocXYZ(tempY2);
    c2d->deallocXYZ(tempY3);
    c2d->deallocXYZ(tempY4);
    c2d->deallocXYZ(tempY5);
*/
}

void SpongeBC::readSpongeAvgFromRestart(){
/*
    string filename = opt->sponge_filename;

    MPI_File fh;
    MPI_Offset disp, filesize;

    IF_RANK0{
	cout << " > Initializing sponge average field from restart file: " << filename << endl;
    }

    int err = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    double *rho_in,
           *rhoU_in, 
           *rhoV_in,
           *rhoW_in,
           *rhoE_in;

    c2d->allocY(rho_in);
    c2d->allocY(rhoU_in);
    c2d->allocY(rhoV_in);
    c2d->allocY(rhoW_in);
    c2d->allocY(rhoE_in);

    disp = 0;

    c2d->readVar(fh, disp, 1, rho_in);
    c2d->readVar(fh, disp, 1, rhoU_in);
    c2d->readVar(fh, disp, 1, rhoV_in);
    c2d->readVar(fh, disp, 1, rhoW_in);
    c2d->readVar(fh, disp, 1, rhoE_in);

    MPI_File_close(&fh);

    FOR_Z_YPEN{
	FOR_Y_YPEN{
	    FOR_X_YPEN{
		int ip = GETMAJIND_YPEN;
		int jp = GETIND_YPEN;

		spongeRhoAvg[ip]  = rho_in[jp];
		spongeRhoUAvg[ip] = rhoU_in[jp];
		spongeRhoVAvg[ip] = rhoV_in[jp];
		spongeRhoWAvg[ip] = rhoW_in[jp];
		spongeRhoEAvg[ip] = rhoE_in[jp];
	    }
	}
    }


    c2d->deallocXYZ(rho_in);
    c2d->deallocXYZ(rhoU_in);
    c2d->deallocXYZ(rhoV_in);
    c2d->deallocXYZ(rhoW_in);
    c2d->deallocXYZ(rhoE_in);

    IF_RANK0 cout << " > Sponge Average Conditions initialized from file! " << endl;
    
*/
}

