#include "AbstractSingleBlockMesh.hpp"


void AbstractSingleBlockMesh::allocateForMesh(){

        //Doing this for base y-pencil solvers...
        x = new double[Nx*Ny];
        y = new double[Nx*Ny];
        J = new double[Nx*Ny];
        J11 = new double[Nx*Ny];
        J12 = new double[Nx*Ny];
        J21 = new double[Nx*Ny];
        J22 = new double[Nx*Ny];

}

void AbstractSingleBlockMesh::solveForJacobians(){


	 cout << " > Solving for Jacobian Matrix... ";

	double *dxde1 = new double[Nx*Ny];
	double *dxde2 = new double[Nx*Ny];
	double *dyde1 = new double[Nx*Ny];
	double *dyde2 = new double[Nx*Ny];

        double *x_t, *y_t, *dxde2_t, *dyde2_t;
	x_t = new double[Nx*Ny];
	y_t = new double[Nx*Ny];
	dxde2_t = new double[Nx*Ny];
	dyde2_t = new double[Nx*Ny];

	transposeMatrix_Fast2(x, Nx, Ny, x_t, opt->blocksize);
	transposeMatrix_Fast2(y, Nx, Ny, y_t, opt->blocksize);


	//Do the y derivatives first
	if(periodicBCY){

	    double *Nm3x, *Nm2x, *Nm1x, *Np1x, *Np2x, *Np3x; 
	    Nm3x = new double[Nx];
	    Nm2x = new double[Nx];
	    Nm1x = new double[Nx];
	    Np1x = new double[Nx];
	    Np2x = new double[Nx];
	    Np3x = new double[Nx];

	    double *Nm3y, *Nm2y, *Nm1y, *Np1y, *Np2y, *Np3y;
	    Nm3y = new double[Nx];
	    Nm2y = new double[Nx];
	    Nm1y = new double[Nx];
	    Np1y = new double[Nx];
	    Np2y = new double[Nx];
	    Np3y = new double[Nx];

	    FOR_X{
		    int iim3 = (Ny-3)*Nx + i; 
		    int iim2 = (Ny-2)*Nx + i;		
		    int iim1 = (Ny-1)*Nx + i;		
		    int iip1 = (0)*Nx + i;		
		    int iip2 = (1)*Nx + i;		
		    int iip3 = (2)*Nx + i;		

		    Nm3x[i] = x[iim3]-periodicYTranslation[0];
		    Nm2x[i] = x[iim2]-periodicYTranslation[0];
		    Nm1x[i] = x[iim1]-periodicYTranslation[0];
		    Np1x[i] = x[iip1]+periodicYTranslation[0];
		    Np2x[i] = x[iip2]+periodicYTranslation[0];
		    Np3x[i] = x[iip3]+periodicYTranslation[0];
	
		    Nm3y[i] = y[iim3]-periodicYTranslation[1];
		    Nm2y[i] = y[iim2]-periodicYTranslation[1];
		    Nm1y[i] = y[iim1]-periodicYTranslation[1];
		    Np1y[i] = y[iip1]+periodicYTranslation[1];
		    Np2y[i] = y[iip2]+periodicYTranslation[1];
		    Np3y[i] = y[iip3]+periodicYTranslation[1];
	
	    }

	    if(derivY->rhsBandwidth == AbstractDerivatives::BW5){

	        derivY->calc1stDerivField_TPB(x_t, dxde2_t, Nm2x, Nm1x, Np1x, Np2x);
	        derivY->calc1stDerivField_TPB(y_t, dyde2_t, Nm2y, Nm1y, Np1y, Np2y);

	    }else if(derivY->rhsBandwidth == AbstractDerivatives::BW7){
	        derivY->calc1stDerivField_TPB(x_t, dxde2_t, Nm3x, Nm2x, Nm1x, Np1x, Np2x, Np3x);
	        derivY->calc1stDerivField_TPB(y_t, dyde2_t, Nm3y, Nm2y, Nm1y, Np1y, Np2y, Np3y);
	    }else{
		cout << "ERROR: Unknown derivY rhs bandwidth" << endl;
		abort();
	    }

	    delete[] Nm3x;
	    delete[] Nm2x;
	    delete[] Nm1x;
	    delete[] Np1x;
	    delete[] Np2x;
	    delete[] Np3x;

	    delete[] Nm3y;
	    delete[] Nm2y;
	    delete[] Nm1y;
	    delete[] Np1y;
	    delete[] Np2y;
	    delete[] Np3y;



	}else{
	    derivY->calc1stDerivField(x_t, dxde2_t);
	    derivY->calc1stDerivField(y_t, dyde2_t);
	}

	transposeMatrix_Fast2(dxde2_t, Ny, Nx, dxde2, opt->blocksize);
	transposeMatrix_Fast2(dyde2_t, Ny, Nx, dyde2, opt->blocksize);



	if(periodicBCX){

    	    double *Nm3x, *Nm2x, *Nm1x, *Np1x, *Np2x, *Np3x;
            Nm3x = new double[Ny];
            Nm2x = new double[Ny];
            Nm1x = new double[Ny];
            Np1x = new double[Ny];
            Np2x = new double[Ny];
            Np3x = new double[Ny];

	    double *Nm3y, *Nm2y, *Nm1y, *Np1y, *Np2y, *Np3y;
            Nm3y = new double[Ny];
            Nm2y = new double[Ny];
            Nm1y = new double[Ny];
            Np1y = new double[Ny];
            Np2y = new double[Ny];
            Np3y = new double[Ny];


            FOR_Y{
                    int iim3 = j*Nx + Nx - 3;
                    int iim2 = j*Nx + Nx - 2;
                    int iim1 = j*Nx + Nx - 1;
                    int iip1 = j*Nx + 0;
                    int iip2 = j*Nx + 1;
                    int iip3 = j*Nx + 2;

                    Nm3x[j] = x[iim3]-periodicXTranslation[0];
                    Nm2x[j] = x[iim2]-periodicXTranslation[0];
                    Nm1x[j] = x[iim1]-periodicXTranslation[0];
                    Np1x[j] = x[iip1]+periodicXTranslation[0];
                    Np2x[j] = x[iip2]+periodicXTranslation[0];
                    Np3x[j] = x[iip3]+periodicXTranslation[0];

                    Nm3y[j] = y[iim3]-periodicXTranslation[1];
                    Nm2y[j] = y[iim2]-periodicXTranslation[1];
                    Nm1y[j] = y[iim1]-periodicXTranslation[1];
                    Np1y[j] = y[iip1]+periodicXTranslation[1];
                    Np2y[j] = y[iip2]+periodicXTranslation[1];
                    Np3y[j] = y[iip3]+periodicXTranslation[1];

            }

	    if(derivX->rhsBandwidth == AbstractDerivatives::BW5){
                derivX->calc1stDerivField_TPB(x, dxde1, Nm2x, Nm1x, Np1x, Np2x);
                derivX->calc1stDerivField_TPB(y, dyde1, Nm2y, Nm1y, Np1y, Np2y);
	    }else if(derivX->rhsBandwidth == AbstractDerivatives::BW7){
                derivX->calc1stDerivField_TPB(x, dxde1, Nm3x, Nm2x, Nm1x, Np1x, Np2x, Np3x);
                derivX->calc1stDerivField_TPB(y, dyde1, Nm3y, Nm2y, Nm1y, Np1y, Np2y, Np3y);
	    }else{
		cout << "ERROR: Unknown derivY rhs bandwidth" << endl;
		abort();
	    }


            delete[] Nm3x;
            delete[] Nm2x;
            delete[] Nm1x;
            delete[] Np1x;
            delete[] Np2x;
            delete[] Np3x;

            delete[] Nm3y;
            delete[] Nm2y;
            delete[] Nm1y;
            delete[] Np1y;
            delete[] Np2y;
            delete[] Np3y;

	}else{
	    derivX->calc1stDerivField(x, dxde1);
	    derivX->calc1stDerivField(y, dyde1);
	}

	//For some reason, in optimized compilers this loop can break things...
	#pragma omp parallel for
	FOR_XY{
	    double Jt   =  dxde1[ip]*dyde2[ip] - dxde2[ip]*dyde1[ip];
	    J11[ip] =  dyde2[ip]/Jt;
	    J12[ip] = -dxde2[ip]/Jt;
	    J21[ip] = -dyde1[ip]/Jt;
	    J22[ip] =  dxde1[ip]/Jt;
	    J[ip]   =  1.0/Jt; 	
  	}

	//Copy data over to solver to make simpler to access
	cs->J = J;
	cs->J11 = J11;
	cs->J12 = J12;
	cs->J21 = J21;
	cs->J22 = J22;

	 cout << "done!" << endl;

	getRange(dxde1, "dxde1", Nx, Ny);	
	getRange(dyde1, "dyde1", Nx, Ny);	
	getRange(dxde2, "dxde2", Nx, Ny);	
	getRange(dyde2, "dyde2", Nx, Ny);	
        getRange(J,     "J", Nx, Ny);
        getRange(J11, "J11", Nx, Ny);
        getRange(J12, "J12", Nx, Ny);
        getRange(J21, "J21", Nx, Ny);
        getRange(J22, "J22", Nx, Ny);
	
    delete[] x_t;
    delete[] y_t;
    delete[] dxde1;
    delete[] dxde2;
    delete[] dyde1;
    delete[] dyde2;

    delete[] dxde2_t; 
    delete[] dyde2_t; 

}

//void AbstractSingleBlockMesh::dumpGrid(){
/*
	double time1 = MPI_Wtime();

	{
	    cout << endl;
	    cout << " > Dumping grid " << endl;
	}

	double *x_temp,
	       *y_temp,
	       *z_temp;

	c2d->allocY(x_temp);
	c2d->allocY(y_temp);
	c2d->allocY(z_temp);

	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{
		    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

		    x_temp[jp] = x[ip];
		    y_temp[jp] = y[ip];
		    z_temp[jp] = z[ip];
		}
	    }
	}

	ofstream outfile;
	outfile.precision(17);
	string outputFileName;
	outputFileName = "GridDump.XYZ";

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;
	c2d->writeVar(fh, disp, 1, x_temp); 
	c2d->writeVar(fh, disp, 1, y_temp); 
	c2d->writeVar(fh, disp, 1, z_temp); 

	MPI_File_close(&fh);

	c2d->deallocXYZ(x_temp);
	c2d->deallocXYZ(y_temp);
	c2d->deallocXYZ(z_temp);

	double time2 = MPI_Wtime();

	{
	    cout << " > Grid dump took " << time2-time1 << " seconds." << endl;   
	}
*/
//}

//void AbstractSingleBlockMesh::getOrderedBlockXiCoordinates(int ip, int jp, int kp, double box_pxi[8][3]){
/*
        double dxi, deta, dzta;

        //This is correct if the domain bounds have been properly set to their max xi, eta, and zeta values in the
        //MPI_Compact.cpp file
        dxi  = d->dx;
        deta = d->dy;
        dzta = d->dz;

        box_pxi[0][0] = dxi *(double)ip;
        box_pxi[1][0] = dxi *(double)ip;
        box_pxi[2][0] = dxi *(double)ip;
        box_pxi[3][0] = dxi *(double)ip;
        box_pxi[4][0] = dxi *(double)(ip+1);
        box_pxi[5][0] = dxi *(double)(ip+1);
        box_pxi[6][0] = dxi *(double)(ip+1);
        box_pxi[7][0] = dxi *(double)(ip+1);

        box_pxi[0][1] = deta*(double)jp;
        box_pxi[1][1] = deta*(double)jp;
        box_pxi[2][1] = deta*(double)(jp+1);
        box_pxi[3][1] = deta*(double)(jp+1);
        box_pxi[4][1] = deta*(double)jp;
        box_pxi[5][1] = deta*(double)jp;
        box_pxi[6][1] = deta*(double)(jp+1);
        box_pxi[7][1] = deta*(double)(jp+1);

        box_pxi[0][2] = dzta*(double)kp;
        box_pxi[1][2] = dzta*(double)(kp+1);
        box_pxi[2][2] = dzta*(double)kp;
        box_pxi[3][2] = dzta*(double)(kp+1);
        box_pxi[4][2] = dzta*(double)kp;
        box_pxi[5][2] = dzta*(double)(kp+1);
        box_pxi[6][2] = dzta*(double)kp;
        box_pxi[7][2] = dzta*(double)(kp+1);
*/
//};


int AbstractSingleBlockMesh::findCVForPoint(double p[2]){

    int cvListSize, cvList[8192];

    adt->buildListForPoint(cvListSize, cvList, p);

    //This function will return -1 for a point not found and the local control volume index in y-major
    //indexing if it is found
    int base_index = -1;
    for(int ii = 0; ii < cvListSize; ii++){

	base_index = cvList[ii];

	//Back out the ip, jp, kp coordinates based off single major index int
	int i =  base_index%Nx;
	int j =  (base_index-i)/Nx;

	int ip = base_index;

        if((i == Nx-1 && !periodicBCX) ||
	   (j == Ny-1 && !periodicBCY)){
	     break;
        }

	double a[2], b[2], c[2], d[2];

	a[0] = x[ip];
	a[1] = y[ip];

	int iii;
	if(i == Nx-1 && periodicBCX){
	  iii = j*Nx + 0;
	  b[0] = x[iii] + periodicXTranslation[0];
	  b[1] = y[iii] + periodicXTranslation[1];
	}else{
	  iii = j*Nx + i + 1;
	  b[0] = x[iii];
	  b[1] = y[iii];
	}

	if(j == Ny-1 && periodicBCY){
	  iii =   (0)*Nx + i;
	  c[0] = x[iii] + periodicYTranslation[0];
	  c[1] = y[iii] + periodicYTranslation[1];
	}else{
    	  iii = (j+1)*Nx + i;
	  c[0] = x[iii];
	  c[1] = y[iii];
	}

	if((j==Ny-1) && (i==Nx-1) && periodicBCY && periodicBCX){
	  iii = 0;
	  d[0] = x[iii]+periodicXTranslation[0]+periodicYTranslation[0];
	  d[1] = y[iii]+periodicXTranslation[1]+periodicYTranslation[1];
	}else if(i==Nx-1 && periodicBCX){
	  iii = (j+1)*Nx + 0;
	  d[0] = x[iii]+periodicXTranslation[0];
	  d[1] = y[iii]+periodicXTranslation[1];
	}else if(j==Ny-1 && periodicBCY){
	  iii =  (0)*Nx + i + 1;
	  d[0] = x[iii]+periodicYTranslation[0];
	  d[1] = y[iii]+periodicYTranslation[1];
	}else{
	  iii = (j+1)*Nx + i + 1;
	  d[0] = x[iii];
	  d[1] = y[iii];
	}

	if(isPointInSquare(p, a, b, c, d)){
	    base_index = cvList[ii];
	    break;
	}else{
	    base_index = -1;
	}
    }

    return base_index;    

}

void AbstractSingleBlockMesh::initMeshADT(){

	    //Initialize the ADT object for interpolation and locating points in the grid...
	     cout << " > Initializing ADT..." << endl;

            //Get our coordinates in neighbors across partitions using halo updates
  
            double (*boundBoxMin)[2] = new double[Nx*Ny][2];
            double (*boundBoxMax)[2] = new double[Nx*Ny][2];

	    FOR_Y{
		FOR_X{
		    int ip = GET2DINDEX_XY;

		    bool noBBFlag = false;	
                    if((i == Nx-1 && !periodicBCX) ||
		       (j == Ny-1 && !periodicBCY)){
                    	noBBFlag = true;
                    }

		    double x_max = -1.0e100; 
	   	    double y_max = -1.0e100; 
		    double x_min =  1.0e100;
		    double y_min =  1.0e100;
		
		    if(!noBBFlag){
			double xp[4][2];

			xp[0][0] = x[ip];
			xp[0][1] = y[ip];

			int ii;

			if(i == Nx-1 && periodicBCX){
			  ii = j*Nx + 0;
			  xp[1][0] = x[ii] + periodicXTranslation[0];
			  xp[1][1] = y[ii] + periodicXTranslation[1];
			}else{
			  ii = j*Nx + i + 1;
			  xp[1][0] = x[ii];
			  xp[1][1] = y[ii];
			}

			if(j == Ny-1 && periodicBCY){
			  ii =   (0)*Nx + i;
			  xp[2][0] = x[ii] + periodicYTranslation[0];
			  xp[2][1] = y[ii] + periodicYTranslation[1];
			}else{
    			  ii = (j+1)*Nx + i;
			  xp[2][0] = x[ii];
			  xp[2][1] = y[ii];
			}

			if((j==Ny-1) && (i==Nx-1) && periodicBCY && periodicBCX){
			  ii = 0;
			  xp[3][0] = x[ii]+periodicXTranslation[0]+periodicYTranslation[0];
			  xp[3][1] = y[ii]+periodicXTranslation[1]+periodicYTranslation[1];
			}else if(i==Nx-1 && periodicBCX){
			  ii = (j+1)*Nx + 0;
			  xp[3][0] = x[ii]+periodicXTranslation[0];
			  xp[3][1] = y[ii]+periodicXTranslation[1];
			}else if(j==Ny-1 && periodicBCY){
			  ii =   (0)*Nx + i + 1;
			  xp[3][0] = x[ii]+periodicYTranslation[0];
			  xp[3][1] = y[ii]+periodicYTranslation[1];
			}else{
			  ii = (j+1)*Nx + i + 1;
			  xp[3][0] = x[ii];
			  xp[3][1] = y[ii];
			}

			for(int iip = 0; iip < 4; iip++){
			    x_max = fmax(x_max, xp[iip][0]);
			    x_min = fmin(x_min, xp[iip][0]);

			    y_max = fmax(y_max, xp[iip][1]);
			    y_min = fmin(y_min, xp[iip][1]);
			}

		    }else{
			x_max = x[ip];	
			x_min = x[ip];	
			y_max = y[ip];
			y_min = y[ip];
		    }

		    boundBoxMin[ip][0] = x_min;
		    boundBoxMin[ip][1] = y_min;
		    boundBoxMax[ip][0] = x_max;
		    boundBoxMax[ip][1] = y_max;
		}
	    }

	    FOR_XY{
		FOR_I2{
		    double delta = 1.0E-6*(boundBoxMax[ip][i] - boundBoxMin[ip][i]);
		    boundBoxMax[ip][i] += delta;
		    boundBoxMin[ip][i] -= delta;
		}
	    }

	     cout << " > Done getting bounding boxes for the CV's, initializing ADT... " << endl;

	    adt = new Adt<double>(Nx*Ny, boundBoxMin, boundBoxMax);

	     cout << " > Done!" << endl;


	    //Dump the solution grid
	    //dumpGrid();

	    delete[] boundBoxMin;
	    delete[] boundBoxMax;


}


void AbstractSingleBlockMesh::handlePeriodicStuff(){
	    //Initiallize all of the flags to false
	    periodicBCX = false;	    
	    periodicBCY = false;	    
	    transPeriodicX = false;
	    transPeriodicY = false;
	    interPeriodicX = false;
	    interPeriodicY = false;

	    //Handle periodic boundary condition flags and displacements
	    if(cs->bc->bcXType == Options::PERIODIC_SOLVE){
		periodicBCX = true;
		if(cs->bc->bcY0 == Options::PERIODIC && cs->bc->bcY1 == Options::PERIODIC){
		    transPeriodicX = true;
		    periodicXTranslation[0] = cs->bc->periodicDisp[0][0];
		    periodicXTranslation[1] = cs->bc->periodicDisp[0][1];
		     cout << "  Periodic x-face translation = {" << periodicXTranslation[0] << ", " << periodicXTranslation[1] <<  "}" << endl;;
		}else if(cs->bc->bcX0 == Options::INTERNALLY_PERIODIC && cs->bc->bcX1 == Options::INTERNALLY_PERIODIC){
		    interPeriodicX = true;
		    periodicXTranslation[0] = 0.0;
		    periodicXTranslation[1] = 0.0;
		}else{
		     cout << "  Incompatable periodic bc inputs in x-direction! " << endl;
		}
	    }

	    if(cs->bc->bcYType == Options::PERIODIC_SOLVE){
		periodicBCY = true;
		if(cs->bc->bcY0 == Options::PERIODIC && cs->bc->bcY1 == Options::PERIODIC){
		    transPeriodicY = true;
		    periodicYTranslation[0] = cs->bc->periodicDisp[1][0]; 
		    periodicYTranslation[1] = cs->bc->periodicDisp[1][1];
		     cout << "  Periodic y-face translation = {" << periodicYTranslation[0] << ", " << periodicYTranslation[1] << "}" << endl;
		}else if(cs->bc->bcY0 == Options::INTERNALLY_PERIODIC && cs->bc->bcY1 == Options::INTERNALLY_PERIODIC){
		    interPeriodicY = true;			
		    periodicYTranslation[0] = 0.0; 
		    periodicYTranslation[1] = 0.0;
		}else{
		     cout << "  Incompatable periodic bc inputs in y-direction! " << endl;
		} 
	    }

}
