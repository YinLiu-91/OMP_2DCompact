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
	    derivX->calc1stDerivField(y, dxde1);
	}

	#pragma omp parallel for
	FOR_XY{
	    J[ip]   =  dxde1[ip]*dyde2[ip] - dxde2[ip]*dyde1[ip];
	    J11[ip] =  dyde2[ip]/J[ip];
	    J12[ip] = -dxde2[ip]/J[ip];
	    J21[ip] = -dyde1[ip]/J[ip];
	    J22[ip] =  dxde1[ip]/J[ip];
	    J[ip]   =  1.0/J[ip]; 	
  	}

	//Copy data over to solver to make simpler to access
	cs->J = J;
	cs->J11 = J11;
	cs->J12 = J12;
	cs->J21 = J21;
	cs->J22 = J22;

	 cout << "done!" << endl;
	
        getRange(J,     "J", Nx, Ny);
        getRange(J11, "J11", Nx, Ny);
        getRange(J12, "J12", Nx, Ny);
        getRange(J21, "J21", Nx, Ny);
        getRange(J22, "J22", Nx, Ny);
	
/*
	 cout << " > Checking the values of the metric indentities, values should be small" << endl;

	double *I1_1, *I1_2, *I1_3;
	double *I2_1, *I2_2, *I2_3;
	double *I3_1, *I3_2, *I3_3;

	c2d->allocY(I1_1);
	c2d->allocY(I1_2);
	c2d->allocY(I1_3);
	c2d->allocY(I2_1);
	c2d->allocY(I2_2);
	c2d->allocY(I2_3);
	c2d->allocY(I3_1);
	c2d->allocY(I3_2);
	c2d->allocY(I3_3);

	//In calculating these identities does there need to be a transformation at the periodic boundaries??

	derivY->calc1stDerivField(J21, I1_2);
	derivY->calc1stDerivField(J22, I2_2);
	derivY->calc1stDerivField(J23, I3_2);


        c2d->transposeY2X_MajorIndex(J11, tempX1);
        c2d->transposeY2X_MajorIndex(J12, tempX2);
        c2d->transposeY2X_MajorIndex(J13, tempX3);

	derivX->calc1stDerivField(tempX1, tempX4);
	derivX->calc1stDerivField(tempX2, tempX5);
	derivX->calc1stDerivField(tempX3, tempX6);

	c2d->transposeX2Y_MajorIndex(tempX4, I1_1);
	c2d->transposeX2Y_MajorIndex(tempX5, I2_1);
	c2d->transposeX2Y_MajorIndex(tempX6, I3_1);


        c2d->transposeY2Z_MajorIndex(J31, tempZ1);
        c2d->transposeY2Z_MajorIndex(J32, tempZ2);
        c2d->transposeY2Z_MajorIndex(J33, tempZ3);

	derivZ->calc1stDerivField(tempZ1, tempZ4);
	derivZ->calc1stDerivField(tempZ2, tempZ5);
	derivZ->calc1stDerivField(tempZ3, tempZ6);

	c2d->transposeZ2Y_MajorIndex(tempZ4, I1_3);
	c2d->transposeZ2Y_MajorIndex(tempZ5, I2_3);
	c2d->transposeZ2Y_MajorIndex(tempZ6, I3_3);

	FOR_XYZ_YPEN{
	    I1_1[ip] = I1_1[ip] + I1_2[ip] + I1_3[ip];
	    I2_1[ip] = I2_1[ip] + I2_2[ip] + I2_3[ip];
	    I3_1[ip] = I3_1[ip] + I3_2[ip] + I3_3[ip];
	}

	getRange(I1_1, "Metric Identity 1", pySize[0], pySize[1], pySize[2], mpiRank);
	getRange(I2_1, "Metric Identity 2", pySize[0], pySize[1], pySize[2], mpiRank);
	getRange(I3_1, "Metric Identity 3", pySize[0], pySize[1], pySize[2], mpiRank);

	c2d->deallocXYZ(I1_1);
	c2d->deallocXYZ(I2_1);
	c2d->deallocXYZ(I3_1);
	c2d->deallocXYZ(I1_2);
	c2d->deallocXYZ(I2_2);
	c2d->deallocXYZ(I3_2);
	c2d->deallocXYZ(I1_3);
	c2d->deallocXYZ(I2_3);
	c2d->deallocXYZ(I3_3);

	//Free up all of the spaces we've been using...
	c2d->deallocXYZ(tempX1);
	c2d->deallocXYZ(tempX2);
	c2d->deallocXYZ(tempX3);
	c2d->deallocXYZ(tempX4);
	c2d->deallocXYZ(tempX5);
	c2d->deallocXYZ(tempX6);
	c2d->deallocXYZ(tempX7);
	c2d->deallocXYZ(tempX8);
	c2d->deallocXYZ(tempX9);
	c2d->deallocXYZ(tempX10);
	c2d->deallocXYZ(tempX11);
	c2d->deallocXYZ(tempX12);
	c2d->deallocXYZ(tempZ1);
	c2d->deallocXYZ(tempZ2);
	c2d->deallocXYZ(tempZ3);
	c2d->deallocXYZ(tempZ4);
	c2d->deallocXYZ(tempZ5);
	c2d->deallocXYZ(tempZ6);
	c2d->deallocXYZ(tempZ7);
	c2d->deallocXYZ(tempZ8);
	c2d->deallocXYZ(tempZ9);
	c2d->deallocXYZ(tempZ10);
	c2d->deallocXYZ(tempZ11);
	c2d->deallocXYZ(tempZ12);

	c2d->deallocXYZ(preJdet1);
	c2d->deallocXYZ(preJdet2);
	c2d->deallocXYZ(preJdet3);

	c2d->deallocXYZ(Jdet1);
	c2d->deallocXYZ(Jdet2);
	c2d->deallocXYZ(Jdet3);
*/


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

//void AbstractSingleBlockMesh::getOrderedBlockCoordinates(int ip, int jp, int kp, double *x_halo, double *y_halo, double *z_halo, double box_p[8][3]){
/*
	int iih_0_0_0;
	int iih_0_0_1;
	int iih_0_1_0;
	int iih_0_1_1;
	int iih_1_0_0;
	int iih_1_0_1;
	int iih_1_1_0;
	int iih_1_1_1;

	//What if we're trying to access *p+1 and we're not periodic in that direction? what happens now in x and z? What should 
	//we return for y? just the origin coordinate?
	//No its fine because C2Decomp will allocate the array with padding even if its not periodic and the *p+1 point should never
	//be returned since we've zero'd the bounding box from above

	int maxIndex = (pySize[0]+2)*(pySize[1]+2)*(pySize[2]+2);


        //This is the halo array index for the same point
        iih_0_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 1;

        //Halo array index for i, j, k+1
        iih_0_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp)*(pySize[0]+2) + ip + 1;

        //Halo array index for i, j+1, k
	if(periodicBCY && jp == (Ny-1)){
            iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
	}else{
            iih_0_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
	}

        //Halo array index for i, j+1, k+1
	if(periodicBCY && jp == (Ny-1)){
            iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 1;
	}else{
            iih_0_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 1;
	}

        //Halo array index for i+1, j, k
        iih_1_0_0 = (kp+1)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

        //Halo array index for i+1, j, k+1
        iih_1_0_1 = (kp+2)*pySize[1]*(pySize[0]+2) + jp*(pySize[0]+2) + ip + 2;

        //Halo array index for i+1, j+1, k
	if(periodicBCY && jp == (Ny-1)){
            iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
	}else{
            iih_1_1_0 = (kp+1)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
	}

        //Halo array index for i+1, j+1, k+1
	if(periodicBCY && jp == (Ny-1)){
	    iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (0)*(pySize[0]+2) + ip + 2;
	}else{
	    iih_1_1_1 = (kp+2)*pySize[1]*(pySize[0]+2) + (jp+1)*(pySize[0]+2) + ip + 2;
	}
 
        bool xEndFlag = false;
        if(pyStart[0] + ip == Nx-1){
            xEndFlag = true;
        }

        bool yEndFlag = false;
        if(jp == Ny-1){
            yEndFlag = true;
        }

        bool zEndFlag = false;
        if(pyStart[2] + kp == Nz-1){
            zEndFlag = true;
        }      

	/////////////////////////////////
	//Definite local point first...//
	//Do 0 0 0 first/////////////////
	//////////////////
	box_p[0][0] = x_halo[iih_0_0_0];
	box_p[0][1] = y_halo[iih_0_0_0];
	box_p[0][2] = z_halo[iih_0_0_0];


	/////////////////
	//Next do 0 0 1// 
	/////////////////

	if(zEndFlag && periodicBCZ){
	    box_p[1][0] = x_halo[iih_0_0_1] + periodicZTranslation[0];
	    box_p[1][1] = y_halo[iih_0_0_1] + periodicZTranslation[1];
	    box_p[1][2] = z_halo[iih_0_0_1] + periodicZTranslation[2];
	}else{ //If here, we're an interior point
	    box_p[1][0] = x_halo[iih_0_0_1];
	    box_p[1][1] = y_halo[iih_0_0_1];
	    box_p[1][2] = z_halo[iih_0_0_1];
	}



	/////////////////
	//Next do 0 1 0//
	/////////////////

	if(yEndFlag && periodicBCY){
	    box_p[2][0] = x_halo[iih_0_1_0] + periodicYTranslation[0];
	    box_p[2][1] = y_halo[iih_0_1_0] + periodicYTranslation[1];
	    box_p[2][2] = z_halo[iih_0_1_0] + periodicYTranslation[2];
	}else{// If, here we're an interior point
	    box_p[2][0] = x_halo[iih_0_1_0];
	    box_p[2][1] = y_halo[iih_0_1_0];
	    box_p[2][2] = z_halo[iih_0_1_0];
	}

	/////////////////
	//Next do 0 1 1//
	/////////////////

	if(zEndFlag && periodicBCZ){
	    if(yEndFlag && periodicBCY){
	        box_p[3][0] = x_halo[iih_0_1_1] + periodicZTranslation[0] + periodicYTranslation[0];
	        box_p[3][1] = y_halo[iih_0_1_1] + periodicZTranslation[1] + periodicYTranslation[1];
	        box_p[3][2] = z_halo[iih_0_1_1] + periodicZTranslation[2] + periodicYTranslation[2];
	    }else{
		box_p[3][0] = x_halo[iih_0_1_1] + periodicZTranslation[0];
		box_p[3][1] = y_halo[iih_0_1_1] + periodicZTranslation[1];
		box_p[3][2] = z_halo[iih_0_1_1] + periodicZTranslation[2];
	    }
	}else{ // in interior domain in z-direction
	    if(yEndFlag && periodicBCY){
		box_p[3][0] = x_halo[iih_0_1_1] + periodicYTranslation[0];
		box_p[3][1] = y_halo[iih_0_1_1] + periodicYTranslation[1];
		box_p[3][2] = z_halo[iih_0_1_1] + periodicYTranslation[2];
	    }else{//If we're here, we're interior 
		box_p[3][0] = x_halo[iih_0_1_1];
		box_p[3][1] = y_halo[iih_0_1_1];
		box_p[3][2] = z_halo[iih_0_1_1];
	    }
	}


	/////////////////
	//Next do 1 0 0//
	/////////////////

	if(xEndFlag && periodicBCX){
	    box_p[4][0] = x_halo[iih_1_0_0] + periodicXTranslation[0];
	    box_p[4][1] = y_halo[iih_1_0_0] + periodicXTranslation[1];
	    box_p[4][2] = z_halo[iih_1_0_0] + periodicXTranslation[2];
	}else{
	    box_p[4][0] = x_halo[iih_1_0_0];
	    box_p[4][1] = y_halo[iih_1_0_0];
	    box_p[4][2] = z_halo[iih_1_0_0];
	}

	/////////////////
	//Next do 1 0 1//
	/////////////////

	if(xEndFlag && periodicBCX){
	    if(zEndFlag && periodicBCZ){
	        box_p[5][0] = x_halo[iih_1_0_1] + periodicXTranslation[0] + periodicZTranslation[0];
	        box_p[5][1] = y_halo[iih_1_0_1] + periodicXTranslation[1] + periodicZTranslation[1];
	        box_p[5][2] = z_halo[iih_1_0_1] + periodicXTranslation[2] + periodicZTranslation[2];
	    }else{
		box_p[5][0] = x_halo[iih_1_0_1] + periodicXTranslation[0];
		box_p[5][1] = y_halo[iih_1_0_1] + periodicXTranslation[1];
		box_p[5][2] = z_halo[iih_1_0_1] + periodicXTranslation[2];
	    }
	}else{
	    if(zEndFlag && periodicBCZ){
		box_p[5][0] = x_halo[iih_1_0_1] + periodicZTranslation[0];
		box_p[5][1] = y_halo[iih_1_0_1] + periodicZTranslation[1];
		box_p[5][2] = z_halo[iih_1_0_1] + periodicZTranslation[2];
	    }else{
		box_p[5][0] = x_halo[iih_1_0_1];
		box_p[5][1] = y_halo[iih_1_0_1];
		box_p[5][2] = z_halo[iih_1_0_1];
	    }
	}

        /////////////////
	//Next do 1 1 0// 
	/////////////////   

	if(xEndFlag && periodicBCX){
	    if(yEndFlag && periodicBCY){
		box_p[6][0] = x_halo[iih_1_1_0] + periodicXTranslation[0] + periodicYTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicXTranslation[1] + periodicYTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicXTranslation[2] + periodicYTranslation[2];
	    }else{
		box_p[6][0] = x_halo[iih_1_1_0] + periodicXTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicXTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicXTranslation[2];
	    }
	}else{
	    if(yEndFlag && periodicBCY){
		box_p[6][0] = x_halo[iih_1_1_0] + periodicYTranslation[0];
		box_p[6][1] = y_halo[iih_1_1_0] + periodicYTranslation[1];
		box_p[6][2] = z_halo[iih_1_1_0] + periodicYTranslation[2];
	    }else{
		box_p[6][0] = x_halo[iih_1_1_0];
		box_p[6][1] = y_halo[iih_1_1_0];
		box_p[6][2] = z_halo[iih_1_1_0];
	    }
	}


        ////////////////////
        //Finally do 1 1 1//
        ////////////////////

        if(xEndFlag && periodicBCX){
	    if(yEndFlag && periodicBCY){
	        if(zEndFlag && periodicBCZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicYTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicYTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicYTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicYTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicYTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicYTranslation[2];
	        }
	    }else{
	        if(zEndFlag && periodicBCZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicXTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicXTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicXTranslation[2];
	        }
	    }
        }else{
	    if(yEndFlag && periodicBCY){
	        if(zEndFlag && periodicBCZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicYTranslation[0] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicYTranslation[1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicYTranslation[2] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicYTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicYTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicYTranslation[2];
	        }
	    }else{
	        if(zEndFlag && periodicBCZ){
		    box_p[7][0] = x_halo[iih_1_1_1] + periodicZTranslation[0];
		    box_p[7][1] = y_halo[iih_1_1_1] + periodicZTranslation[1];
		    box_p[7][2] = z_halo[iih_1_1_1] + periodicZTranslation[2];
	        }else{
		    box_p[7][0] = x_halo[iih_1_1_1];
		    box_p[7][1] = y_halo[iih_1_1_1];
		    box_p[7][2] = z_halo[iih_1_1_1];
	        }
	    }
        }
*/
//};

//void AbstractSingleBlockMesh::generateCoordinateHaloArrays(double *&x_halo, double *&y_halo, double *&z_halo){
/*
	//Really should implement halo transfers for major indexed arrays
        double *x_temp1, *y_temp1, *z_temp1;
        c2d->allocX(x_temp1);
        c2d->allocX(y_temp1);
        c2d->allocX(z_temp1);

        //So move to x-pencil then to non y-major y_pencil...
        c2d->transposeY2X_MajorIndex(x, x_temp1);
        c2d->transposeY2X_MajorIndex(y, y_temp1);
        c2d->transposeY2X_MajorIndex(z, z_temp1);

        //Then back over to y-pencil in x-major array...
        double *x_temp2, *y_temp2, *z_temp2;
        c2d->allocY(x_temp2); c2d->allocY(y_temp2); c2d->allocY(z_temp2);

        c2d->transposeX2Y(x_temp1, x_temp2);
        c2d->transposeX2Y(y_temp1, y_temp2);
        c2d->transposeX2Y(z_temp1, z_temp2);

        delete[] x_temp1;
        delete[] y_temp1;
        delete[] z_temp1;

        c2d->updateHalo(x_temp2, x_halo, 1, 1);
        c2d->updateHalo(y_temp2, y_halo, 1, 1);
        c2d->updateHalo(z_temp2, z_halo, 1, 1);

        delete[] x_temp2;
        delete[] y_temp2;
        delete[] z_temp2;
*/
//};

//int AbstractSingleBlockMesh::findCVForPoint(double p[3], double *x_halo, double *y_halo, double *z_halo){

   /* 
    int base_index = -1;
    FOR_XYZ_YPEN{
	
	base_index = ip;

	//Back out the ip, jp, kp coordinates based off single major index int
        int jpp =  base_index%pySize[1];
        int kpp = (base_index/pySize[1])%pySize[2];
        int ipp =  base_index/(pySize[2]*pySize[1]);

        double box_p[8][3];

        getOrderedBlockCoordinates(ipp, jpp, kpp, x_halo, y_halo, z_halo, box_p);

	if(isPointInHexa(p, box_p)){
	    break;
	}else{
	    base_index = -1;
	}

    }
*/
/*
    int cvListSize, cvList[8192];

    adt->buildListForPoint(cvListSize, cvList, p);

    //This function will return -1 for a point not found and the local control volume index in y-major
    //indexing if it is found
    int base_index = -1;
    for(int ii = 0; ii < cvListSize; ii++){

	base_index = cvList[ii];

	//Back out the ip, jp, kp coordinates based off single major index int
	int jp =  base_index%pySize[1];
	int kp = (base_index/pySize[1])%pySize[2];
	int ip =  base_index/(pySize[2]*pySize[1]);

	double box_p[8][3];

	getOrderedBlockCoordinates(ip, jp, kp, x_halo, y_halo, z_halo, box_p);	

	if(isPointInHexa(p, box_p)){
	    base_index = cvList[ii];
	    break;
	}else{
	    base_index = -1;
	}
    }

    return base_index;    
*/
//}

//void AbstractSingleBlockMesh::initMeshADT(){
/*
	    //Initialize the ADT object for interpolation and locating points in the grid...
	     cout << " > Initializing ADT..." << endl;

            //Get our coordinates in neighbors across partitions using halo updates
            double *x_halo = NULL;
            double *y_halo = NULL;
            double *z_halo = NULL;

	    generateCoordinateHaloArrays(x_halo, y_halo, z_halo);
  
	    int Nlocal = pySize[0]*pySize[1]*pySize[2];
            double (*boundBoxMin)[3] = new double[Nlocal][3];
            double (*boundBoxMax)[3] = new double[Nlocal][3];
	   
	    
 
            //Cycle through the halo arrays of coordinates
            for(int kp = 0; kp < pySize[2]; kp++){
                for(int jp = 0; jp < pySize[1]; jp++){
                    for(int ip = 0; ip < pySize[0]; ip++){

                        //This is the non-halo array index
                        int ii = kp*pySize[1]*pySize[0] + jp*pySize[0] + ip;

                        //This is the non-halo array index, y-major
                        int ii_major = ip*pySize[2]*pySize[1] + kp*pySize[1] + jp;


		        bool xEndFlag = false;
			if(pyStart[0] + ip == Nx-1){
			    xEndFlag = true;
			}

			bool yEndFlag = false;
			if(jp == Ny-1){
			    yEndFlag = true;
			}

			bool zEndFlag = false;
			if(pyStart[2] + kp == Nz-1){
			    zEndFlag = true;
			}      

			bool noBBFlag = false;	
                        if((xEndFlag && !periodicBCX) ||
			   (yEndFlag && !periodicBCY) ||
			   (zEndFlag && !periodicBCZ)){
                            noBBFlag = true;
                        }


			double box_p[8][3];
			if(!noBBFlag){
			    getOrderedBlockCoordinates(ip, jp, kp, x_halo, y_halo, z_halo, box_p);
			}

			double x_max = -1.0e100; 
			double y_max = -1.0e100; 
			double z_max = -1.0e100; 
			double x_min =  1.0e100;
			double y_min =  1.0e100;
			double z_min =  1.0e100;

			if(!noBBFlag){

			    for(int iip = 0; iip < 8; iip++){
    			       x_max = fmax(x_max, box_p[iip][0]);
    			       x_min = fmin(x_min, box_p[iip][0]);

    			       y_max = fmax(y_max, box_p[iip][1]);
    			       y_min = fmin(y_min, box_p[iip][1]);

    			       z_max = fmax(z_max, box_p[iip][2]);
    			       z_min = fmin(z_min, box_p[iip][2]);
			    }

		  	}else{
	 		    x_max = x[ii];
			    x_min = x[ii];
			    y_max = y[ii];
			    y_min = y[ii];
			    z_min = z[ii];
			    z_min = z[ii];
			}

			//We'll usually be accessing this in the major indexing fashion
             	        boundBoxMin[ii_major][0] = x_min;
             	        boundBoxMin[ii_major][1] = y_min;
             	        boundBoxMin[ii_major][2] = z_min;
			
;
		        boundBoxMax[ii_major][0] = x_max; 
		        boundBoxMax[ii_major][1] = y_max; 
		        boundBoxMax[ii_major][2] = z_max; 
                    }
                }
            }

	    FOR_XYZ_YPEN{
		FOR_I3{
		    double delta = 1.0E-6*(boundBoxMax[ip][i] - boundBoxMin[ip][i]);
		    boundBoxMax[ip][i] += delta;
		    boundBoxMin[ip][i] -= delta;
		}
	    }

	     cout << " > Done getting bounding boxes for the CV's, initializing ADT... " << endl;

	    adt = new Adt<double>(Nlocal, boundBoxMin, boundBoxMax);

	     cout << " > Done!" << endl;


	    //Dump the solution grid
	    //dumpGrid();

	    delete[] boundBoxMin;
	    delete[] boundBoxMax;
	    delete[] x_halo;
	    delete[] y_halo;
	    delete[] z_halo;
*/

//}


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
