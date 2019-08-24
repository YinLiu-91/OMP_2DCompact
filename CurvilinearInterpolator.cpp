#include "CurvilinearInterpolator.hpp"

void CurvilinearInterpolator::interpolateData(double *dataIn, double *interpedDataOut){

    //Get the halo of the data
    for(int ii = 0; ii < pointFoundCount; ii++){

        //Get the ip, jp, kp based off index
        int icv = icvList[ii];

        int i = icv%Nx;
	int j = (icv-i)/Nx;

	//Pull the required data from the halo array
	double box_data[4];

	box_data[0] = dataIn[icv];

	int iii;
	if(i == Nx-1 && cs->msh->periodicBCX){
	    iii = j*Nx + 0;
	}else{
	    iii = j*Nx + i + 1;
	}
	box_data[1] = dataIn[iii];

	if(j == Ny-1 && cs->msh->periodicBCY){
	    iii =   (0)*Nx + i;
	}else{
    	    iii = (j+1)*Nx + i;
	}
	box_data[2] = dataIn[iii];

	if((j==Ny-1) && (i==Nx-1) && cs->msh->periodicBCY && cs->msh->periodicBCX){
	    iii = 0;
	}else if(i==Nx-1 && cs->msh->periodicBCX){
	    iii = (j+1)*Nx + 0;
	}else if(j==Ny-1 && cs->msh->periodicBCY){
	    iii =  (0)*Nx + i + 1;
	}else{
	    iii = (j+1)*Nx + i + 1;
	}
	box_data[3] = dataIn[iii];

	//Interpolate using the weights for this point...
	interpedDataOut[ii] = 0;
	for(int jj = 0; jj < 4; jj++){
	    interpedDataOut[ii] += Ni[ii][jj]*box_data[jj]; 
	}

	//For now just give us something
	//interpedDataOut[ii] = dataIn[icv];
//	cout << icv << " " << interpedDataOut[ii] << " " << cs->msh->x[icv] << endl;
    }
  
};
