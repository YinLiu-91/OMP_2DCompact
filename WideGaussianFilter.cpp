#include "WideGaussianFilter.hpp"

using namespace std;

void WideGaussianFilter::filterPeriodic(double *phi, double *phiF){

    for(int ip = 0; ip < N; ip++){
	if(ip == 0){
	    phiF[ip] = a4*(phi[N-4] + phi[4]) + \
		       a3*(phi[N-3] + phi[3]) + \
		       a2*(phi[N-2] + phi[2]) + \
		       a1*(phi[N-1] + phi[1]) + \
		       a0*(phi[0]);
	}else if(ip == 1){
	    phiF[ip] = a4*(phi[N-3] + phi[5]) + \
		       a3*(phi[N-2] + phi[4]) + \
		       a2*(phi[N-1] + phi[3]) + \
		       a1*(phi[0]   + phi[2]) + \
		       a0*(phi[1]);
	}else if(ip == 2){
	    phiF[ip] = a4*(phi[N-2] + phi[6]) + \
		       a3*(phi[N-1] + phi[5]) + \
		       a2*(phi[0]   + phi[4]) + \
		       a1*(phi[1]   + phi[3]) + \
		       a0*(phi[2]);
	}else if(ip == 3){
	    phiF[ip] = a4*(phi[N-1] + phi[7]) + \
		       a3*(phi[0]   + phi[6]) + \
		       a2*(phi[1]   + phi[5]) + \
		       a1*(phi[2]   + phi[4]) + \
		       a0*(phi[3]);
	}else if(ip == N-4){
	    phiF[ip] = a4*(phi[N-8] + phi[0])   + \
		       a3*(phi[N-7] + phi[N-1]) + \
		       a2*(phi[N-6] + phi[N-2]) + \
		       a1*(phi[N-5] + phi[N-3]) + \
		       a0*(phi[N-4]);
	}else if(ip == N-3){
	    phiF[ip] = a4*(phi[N-7] + phi[1])   + \
		       a3*(phi[N-6] + phi[0])   + \
		       a2*(phi[N-5] + phi[N-1]) + \
		       a1*(phi[N-4] + phi[N-2]) + \
		       a0*(phi[N-3]);
	}else if(ip == N-2){
	    phiF[ip] = a4*(phi[N-6] + phi[2])   + \
		       a3*(phi[N-5] + phi[1])   + \
		       a2*(phi[N-4] + phi[0])   + \
		       a1*(phi[N-3] + phi[N-1]) + \
		       a0*(phi[N-2]);
	}else if(ip == N-1){
	    phiF[ip] = a4*(phi[N-5] + phi[3]) + \
		       a3*(phi[N-4] + phi[2]) + \
		       a2*(phi[N-3] + phi[1]) + \
		       a1*(phi[N-2] + phi[0]) + \
		       a0*(phi[N-1]);
	}else{
	    phiF[ip] = a4*(phi[ip-4] + phi[ip+4]) + \
		       a3*(phi[ip-3] + phi[ip+3]) + \
		       a2*(phi[ip-2] + phi[ip+2]) + \
		       a1*(phi[ip-1] + phi[ip+1]) + \
		       a0*(phi[ip]);
	}	
    }
}

void WideGaussianFilter::filterDirichlet(double *phi, double *phiF){


    //Probably a better small support gaussian filter than can be used here
    //instead of reflecting the filter
    for(int ip = 0; ip < N; ip++){
	if(ip == 0){
	    phiF[ip] = a4*(phi[4] + phi[4]) + \
		       a3*(phi[3] + phi[3]) + \
		       a2*(phi[2] + phi[2]) + \
		       a1*(phi[1] + phi[1]) + \
		       a0*(phi[0]);
	}else if(ip == 1){
	    phiF[ip] = a4*(phi[3] + phi[5]) + \
		       a3*(phi[2] + phi[4]) + \
		       a2*(phi[1] + phi[3]) + \
		       a1*(phi[0] + phi[2]) + \
		       a0*(phi[1]);
	}else if(ip == 2){
	    phiF[ip] = a4*(phi[2] + phi[6]) + \
		       a3*(phi[1] + phi[5]) + \
		       a2*(phi[0] + phi[4]) + \
		       a1*(phi[1] + phi[3]) + \
		       a0*(phi[2]);
	}else if(ip == 3){
	    phiF[ip] = a4*(phi[1] + phi[7]) + \
		       a3*(phi[0] + phi[6]) + \
		       a2*(phi[1] + phi[5]) + \
		       a1*(phi[2] + phi[4]) + \
		       a0*(phi[3]);
	}else if(ip == N-4){
	    phiF[ip] = a4*(phi[N-8] + phi[N-2])   + \
		       a3*(phi[N-7] + phi[N-1]) + \
		       a2*(phi[N-6] + phi[N-2]) + \
		       a1*(phi[N-5] + phi[N-3]) + \
		       a0*(phi[N-4]);
	}else if(ip == N-3){
	    phiF[ip] = a4*(phi[N-7] + phi[N-3])   + \
		       a3*(phi[N-6] + phi[N-2])   + \
		       a2*(phi[N-5] + phi[N-1]) + \
		       a1*(phi[N-4] + phi[N-2]) + \
		       a0*(phi[N-3]);
	}else if(ip == N-2){
	    phiF[ip] = a4*(phi[N-6] + phi[N-4])   + \
		       a3*(phi[N-5] + phi[N-3])   + \
		       a2*(phi[N-4] + phi[N-2])   + \
		       a1*(phi[N-3] + phi[N-1]) + \
		       a0*(phi[N-2]);
	}else if(ip == N-1){
	    phiF[ip] = a4*(phi[N-5] + phi[N-5]) + \
		       a3*(phi[N-4] + phi[N-4]) + \
		       a2*(phi[N-3] + phi[N-3]) + \
		       a1*(phi[N-2] + phi[N-2]) + \
		       a0*(phi[N-1]);
	}else{
	    phiF[ip] = a4*(phi[ip-4] + phi[ip+4]) + \
		       a3*(phi[ip-3] + phi[ip+3]) + \
		       a2*(phi[ip-2] + phi[ip+2]) + \
		       a1*(phi[ip-1] + phi[ip+1]) + \
		       a0*(phi[ip]);
	}	
    }

}

void WideGaussianFilter::filter(double *phi, double *phiF){

    if(bcType == Options::PERIODIC_SOLVE){
	filterPeriodic(phi, phiF);
    }else if(bcType == Options::DIRICHLET_SOLVE){
	filterDirichlet(phi, phiF);
    }

}

void WideGaussianFilter::filterField(double *dataIn, double *dataOut){


    const int numThreads = NUMTHREADSNEST;

    if(currentDir == AbstractDerivatives::DIRX){

	#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    	FOR_Y{

	    bool passFlag = false;
	    if(bc->bcYType == Options::DIRICHLET_SOLVE && (j == 0 || j == Ny-1)){
		passFlag = true;
	    }

    	    double *dataInLocal, *dataOutLocal;
            int ii = j*Nx;
	    dataInLocal  = &dataIn[ii];
            dataOutLocal = &dataOut[ii];
	    if(passFlag){
	        memcpy(dataOutLocal, dataInLocal, sizeof(double)*Nx);
	    }else{
		filter(dataInLocal, dataOutLocal);
	    }

	}

    }else if(currentDir == AbstractDerivatives::DIRY){

	#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
        FOR_X{
	    bool passFlag = false;
	    if(bc->bcXType == Options::DIRICHLET_SOLVE && (i == 0 || i == Nx-1))
		passFlag = true;

	    double *dataInLocal, *dataOutLocal;
            int ii = i*Ny;
            dataInLocal  = &dataIn[ii];
            dataOutLocal = &dataOut[ii];
	    if(passFlag){
		memcpy(dataOutLocal, dataInLocal, sizeof(double)*Ny);
	    }else{
                filter(dataInLocal, dataOutLocal);
	    }
        }

    }

}
 
