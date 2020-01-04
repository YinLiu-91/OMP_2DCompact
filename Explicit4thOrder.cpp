#include "Explicit4thOrder.hpp"

void Explicit4thOrder::calc1stDerivField(double *dataIn, double *dataOut){

    const int numThreads = NUMTHREADSNEST;

    if(currentDir == DIRX){

	#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
	FOR_Y{
 	    double *dataInLocal, *dataOutLocal;
	    int ii = j*Nx;
	    dataInLocal  = &dataIn[ii];
	    dataOutLocal = &dataOut[ii];
	    calc1stDeriv(dataInLocal, dataOutLocal);
	}


    }else if(currentDir == DIRY){
	#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
	FOR_X{
	    double *dataInLocal, *dataOutLocal;
	    int ii = i*Ny;
	    dataInLocal  = &dataIn[ii];
	    dataOutLocal = &dataOut[ii];
	    calc1stDeriv(dataInLocal, dataOutLocal);
	}
    }

}

void Explicit4thOrder::calc1stDeriv(double *phi, double *dphi){

    if(bcType == Options::PERIODIC_SOLVE){
	Calc4thPeriodic(phi, dphi);	
    }else if(bcType == Options::DIRICHLET_SOLVE){
	Calc4thDirichlet(phi, dphi);	
    }

}

void Explicit4thOrder::Calc4thPeriodic(double *phi, double *dphidx){



    dphidx[0]  = (b/6.0)*(phi[3] - 9.0*phi[1] + 16.0*phi[0] - 9.0*phi[N-1] + phi[N-3]);
    dphidx[0] +=       a*(phi[2] - 4.0*phi[1] +  6.0*phi[0] - 4.0*phi[N-1] + phi[N-2]);                                              

    dphidx[1]  = (b/6.0)*(phi[4] - 9.0*phi[2] + 16.0*phi[1] - 9.0*phi[0]   + phi[N-2]);
    dphidx[1] +=       a*(phi[3] - 4.0*phi[2] +  6.0*phi[1] - 4.0*phi[0]   + phi[N-1]);                                              

    dphidx[2]  = (b/6.0)*(phi[5] - 9.0*phi[3] + 16.0*phi[2] - 9.0*phi[1]   + phi[N-1]);
    dphidx[2] +=       a*(phi[4] - 4.0*phi[3] +  6.0*phi[2] - 4.0*phi[1]   + phi[0]);                                              

    for(int ip = 3; ip < N-3; ip++){                                           
        dphidx[ip]  = (b/6.0)*(phi[ip+3] - 9.0*phi[ip+1] + 16.0*phi[ip] - 9.0*phi[ip-1] + phi[ip-3]);
        dphidx[ip] +=       a*(phi[ip+2] - 4.0*phi[ip+1] +  6.0*phi[ip] - 4.0*phi[ip-1] + phi[ip-2]);                                              
    }   

    dphidx[N-1]  = (b/6.0)*(phi[2] - 9.0*phi[0] + 16.0*phi[N-1] - 9.0*phi[N-2]   + phi[N-4]);
    dphidx[N-1] +=       a*(phi[1] - 4.0*phi[0] +  6.0*phi[N-1] - 4.0*phi[N-2]   + phi[N-3]);                                              

    dphidx[N-2]  = (b/6.0)*(phi[1] - 9.0*phi[N-1] + 16.0*phi[N-2] - 9.0*phi[N-3]   + phi[N-5]);
    dphidx[N-2] +=       a*(phi[0] - 4.0*phi[N-1] +  6.0*phi[N-2] - 4.0*phi[N-3]   + phi[N-4]);                                              

    dphidx[N-3]  = (b/6.0)*(phi[0]   - 9.0*phi[N-2] + 16.0*phi[N-3] - 9.0*phi[N-4]   + phi[N-6]);
    dphidx[N-3] +=       a*(phi[N-1] - 4.0*phi[N-2] +  6.0*phi[N-3] - 4.0*phi[N-4]   + phi[N-5]);                                              

    for(int ip = 0; ip < N; ip++){
        dphidx[ip] /= pow(dd,4.0);
    }


}

void Explicit4thOrder::Calc4thDirichlet(double *phi, double *dphidx){


    dphidx[0] = a11*phi[0]   + a21*phi[1] + a31*phi[2] +             
                     a41*phi[3] + a51*phi[4] + a61*phi[5];              
    dphidx[1] = a12*phi[0] + a22*phi[1]   + a32*phi[2] +             
                     a42*phi[3] + a52*phi[4] + a62*phi[5];              
    dphidx[2] = a13*phi[0] + a23*phi[1] + a33*phi[2] +               
                     a43*phi[3] + a53*phi[4];  

    for(int ip = 3; ip < N-4; ip++){                                           
        dphidx[ip]  = (b/6.0)*(phi[ip+3] - 9.0*phi[ip+1] + 16.0*phi[ip] - 9.0*phi[ip-1] + phi[ip-3]);
        dphidx[ip] +=       a*(phi[ip+2] - 4.0*phi[ip+1] +  6.0*phi[ip] - 4.0*phi[ip-1] + phi[ip-2]);                                              
    }   

    dphidx[N-1] = a11*phi[N-1]   + a21*phi[N-2] + a31*phi[N-3] +             
                     a41*phi[N-4] + a51*phi[N-5] + a61*phi[N-6];
    dphidx[N-2] = a12*phi[N-1] + a22*phi[N-2]   + a32*phi[N-3] +             
                     a42*phi[N-4] + a52*phi[N-5] + a62*phi[N-6];
    dphidx[N-3] = a13*phi[N-1] + a23*phi[N-2] + a33*phi[N-3] +
                     a43*phi[N-4] + a53*phi[N-5];

    for(int ip = 0; ip < N; ip++){
        dphidx[ip] /= pow(dd,4.0);
    }

}


