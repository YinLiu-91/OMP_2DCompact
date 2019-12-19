#include "LADKawai.hpp"

void LADKawai::calcVelocityTensorStuff(double *gradU[2][2]){
    FOR_XY{
	double S00, S01, S11;
	S00 = gradU[0][0][ip];
	S11 = gradU[1][1][ip];
	S01 = 0.5*(gradU[0][1][ip] + gradU[1][0][ip]);
	S[ip] = sqrt(2*S00*S00 + 2*S11*S11 + 4*S01);

	dil[ip] = S00 + S11;
    }
}

void LADKawai::calcLADViscosity(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE){

    //compute the spacing component of the function
    double (*deltaComp)[2] = new double[Nx*Ny][2];
    
    FOR_Y{
	FOR_X{
	    int ip = GET2DINDEX_XY;
	
	    double delta_xl[2][2]; 
	    double delta2_l_mu[2];

	    if(i == 0){
	        delta_xl[0][0] = cs->msh->x[i+1]-cs->msh->x[i];
	        delta_xl[0][1] = cs->msh->y[i+1]-cs->msh->y[i];
	    }else if(Nx-1){
		delta_xl[0][0] = cs->msh->x[Nx-1]-cs->msh->x[Nx-2];
		delta_xl[0][1] = cs->msh->y[Nx-1]-cs->msh->y[Nx-2];
	    }else{
		delta_xl[0][0] = 0.5*(cs->msh->x[i+1]-cs->msh->x[i-1]);
		delta_xl[0][1] = 0.5*(cs->msh->y[i+1]-cs->msh->y[i-1]);
	    }

	    if(j == 0){
	        delta_xl[1][0] = cs->msh->x[j+1]-cs->msh->x[j];
	        delta_xl[1][1] = cs->msh->y[j+1]-cs->msh->y[j];
	    }else if(Ny-1){
		delta_xl[1][0] = cs->msh->x[Ny-1]-cs->msh->x[Ny-2];
		delta_xl[1][1] = cs->msh->y[Ny-1]-cs->msh->y[Ny-2];
	    }else{
		delta_xl[1][0] = 0.5*(cs->msh->x[j+1]-cs->msh->x[j-1]);
		delta_xl[1][1] = 0.5*(cs->msh->y[j+1]-cs->msh->y[j-1]);
	    }

	    delta2_l_mu[0] = (delta_xl[0][0]*delta_xl[0][0] + delta_xl[0][1]*delta_xl[0][1]); 
	    delta2_l_mu[1] = (delta_xl[1][0]*delta_xl[1][0] + delta_xl[1][1]*delta_xl[1][1]); 

	    //Multiply by the generalized coordinate spacing...
	    //Assumine its constant, which is like, the point
	    double dxi[2];
	    dxi[0] = cs->dom->x[1]-cs->dom->x[0]; 
	    dxi[1] = cs->dom->y[1]-cs->dom->y[0]; 

	    //Doing r=4
	    deltaComp[ip][0] = pow(dxi[0],4.0)*delta2_l_mu[0];
	    deltaComp[ip][1] = pow(dxi[1],4.0)*delta2_l_mu[1];
	}
    }

    //Doing the derivative four times for now, need to implement just explicit...
    double *temp1 = new double[Nx*Ny];  
    double *temp2 = new double[Nx*Ny];  
    double *dFmu4dx04 = new double[Nx*Ny];
    double *dFmu4dx14 = new double[Nx*Ny];

    derivX->calc1stDerivField(S,     temp1);
    derivX->calc1stDerivField(temp1, temp2);
    derivX->calc1stDerivField(temp2, temp1);
    derivX->calc1stDerivField(temp1, dFmu4dx04);

    transposeMatrix_Fast2(S, Nx, Ny, dFmu4dx14, cs->opt->blocksize);
    derivY->calc1stDerivField(dFmu4dx14, temp1);
    derivY->calc1stDerivField(temp1,     temp2);
    derivY->calc1stDerivField(temp2,     temp1);
    derivY->calc1stDerivField(temp1,     temp2);
    transposeMatrix_Fast2(temp2, Ny, Nx, dFmu4dx14, cs->opt->blocksize);
 
    //Calculating the thing thats going to get filtered...
    FOR_XY{
	temp1[ip] = rho[ip]*fabs(dFmu4dx04[ip]*deltaComp[ip][0] + dFmu4dx14[ip]*deltaComp[ip][1]);//Should be multiplied by the wall damping function...
    }

    //Filter in both direcitons
    filtX->filterField(temp1, temp2);
    
    transposeMatrix_Fast2(temp2, Nx, Ny, temp1, cs->opt->blocksize);
    filtY->filterField(temp1, temp2);
    transposeMatrix_Fast2(temp2, Ny, Nx, temp1, cs->opt->blocksize);

    FOR_XY{
	mu_star[ip] = C_mu*temp1[ip];
    }

    delete[] temp1;
    delete[] temp2;
    delete[] dFmu4dx04;
    delete[] dFmu4dx14;
    delete[] (deltaComp);

}


void LADKawai::calcLADBeta(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE){

};


void LADKawai::calcLADK(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE){

};
