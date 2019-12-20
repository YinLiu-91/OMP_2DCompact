#include "LADKawai.hpp"

void LADKawai::calcVelocityTensorStuff(double *gradU[2][2]){

    #pragma omp parallel for
    FOR_XY{
	double S00, S01, S11;
	S00 = gradU[0][0][ip];
	S11 = gradU[1][1][ip];
	S01 = 0.5*(gradU[0][1][ip] + gradU[1][0][ip]);
	S[ip] = sqrt(2*S00*S00 + 2*S11*S11 + 4*S01);

	dil[ip] = S00 + S11;

	vort[ip] = gradU[1][0][ip] - gradU[0][1][ip];
    }
}

void LADKawai::calc4thOrderDerivative(double *phi, double *d4phi0, double *d4phi1, double *work1, double *work2){

    derivX->calc1stDerivField(phi,   work1);
    derivX->calc1stDerivField(work1, work2);
    derivX->calc1stDerivField(work2, work1);
    derivX->calc1stDerivField(work1, d4phi0);

    transposeMatrix_Fast2(phi, Nx, Ny, d4phi1, cs->opt->blocksize);
    derivY->calc1stDerivField(d4phi1, work1);
    derivY->calc1stDerivField(work1,  work2);
    derivY->calc1stDerivField(work2,  work1);
    derivY->calc1stDerivField(work1,  work2);
    transposeMatrix_Fast2(work2, Ny, Nx, d4phi1, cs->opt->blocksize);
 
}

void LADKawai::calcLADViscosity(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE){

    //compute the spacing component of the function
    double (*deltaComp)[2] = new double[Nx*Ny][2];
    
    #pragma omp parallel for collapse(2)
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

    calc4thOrderDerivative(S, dFmu4dx04, dFmu4dx14, temp1, temp2);

    getRange(S, "S", Nx, Ny);
    getRange(dFmu4dx04, "dFmu4dx04", Nx, Ny);
    getRange(dFmu4dx14, "dFmu4dx14", Nx, Ny);

    //Calculating the thing thats going to get filtered...
    
    #pragma omp parallel for
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


    //Calculate the unit density gradient...
    double *vi[] = {rho};
    vector<double *> vecIn(vi, vi+sizeof(vi)/sizeof(vi[0]));

    double *drhod0 = new double[Nx*Ny]; 
    double *drhod1 = new double[Nx*Ny]; 
    double *vo[] = {drhod0, drhod1};
    vector<double *> vecOut(vo, vo+sizeof(vo)/sizeof(vo[0]));

    cs->computeGradient(vecIn, vecOut); 

    double *drhodx = new double[Nx*Ny];
    double *drhody = new double[Nx*Ny];

    #pragma omp parallel for
    FOR_XY{
	drhodx[ip] = cs->J11[ip]*drhod0[ip] + cs->J21[ip]*drhod1[ip];
	drhody[ip] = cs->J12[ip]*drhod0[ip] + cs->J22[ip]*drhod1[ip];
	
	double drhomag = sqrt(drhodx[ip]*drhodx[ip] + drhody[ip]*drhody[ip]);
	drhodx[ip] /= drhomag;
	drhody[ip] /= drhomag;
    } 

    delete[] drhod0;
    delete[] drhod1;

    //computing the dilatation and shock sensor...
    double *temp1	= new double[Nx*Ny];
    double *temp2	= new double[Nx*Ny];
    double *dFbeta4dx04 = new double[Nx*Ny];
    double *dFbeta4dx14 = new double[Nx*Ny];
    calc4thOrderDerivative(dil, dFbeta4dx04, dFbeta4dx14, temp1, temp2);

    #pragma omp parallel for
    FOR_XY{
	double eps = 1E-16;
	double ducros = (dil[ip]*dil[ip])/(dil[ip]*dil[ip] + vort[ip]*vort[ip] + eps);	
	
	double heavyside;
	if(-dil[ip] > 0){
	    heavyside = 1.0;
	}else{
	    heavyside = 0.0;
	}
	
	fsw[ip] = ducros*heavyside;	
    }
    
    //compute the spacing component of the function
    double (*deltaComp)[2] = new double[Nx*Ny][2];
 
    #pragma omp parallel for collapse(2)
    FOR_Y{
	FOR_X{
	    int ip = GET2DINDEX_XY;
	
	    double delta_xl[2][2]; 
	    double delta2_l_beta[2];

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

	    delta2_l_beta[0] = delta_xl[0][0]*drhodx[ip] + delta_xl[0][1]*drhody[ip]; 
	    delta2_l_beta[1] = delta_xl[1][0]*drhodx[ip] + delta_xl[1][1]*drhody[ip]; 

	    //Square it...
	    delta2_l_beta[0] *= delta2_l_beta[0];
	    delta2_l_beta[1] *= delta2_l_beta[1];
	    

	    //Multiply by the generalized coordinate spacing...
	    //Assumine its constant, which is like, the point
	    double dxi[2];
	    dxi[0] = cs->dom->x[1]-cs->dom->x[0]; 
	    dxi[1] = cs->dom->y[1]-cs->dom->y[0]; 

	    //Doing r=4
	    deltaComp[ip][0] = pow(dxi[0],4.0)*delta2_l_beta[0];
	    deltaComp[ip][1] = pow(dxi[1],4.0)*delta2_l_beta[1];
	}
    }

    //calculating the whole thing thats getting filtered...
    #pragma omp parallel for
    FOR_XY{
	temp1[ip] = rho[ip]*fsw[ip]*(dFbeta4dx04[ip]*deltaComp[ip][0] + dFbeta4dx14[ip]*deltaComp[ip][1]); //may need to be multiplied by wall damping function
    }

    //Filter in both direcitons
    filtX->filterField(temp1, temp2);
    
    transposeMatrix_Fast2(temp2, Nx, Ny, temp1, cs->opt->blocksize);
    filtY->filterField(temp1, temp2);
    transposeMatrix_Fast2(temp2, Ny, Nx, temp1, cs->opt->blocksize);

    #pragma omp parallel for
    FOR_XY{
	beta_star[ip] = C_beta*temp1[ip];
    }



    delete[] temp1;
    delete[] temp2;
    delete[] dFbeta4dx04;
    delete[] dFbeta4dx14;
    delete[] drhodx;
    delete[] drhody;
    delete[] deltaComp;

};


void LADKawai::calcLADK(double *gradU[2][2], double *rho, double *rhoU, double *rhoV, double *rhoE){

    //Calculate e, could probably just pull p from solver but with the inputs we'll just do it here...
    double *e = new double[Nx*Ny];
    #pragma omp parallel for
    FOR_XY{
	double U = rhoU[ip]/rho[ip];
	double V = rhoV[ip]/rho[ip];
	double p = cs->ig->solvep(rho[ip], rhoE[ip], U, V);	
	e[ip] = p/((cs->ig->gamma - 1.0)*rho[ip]);
    }

    double *vi[] = {e};
    vector<double *> vecIn(vi, vi+sizeof(vi)/sizeof(vi[0]));

    double *ded0 = new double[Nx*Ny]; 
    double *ded1 = new double[Nx*Ny]; 
    double *vo[] = {ded0, ded1};
    vector<double *> vecOut(vo, vo+sizeof(vo)/sizeof(vo[0]));

    cs->computeGradient(vecIn, vecOut); 

    double *dedx = new double[Nx*Ny];
    double *dedy = new double[Nx*Ny];
    #pragma omp parallel for
    FOR_XY{
	dedx[ip] = cs->J11[ip]*ded0[ip] + cs->J21[ip]*ded1[ip];
	dedy[ip] = cs->J12[ip]*ded0[ip] + cs->J22[ip]*ded1[ip];
	
	double demag = sqrt(dedx[ip]*dedx[ip] + dedy[ip]*dedy[ip]);
	dedx[ip] /= demag;
	dedy[ip] /= demag;
    } 

    delete[] ded0;
    delete[] ded1;

    //computing the internal energy 4th derivative...
    double *temp1    = new double[Nx*Ny];
    double *temp2    = new double[Nx*Ny];
    double *dFk4dx04 = new double[Nx*Ny];
    double *dFk4dx14 = new double[Nx*Ny];
    calc4thOrderDerivative(e, dFk4dx04, dFk4dx14, temp1, temp2);

    //compute the spacing component of the function
    double (*deltaComp)[2] = new double[Nx*Ny][2];
 
    #pragma omp parallel for collapse(2)
    FOR_Y{
	FOR_X{
	    int ip = GET2DINDEX_XY;
	
	    double delta_xl[2][2]; 
	    double delta_l_k[2];

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

	    delta_l_k[0] = fabs(delta_xl[0][0]*dedx[ip] + delta_xl[0][1]*dedy[ip]); 
	    delta_l_k[1] = fabs(delta_xl[1][0]*dedx[ip] + delta_xl[1][1]*dedy[ip]); 


	    //Multiply by the generalized coordinate spacing...
	    //Assumine its constant, which is like, the point
	    double dxi[2];
	    dxi[0] = cs->dom->x[1]-cs->dom->x[0]; 
	    dxi[1] = cs->dom->y[1]-cs->dom->y[0]; 

	    //Doing r=4
	    deltaComp[ip][0] = pow(dxi[0],4.0)*delta_l_k[0];
	    deltaComp[ip][1] = pow(dxi[1],4.0)*delta_l_k[1];
	}
    }

    //calculating the thing thats going to get filtered...
    #pragma omp parallel for
    FOR_XY{
	double U = rhoU[ip]/rho[ip];
	double V = rhoV[ip]/rho[ip];
	double p = cs->ig->solvep(rho[ip], rhoE[ip], U, V);	

	double c = cs->ig->solveSOS(rho[ip], p);
	double T = cs->ig->solveT(rho[ip], p);

	temp1[ip] = (rho[ip]*c/T)*(dFk4dx04[ip]*deltaComp[ip][0] + dFk4dx14[ip]*deltaComp[ip][1]);
	  
    }

    //Filter in both direcitons
    filtX->filterField(temp1, temp2);
    
    transposeMatrix_Fast2(temp2, Nx, Ny, temp1, cs->opt->blocksize);
    filtY->filterField(temp1, temp2);
    transposeMatrix_Fast2(temp2, Ny, Nx, temp1, cs->opt->blocksize);

    #pragma omp parallel for
    FOR_XY{
	k_star[ip] = C_k*temp1[ip];
    }

    delete[] dedx;
    delete[] dedy;
    delete[] e;
    delete[] temp1;
    delete[] temp2;
    delete[] dFk4dx04;
    delete[] dFk4dx14;
    delete[] deltaComp;

    getRange(mu_star,   "LAD->MU_STAR   ", Nx, Ny);
    getRange(beta_star, "LAD->BETA_STAR ", Nx, Ny);
    getRange(k_star,    "LAD->K_STAR    ", Nx, Ny);

};
