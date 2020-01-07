#include "CurvilinearCSolver.hpp"

void CurvilinearCSolver::initializeSolverData(){

   
//    if(useTiming) ft1 = MPI_Wtime();

        cout << endl;
        cout << " > Allocating Solver Arrays..." << endl;
//        double workSize = 0;
//        workSize = 104.0 * (double)N * 8.0;
//        cout << " > Need " << workSize/1024.0/1024.0/1024.0 << " Gb of memory required to allocate solver arrays " << endl;


    //3
    dU1 = new double[N];
    dU2 = new double[N];


    //6
    dV1 = new double[N];
    dV2 = new double[N];

    //9
    dW1 = new double[N];
    dW2 = new double[N];

    //12
    dT1 = new double[N];
    dT2 = new double[N];

    //18
    Tau11 = new double[N];
    Tau12 = new double[N];
    Tau22 = new double[N];

    Tau21 = Tau12;
   
    //21
    cont_1 = new double[N];
    cont_2 = new double[N];

    //24
    mom1_1 = new double[N];
    mom1_2 = new double[N];

    //27
    mom2_1 = new double[N];
    mom2_2 = new double[N];

    //33
    engy_1 = new double[N];
    engy_2 = new double[N];

    //37
    rho1 = new double[N];
    rhok = new double[N];
    rhok2 = new double[N];
    rho2 = new double[N];

    //41
    rhoU1 = new double[N];
    rhoUk = new double[N];
    rhoUk2 = new double[N];
    rhoU2 = new double[N];

    //45
    rhoV1 = new double[N];
    rhoVk = new double[N];
    rhoVk2 = new double[N];
    rhoV2 = new double[N];

    //53
    rhoE1 = new double[N];
    rhoEk = new double[N];
    rhoEk2 = new double[N];
    rhoE2 = new double[N];

    //58 these will be cleared though...
    rho0 = new double[N];
    U0 = new double[N];
    V0 = new double[N];
    p0 = new double[N];

    //61
    U = new double[N];
    V = new double[N];
    T = new double[N];
    Ucurv = new double[N];
    Vcurv = new double[N];

    //72
    p = new double[N];
    mu = new double[N];
    sos = new double[N];

    mu_eff = new double[N];
    k_eff  = new double[N];

    //84

    temp1 = new double[N]; tempVec.push_back(temp1);
    temp2 = new double[N]; tempVec.push_back(temp2);
    temp3 = new double[N]; tempVec.push_back(temp3);
    temp4 = new double[N]; tempVec.push_back(temp4);
    temp5 = new double[N]; tempVec.push_back(temp5);
    temp6 = new double[N]; tempVec.push_back(temp6);
    temp7 = new double[N]; tempVec.push_back(temp7);
    temp8 = new double[N]; tempVec.push_back(temp8);
    temp9 = new double[N]; tempVec.push_back(temp9);
    temp10= new double[N]; tempVec.push_back(temp10);
    temp11= new double[N]; tempVec.push_back(temp11);
    temp12= new double[N]; tempVec.push_back(temp12);
    temp13= new double[N]; tempVec.push_back(temp13);
    temp14= new double[N]; tempVec.push_back(temp14);
    temp15= new double[N]; tempVec.push_back(temp15);
    temp16= new double[N]; tempVec.push_back(temp16);
    temp17= new double[N]; tempVec.push_back(temp17);
    temp18= new double[N]; tempVec.push_back(temp18);
    temp19= new double[N]; tempVec.push_back(temp19);
    temp20= new double[N]; tempVec.push_back(temp20);
    temp21= new double[N]; tempVec.push_back(temp21);
    temp22= new double[N]; tempVec.push_back(temp22);
    temp23= new double[N]; tempVec.push_back(temp23);
    temp24= new double[N]; tempVec.push_back(temp24);
    temp25= new double[N]; tempVec.push_back(temp25);
    
    /*
    c2d->allocX(tempX1);  tempVec.push_back(tempX1);
    c2d->allocX(tempX2);  tempVec.push_back(tempX2);
    c2d->allocX(tempX3);  tempVec.push_back(tempX3);
    c2d->allocX(tempX4);  tempVec.push_back(tempX4);
    c2d->allocX(tempX5);  tempVec.push_back(tempX5);
    c2d->allocX(tempX6);  tempVec.push_back(tempX6);
    c2d->allocX(tempX7);  tempVec.push_back(tempX7);
    c2d->allocX(tempX8);  tempVec.push_back(tempX8);
    c2d->allocX(tempX9);  tempVec.push_back(tempX9);
    c2d->allocX(tempX10); tempVec.push_back(tempX10);
*/

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > initSolDat Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/
}


void CurvilinearCSolver::calcDtFromCFL(){

//    if(useTiming) ft1 = MPI_Wtime();    

    //Calculate the wave speed over the local spacings...
   // double *UChar_dx, *muchar, *betachar, *kappachar, *delta;
   // muchar = new double[N];
   // betachar = new double[N];
   // kappachar = new double[N];

    //Get the largest value in the domain
    double lmax_UChar_dx = 100000000.0;
    #pragma omp parallel for reduction(min:lmax_UChar_dx)
    FOR_XY{

	//getting the local spacings
	double dx = dom->dx/J11[ip];	//Are these complete? Assumes orthogonal grid
	double dy = dom->dy/J22[ip];	
	double delta = sqrt(dx*dy);

	//calculating the inviscid characteristic
	double UChar_dx = fabs(U[ip])/dx + fabs(V[ip])/dy + sos[ip]*sqrt((1.0/(dx*dx)) + (1.0/(dy*dy)));
	if(UChar_dx < lmax_UChar_dx){
	    lmax_UChar_dx = UChar_dx;
	}

	double muchar, betachar, kappachar;
	if(LADFlag){
	    double mu_eff   = mu[ip] + lad->mu_star[ip];
	    double beta_eff = lad->beta_star[ip];
	    double k_eff = (ig->cp/ig->Pr)*mu[ip] + lad->k_star[ip]; 

	    muchar    = rho1[ip]*delta*delta/mu_eff; 
	    betachar  = rho1[ip]*delta*delta/beta_eff;
	    kappachar = rho1[ip]*sos[ip]*sos[ip]*delta*delta/(k_eff*T[ip]);

	    //TODO finish this up... 
	}

    }
    
    double max_UChar_dx = 1.0/lmax_UChar_dx;

    if(ts->timeSteppingType==Options::CONST_CFL){
	ts->dt = ts->CFL*max_UChar_dx;
    }else if(ts->timeSteppingType==Options::CONST_DT){
	ts->CFL = ts->dt/max_UChar_dx;
    }
  
    if(timeStep == opt->timeStep){
	timeStep++;
	time = opt->time;
	cout << endl;
    }else{
	timeStep++;
	time += ts->dt;
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > calcDtCFL  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/


}

void CurvilinearCSolver::computeGradient(vector<double*> vecIn, vector<double*>vecOut){

    if(vecIn.size() > 4){
	cout << "ERROR:computeGradient: NEED MORE TEMPORARY MEMORY, COMPUTING GRADIENT ON LIST LARGER THAN 4" << endl;
    }else{

	int vecInSize = vecIn.size();

	if(vecOut.size() != vecInSize*2){
	    cout << "ERROR:computeGradient: vecOut size should be 2x the vecIn size!" << endl;
	}

	if(compStyle == VANILLA){

	    //Xi1 Derivatives First
	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi1->calc1stDerivField(vecIn[ip], vecOut[2*ip]);
       	    }

	    //Xi2 Derivatives Next
	    for(int ip = 0; ip < vecInSize; ip++){
	        transposeMatrix_Fast2(vecIn[ip], Nx, Ny, tempVec[ip], opt->blocksize);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        derivXi1->calc1stDerivField(tempVec[ip], tempVec[ip+vecInSize]);
	    }

	    for(int ip = 0; ip < vecInSize; ip++){
	        transposeMatrix_Fast2(tempVec[ip+vecInSize], Ny, Nx, vecOut[2*ip + 1], opt->blocksize);
	    }


	}
    }
};

void CurvilinearCSolver::computeGradDotComponents(vector<double*> vecIn, vector<double*>vecOut){

    if(vecIn.size()/2 > 4){
	cout << "ERROR:computeGradDotComponents: NEED MORE TEMPORARY MEMORY, COMPUTING GRADIENT ON LIST LARGER THAN 4" << endl;
    }else{

	int vecInSize = vecIn.size();
	int N = vecIn.size()/2;

	if(vecOut.size() != vecInSize){
	    cout << "ERROR:computeGradDotComponents: vecOut size should be the same as the vecIn size!" << endl;
	}

	if(N*2 != vecInSize){
	    cout << "ERROR::computeGradDotComponents: missing an awkward number of inputs/outputs" << endl;
	}

        if(compStyle == VANILLA){

	    //Xi1 Derivatives First
	    for(int ip = 0; ip < N; ip++){
	        derivXi1->calc1stDerivField(vecIn[ip], vecOut[ip]);
       	    }

	    //Xi2 Derivatives Next
	    for(int ip = 0; ip < N; ip++){
	        transposeMatrix_Fast2(vecIn[ip+N], Nx, Ny, tempVec[ip], opt->blocksize);
	    }

	    for(int ip = 0; ip < N; ip++){
	        derivXi2->calc1stDerivField(tempVec[ip], tempVec[ip+N]);
	    }

	    for(int ip = 0; ip < N; ip++){
	        transposeMatrix_Fast2(tempVec[ip+N], Ny, Nx, vecOut[ip+N], opt->blocksize);
	    }


	}

    }

};

void CurvilinearCSolver::preStepDerivatives(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP;
    double *rhoUP;
    double *rhoVP;
    double *rhoEP;

    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoEP = rhoE1; 
    }else{
	rhoP  = rhok;
	rhoUP = rhoUk;
	rhoVP = rhoVk;
	rhoEP = rhoEk; 
    }

    double ftt1 = 0, ftt2 = 0;

//    if(useTiming) ftt1 = MPI_Wtime();

    /////////////////////////////
    // Calculate the Gradients //
    /////////////////////////////

    double *vi[] = {U, V};
    vector<double*> vecIn(vi, vi+sizeof(vi)/sizeof(vi[0]));
    
    double *vo[] = {dU1, dU2, dV1, dV2};
    vector<double*> vecOut(vo, vo+sizeof(vo)/sizeof(vo[0]));

    computeGradient(vecIn, vecOut);

/*
    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > xidervtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
	ftt1 = MPI_Wtime();
    }
*/
    /////////////////////////////////////////
    // CALC TAU COMPONENTS & EQN COMPONENTS//
    /////////////////////////////////////////

    double *preCont_1, *preCont_2;
    double *preMom1_1, *preMom1_2; 
    double *preMom2_1, *preMom2_2; 
    double *preEngy_1, *preEngy_2; 

    preMom1_1 = temp8;  preMom1_2 = temp9;  
    preMom2_1 = temp10; preMom2_2 = temp11; 
    preEngy_1 = temp12; preEngy_2 = temp13;
    preCont_1 = temp14; preCont_2 = temp15;
    

    //Do LES/Bulk viscosity model stuff here...
    //we'll pre-calculate the velocity tensor...
    double *gradU[2][2];
    gradU[0][0] = temp16; 
    gradU[0][1] = temp17; 
    gradU[1][0] = temp18; 
    gradU[1][1] = temp19; 

    //Now recalculate properties in the new space
    #pragma omp parallel for
    FOR_XY{
	//Calculate our velocity gradient tensor
   	gradU[0][0][ip] = J11[ip]*dU1[ip] + J21[ip]*dU2[ip];
	gradU[0][1][ip] = J12[ip]*dU1[ip] + J22[ip]*dU2[ip];
	gradU[1][0][ip] = J11[ip]*dV1[ip] + J21[ip]*dV2[ip];
	gradU[1][1][ip] = J12[ip]*dV1[ip] + J22[ip]*dV2[ip];

    }
/*
    //getSGSViscosity...
    if(LESFlag){
	if(opt->lesModel == Options::VREMAN){

	    les->getSGSViscosity(gradU, rhoP, rhoUP, rhoVP, rhoWP, rhoEP);

	}else if(opt->lesModel == Options::DSM){

	    les->calcMuSGSTaukk(gradU, rhoP, rhoUP, rhoVP, rhoWP, rhoEP);

	    //include taukk into pressure
	    FOR_XYZ_YPEN{
		//pressure w/ diagonal term
	        p[ip] -= 0.5*les->taukk[ip];

		//Update temperature w/ new pressure
		T[ip] = ig->solveT(rhoP[ip], p[ip]);

		//update viscosity & SOS w/ new temperature
	        mu[ip]  = ig->solveMu(T[ip]);
		sos[ip] = ig->solveSOS(rhoP[ip], p[ip]); 	
	    }
	}
    }
*/

    //Going to do the locally added dissipation stuff here...
    if(LADFlag){
	if(opt->ladModel == Options::KAWAI){
	    lad->calcVelocityTensorStuff(gradU);
	    lad->calcLADViscosity(gradU, rhoP, rhoUP, rhoVP, rhoEP);
	    lad->calcLADBeta(gradU, rhoP, rhoUP, rhoVP, rhoEP);
	    lad->calcLADK(gradU, rhoP, rhoUP, rhoVP, rhoEP);
	}
    }

    double *viT2[] = {T};
    vector<double*> vecInT(viT2, viT2+sizeof(viT2)/sizeof(viT2[0]));
    
    double *voT2[] = {dT1, dT2};
    vector<double*> vecOutT(voT2, voT2+sizeof(voT2)/sizeof(voT2[0]));

    computeGradient(vecInT, vecOutT);

/*
    if(LESFlag){
	if(opt->lesModel == Options::DSM){
	
	    double *gradT[3];
	    c2d->allocY(gradT[0]);
	    c2d->allocY(gradT[1]);
	    c2d->allocY(gradT[2]);

	    FOR_XYZ_YPEN{
		gradT[0][ip] = J11[ip]*dT1[ip] + J21[ip]*dT2[ip] + J31[ip]*dT3[ip];  
		gradT[1][ip] = J12[ip]*dT1[ip] + J22[ip]*dT2[ip] + J32[ip]*dT3[ip];  
		gradT[2][ip] = J13[ip]*dT1[ip] + J23[ip]*dT2[ip] + J33[ip]*dT3[ip];  
	    }

	    les->calcKSGS(gradT, rhoP, rhoUP, rhoVP, rhoWP, T);
	}
    } 
*/

    #pragma omp parallel for
    FOR_XY{

        double b_1, b_2, b_3;
        double q_1, q_2, q_3;

   	double dUdx = gradU[0][0][ip]; 
	double dUdy = gradU[0][1][ip];
	double dVdx = gradU[1][0][ip];
	double dVdy = gradU[1][1][ip];

	//Viscous stress tensor
	Tau11[ip] =  (4.0/3.0)*dUdx - (2.0/3.0)*dVdy;
	Tau22[ip] = -(2.0/3.0)*dUdx + (4.0/3.0)*dVdy;

	Tau12[ip] =  dUdy + dVdx;

	mu_eff[ip] = 0.0, k_eff[ip] = 0.0;
/*
	if(LESFlag){

	    mu_eff = mu[ip] + les->mu_sgs[ip];

	    if(opt->lesModel == Options::DSM){
		 k_eff  = ig->cp*mu[ip]/ig->Pr  + les->k_sgs[ip];   
	    }else{
	         k_eff  = ig->cp*(mu[ip]/ig->Pr + les->mu_sgs[ip]/Pr_t);
	    }

	    //Apply clipping to mu_eff and k_eff
	    if(mu_eff < 0.0){
		mu_eff = 0.0;
	    }

	    if(k_eff < 0.0){
		k_eff = 0.0;
	    }


	}else{
*/
	    mu_eff[ip] = mu[ip];
	    k_eff[ip]  = (ig->cp/ig->Pr)*mu[ip];
//	}


	if(LADFlag){
	    //Should some clipping or limiting be done here?
	    mu_eff[ip] += lad->mu_star[ip];
	    k_eff[ip]  += lad->k_star[ip];
	}

	Tau11[ip] *= mu_eff[ip];
	Tau22[ip] *= mu_eff[ip];
	Tau12[ip] *= mu_eff[ip];

	//Adding back in the LAD bulk viscosity...
	if(LADFlag){
	    Tau11[ip] += lad->beta_star[ip]*(dUdx + dVdy);
	    Tau22[ip] += lad->beta_star[ip]*(dUdx + dVdy);
	}

	//k_sgs = (ig->cp/Pr_t)*mu_sgs;

	//Thermal conduction terms
	q_1 = k_eff[ip]*(J11[ip]*dT1[ip] + J21[ip]*dT2[ip]);
	q_2 = k_eff[ip]*(J12[ip]*dT1[ip] + J22[ip]*dT2[ip]);

	//Viscous work + thermal conduction terms
	b_1 = U[ip]*Tau11[ip] + V[ip]*Tau12[ip] + q_1;
	b_2 = U[ip]*Tau21[ip] + V[ip]*Tau22[ip] + q_2;


	//Continuity Terms
	preCont_1[ip] = (1.0/J[ip])*(rhoP[ip]*Ucurv[ip]);
	preCont_2[ip] = (1.0/J[ip])*(rhoP[ip]*Vcurv[ip]);

	//Combined Euler + viscous terms
	double F1, G1;
	double F2, G2;

	F1 = rhoUP[ip]*Ucurv[ip] + J11[ip]*p[ip];
	F2 = rhoUP[ip]*Vcurv[ip] + J21[ip]*p[ip];

	G1 = J11[ip]*Tau11[ip] + J12[ip]*Tau21[ip];
	G2 = J21[ip]*Tau11[ip] + J22[ip]*Tau21[ip];

	preMom1_1[ip] = (1.0/J[ip])*(F1 - G1);
	preMom1_2[ip] = (1.0/J[ip])*(F2 - G2);

	F1 = rhoVP[ip]*Ucurv[ip] + J12[ip]*p[ip];
	F2 = rhoVP[ip]*Vcurv[ip] + J22[ip]*p[ip];

	G1 = J11[ip]*Tau12[ip] + J12[ip]*Tau22[ip];
	G2 = J21[ip]*Tau12[ip] + J22[ip]*Tau22[ip];

	preMom2_1[ip] = (1.0/J[ip])*(F1 - G1);
	preMom2_2[ip] = (1.0/J[ip])*(F2 - G2);


	F1 = (rhoEP[ip] + p[ip])*Ucurv[ip];
	F2 = (rhoEP[ip] + p[ip])*Vcurv[ip];

	G1 = J11[ip]*b_1 + J12[ip]*b_2;
	G2 = J21[ip]*b_1 + J22[ip]*b_2;

	preEngy_1[ip] = (1.0/J[ip])*(F1 - G1);
	preEngy_2[ip] = (1.0/J[ip])*(F2 - G2);

    }

/*
    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " > Tau&Components Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }
*/
    ///////////////////////////////////////
    // COMPUTE DERIVATIVES OF COMPONENTS //
    ///////////////////////////////////////

    double *vi2[] = {preCont_1, preMom1_1, preMom2_1, preEngy_1, \
		     preCont_2, preMom1_2, preMom2_2, preEngy_2};
    vector<double*> vecIn2(vi2, vi2+sizeof(vi2)/sizeof(vi2[0]));

 
    double *vo2[] = {cont_1, mom1_1, mom2_1, engy_1, \
		     cont_2, mom1_2, mom2_2, engy_2};
    vector<double*> vecOut2(vo2, vo2+sizeof(vo2)/sizeof(vo2[0]));

    computeGradDotComponents(vecIn2, vecOut2);

/*
    if(useTiming){
	ftt2 = MPI_Wtime();
	IF_RANK0 cout << " >  derivtrans Timing: " << setw(6)  << (int)((ftt2-ftt1)*1000) << "ms" << endl;
    }
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > preStepDer Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/

}


void CurvilinearCSolver::solveContinuity(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP;
    if(rkStep == 1){
        rhoP = rho1;
    }else{
        rhoP = rhok;
    }

    #pragma omp parallel for
    FOR_XY{
        double spgSource;

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoP[ip], spg->spongeRhoAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;
		
	rhok2[ip]  = ts->dt*J[ip]*(-cont_1[ip] - cont_2[ip] + (spgSource+contRHSSource(ip))/J[ip]);
    }
/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveCont  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/
}

void CurvilinearCSolver::solveXMomentum(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoUP;
    if(rkStep == 1){
        rhoUP = rhoU1;
    }else{
        rhoUP = rhoUk;
    }

    #pragma omp parallel for
    FOR_XY{
        double spgSource;

	if(spongeFlag)
	    spgSource = calcSpongeSource(rhoUP[ip], spg->spongeRhoUAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

	rhoUk2[ip] += -mom1_1[ip] - mom1_2[ip] + (spgSource+xmomRHSSource(ip))/J[ip];
	rhoUk2[ip] *= ts->dt*J[ip];
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveXMom  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/


}

void CurvilinearCSolver::solveYMomentum(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoVP;
    if(rkStep == 1){
        rhoVP = rhoV1;
    }else{
        rhoVP = rhoVk;
    }

    #pragma omp parallel for
    FOR_XY{ 
        double spgSource;

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoVP[ip], spg->spongeRhoVAvg[ip], spg->sigma[ip]);
	else
	    spgSource = 0.0;

	rhoVk2[ip] += -mom2_1[ip] -mom2_2[ip] + (spgSource+ymomRHSSource(ip))/J[ip];
        rhoVk2[ip] *= ts->dt*J[ip];
   }


/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveYMom  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/


}

void CurvilinearCSolver::solveEnergy(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoEP;
    if(rkStep == 1){
        rhoEP = rhoE1;
    }else{
        rhoEP = rhoEk;
    }

    #pragma omp parallel for
    FOR_XY{
        double spgSource;

        if(spongeFlag)
            spgSource = calcSpongeSource(rhoEP[ip], spg->spongeRhoEAvg[ip], spg->sigma[ip]);
        else
            spgSource = 0.0;

	rhoEk2[ip] = -engy_1[ip] - engy_2[ip] + (spgSource+engyRHSSource(ip))/J[ip];
	rhoEk2[ip] *= ts->dt*J[ip];
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > solveEngy  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/

}

void CurvilinearCSolver::filterConservedData(){

//    if(useTiming) ft1 = MPI_Wtime();

    double *rhoP_in, *rhoP_out;
    double *rhoUP_in, *rhoUP_out;
    double *rhoVP_in, *rhoVP_out;
    double *rhoEP_in, *rhoEP_out;
    if(!rkLast){
	//in
        rhoP_in  = rhok;
        rhoUP_in = rhoUk;
        rhoVP_in = rhoVk;
        rhoEP_in = rhoEk;

	//out
        rhoP_out  = rhok;
        rhoUP_out = rhoUk;
        rhoVP_out = rhoVk;
        rhoEP_out = rhoEk;

    }else{
        rhoP_in  = rho2;
        rhoUP_in = rhoU2;
        rhoVP_in = rhoV2;
        rhoEP_in = rhoE2;

	//out
        rhoP_out  = rho1;
        rhoUP_out = rhoU1;
        rhoVP_out = rhoV1;
        rhoEP_out = rhoE1;
    }


    //Need to do round robin of filtering directions...
    if(timeStep%ts->filterStep == 0){

        //Advance the filtering time step       
        filterTimeStep++;

        //Going to try and be cute to minimize dmemory allocation
        if(filterTimeStep%2 == 0){

            //Here we'll do Y->X     
	    transposeMatrix_Fast2(rhoP_in,  Nx, Ny, temp1, opt->blocksize);
	    transposeMatrix_Fast2(rhoUP_in, Nx, Ny, temp2, opt->blocksize);
	    transposeMatrix_Fast2(rhoVP_in, Nx, Ny, temp3, opt->blocksize);
	    transposeMatrix_Fast2(rhoEP_in, Nx, Ny, temp4, opt->blocksize);

            filtY->filterField(temp1, temp5);
	    filtY->filterField(temp2, temp6);
	    filtY->filterField(temp3, temp7);
	    filtY->filterField(temp4, temp8);

	    transposeMatrix_Fast2(temp5, Ny, Nx, temp9,  opt->blocksize); 
	    transposeMatrix_Fast2(temp6, Ny, Nx, temp10, opt->blocksize); 
	    transposeMatrix_Fast2(temp7, Ny, Nx, temp11, opt->blocksize); 
	    transposeMatrix_Fast2(temp8, Ny, Nx, temp12, opt->blocksize); 

	    filtX->filterField(temp9,  rhoP_out);
	    filtX->filterField(temp10, rhoUP_out);
	    filtX->filterField(temp11, rhoVP_out);
	    filtX->filterField(temp12, rhoEP_out);

        }else if(filterTimeStep%2 == 1){

            //Here we'll do X->Y     
	    filtX->filterField(rhoP_in, temp1);
	    filtX->filterField(rhoUP_in, temp2);
	    filtX->filterField(rhoVP_in, temp3);
	    filtX->filterField(rhoEP_in, temp4);

	    transposeMatrix_Fast2(temp1, Nx, Ny, temp5, opt->blocksize);
	    transposeMatrix_Fast2(temp2, Nx, Ny, temp6, opt->blocksize);
	    transposeMatrix_Fast2(temp3, Nx, Ny, temp7, opt->blocksize);
	    transposeMatrix_Fast2(temp4, Nx, Ny, temp8, opt->blocksize);

	    filtY->filterField(temp5, temp9);
	    filtY->filterField(temp6, temp10);
	    filtY->filterField(temp7, temp11);
	    filtY->filterField(temp8, temp12);

	    transposeMatrix_Fast2(temp9,  Ny, Nx, rhoP_out, opt->blocksize);
	    transposeMatrix_Fast2(temp10, Ny, Nx, rhoUP_out, opt->blocksize);
	    transposeMatrix_Fast2(temp11, Ny, Nx, rhoVP_out, opt->blocksize);
	    transposeMatrix_Fast2(temp12, Ny, Nx, rhoEP_out, opt->blocksize);
        }
 
    //If not filtering, need to copy the solution over to the *1 variables
    }else{

	if(rkLast){
	    memcpy(rho1,  rho2,  sizeof(double)*Nx*Ny);
	    memcpy(rhoU1, rhoU2, sizeof(double)*Nx*Ny);
	    memcpy(rhoV1, rhoV2, sizeof(double)*Nx*Ny);
	    memcpy(rhoE1, rhoE2, sizeof(double)*Nx*Ny);
	}
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > filterCons Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/


};


void CurvilinearCSolver::updateNonConservedData(){

//    if(useTiming) ft1 = MPI_Wtime();


    if(!rkLast){

	#pragma omp parallel for
	FOR_XY{
	    U[ip]   = ig->solveU(rhok[ip], rhoUk[ip]);
	    V[ip]   = ig->solveU(rhok[ip], rhoVk[ip]);
	    p[ip]   = ig->solvep(rhok[ip], rhoEk[ip], U[ip], V[ip]);
	    T[ip]   = ig->solveT(rhok[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    sos[ip] = ig->solveSOS(rhok[ip], p[ip]);
            Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip];
            Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
	}

    }else if(rkLast){

	#pragma omp parallel for
	FOR_XY{
	    U[ip]   = ig->solveU(rho1[ip], rhoU1[ip]);
	    V[ip]   = ig->solveU(rho1[ip], rhoV1[ip]);
	    p[ip]   = ig->solvep(rho1[ip], rhoE1[ip], U[ip], V[ip]);
	    T[ip]   = ig->solveT(rho1[ip], p[ip]);
	    mu[ip]  = ig->solveMu(T[ip]);
	    sos[ip] = ig->solveSOS(rho1[ip], p[ip]);
            Ucurv[ip] = J11[ip]*U[ip] + J12[ip]*V[ip];
            Vcurv[ip] = J21[ip]*U[ip] + J22[ip]*V[ip];
	}
    }



/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > updNonCons Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }
*/
/*
    cout << "RHO" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhok[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    } 


    cout << "RHOU" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhoUk[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    }
 
    cout << "RHOV" << endl; 
    FOR_Z_YPEN{
        FOR_Y_YPEN{
            FOR_X_YPEN{
                int ip = GETMAJIND_YPEN;
	 	cout << rhoVk[ip] << " ";
            }
	    cout << endl;
        }
	cout << endl;
    } 
*/
}

void CurvilinearCSolver::checkSolution(){


//    if(useTiming) ft1 = MPI_Wtime();

    if(timeStep%ts->checkStep == 0){

	//t2 = MPI_Wtime();

	    cout << endl;
            cout << "-------------------------------------------------" << endl;
            cout << " Step = "<< timeStep << ", time = " << time << ", dt = " << ts->dt << endl;
            cout << "-------------------------------------------------" << endl;
            cout << "  Time since last timestep = " << t2 - t1  << endl;
	

	if(spongeFlag)
	    getRange(spg->sigma, "SIGMA", Nx, Ny);

        getRange(rho1, "RHO", Nx, Ny);
        getRange(U, "U", Nx, Ny);
        getRange(V, "V", Nx, Ny);
        getRange(Ucurv, "Ucurv", Nx, Ny);
        getRange(Vcurv, "Vcurv", Nx, Ny);
        getRange(p, "P", Nx, Ny);
        getRange(T, "T", Nx, Ny);
        getRange(mu, "mu", Nx, Ny);
        getRange(rhoU1, "RHOU", Nx, Ny);
        getRange(rhoV1, "RHOV", Nx, Ny);
        getRange(rhoE1, "RHOE", Nx, Ny);
        getRange(sos, "SOS", Nx, Ny);
/*
	if(LESFlag){
	    getRange(les->mu_sgs, "MU SGS", Nx, Ny, Nz, mpiRank);
	    if(opt->lesModel == Options::DSM){
	    getRange(les->taukk, "TAUKK ", Nx, Ny, Nz, mpiRank);
	    getRange(les->k_sgs, "K SGS", Nx, Ny, Nz, mpiRank);
		if(opt->dumpCoeffRange){
		    getRange(les->C, "LES C", Nx, Ny, Nz, mpiRank);
		    getRange(les->CI, "LES CI", Nx, Ny, Nz, mpiRank);
		    getRange(les->Prt, "LES Prt", Nx, Ny, Nz, mpiRank);
		}	
	    }
	}

	if(statsFlag){
  	    if(stats->velocityStats){
	        getRange(stats->UAVG, "UAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->URMS, "URMS", Nx, Ny, Nz, mpiRank);
	        getRange(stats->VAVG, "VAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->VRMS, "VRMS", Nx, Ny, Nz, mpiRank);
	        getRange(stats->WAVG, "WAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->WRMS, "WRMS", Nx, Ny, Nz, mpiRank);
	        getRange(stats->UVREY, "UVREY", Nx, Ny, Nz, mpiRank);
	        getRange(stats->UWREY, "UWREY", Nx, Ny, Nz, mpiRank);
	        getRange(stats->VWREY, "VWREY", Nx, Ny, Nz, mpiRank);
	    }

	    if(stats->thermoStats){
	        getRange(stats->RHOAVG, "RHOAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->RHORMS, "RHORMS", Nx, Ny, Nz, mpiRank);
	        getRange(stats->PAVG, "PAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->PRMS, "PRMS", Nx, Ny, Nz, mpiRank);
	        getRange(stats->TAVG, "TAVG", Nx, Ny, Nz, mpiRank);
	        getRange(stats->TRMS, "TRMS", Nx, Ny, Nz, mpiRank);
	    }
	}


	IF_RANK0 cout << endl;

	t1 = MPI_Wtime();
*/
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > checkSoln  Timing: " << setw(6)  << (int)((ft2-ft1)*1000) << "ms" << endl;
    }	 
*/

};


void CurvilinearCSolver::dumpSolution(){
/*
//    if(useTiming) ft1 = MPI_Wtime();

	double time1 = MPI_Wtime();

	//Converting over to X-Major Matrices
	FOR_Z_YPEN{
	    FOR_Y_YPEN{
		FOR_X_YPEN{
		    int ip = GETMAJIND_YPEN;
		    int jp = GETIND_YPEN;

		    tempY1[jp] = rho1[ip];
		    tempY2[jp] = rhoU1[ip];
		    tempY3[jp] = rhoV1[ip];
		    tempY4[jp] = rhoW1[ip];
		    tempY5[jp] = rhoE1[ip];

		    tempY6[jp] = msh->x[ip];
		    tempY7[jp] = msh->y[ip];
		    tempY8[jp] = msh->z[ip];

		}
	    }
	}
	
	


	IF_RANK0{
            cout << endl;
            cout << " > ===============" << endl;
            cout << " >  DUMPING FIELD " << endl;
            cout << " > ===============" << endl;
	}

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
	outputFileName = "SolutionDump.";
	
	ostringstream timeStepString;
        timeStepString << timeStep;

	outputFileName.append(timeStepString.str());

	MPI_File fh;
	MPI_Offset disp, filesize;

	MPI_File_open(MPI_COMM_WORLD, outputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	filesize = 0;
	MPI_File_set_size(fh, filesize);

	disp = 0;

	double cTS, cNx, cNy, cNz;
	cTS = (double)timeStep;
	cNx = (double)Nx;	
	cNy = (double)Ny;	
	cNz = (double)Nz;	

	c2d->writeScalar(fh, disp, 1, &cTS);
	c2d->writeScalar(fh, disp, 1, &time);
	c2d->writeScalar(fh, disp, 1, &cNx); 
	c2d->writeScalar(fh, disp, 1, &cNy); 
	c2d->writeScalar(fh, disp, 1, &cNz); 
	c2d->writeVar(fh, disp, baseDirection, tempY6);
	c2d->writeVar(fh, disp, baseDirection, tempY7);
	c2d->writeVar(fh, disp, baseDirection, tempY8);
	c2d->writeVar(fh, disp, baseDirection, tempY1);
	c2d->writeVar(fh, disp, baseDirection, tempY2);
	c2d->writeVar(fh, disp, baseDirection, tempY3);
	c2d->writeVar(fh, disp, baseDirection, tempY4);
	c2d->writeVar(fh, disp, baseDirection, tempY5);

	MPI_File_close(&fh);


	if(spongeFlag){
	    spg->dumpSpongeAvg(timeStep);
	}

        if(statsFlag){
	    stats->dumpStatsFields();
	}

	double time2 = MPI_Wtime();

	IF_RANK0{
	    cout << endl;
	    cout << " > File dump took " << time2-time1 << " seconds." << endl;
	}
/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > dumpSoln Timing: " << setw(6)  << endl;
	
    }	 
*/


}

void CurvilinearCSolver::addImageOutput(PngWriter *png){
    imageList.push_back(png);
}

void CurvilinearCSolver::writeImages(){

//    if(useTiming) ft1 = MPI_Wtime();

    //Get the timestep string just in case...
    ostringstream timeStepStringT;
    timeStepStringT << timeStep;
    int zeroPad = 6;
    string timeStepString = string(zeroPad - (timeStepStringT.str()).length(), '0') + timeStepStringT.str();

    //iterating through the image list
    for(list<PngWriter*>::iterator iter=imageList.begin(); iter != imageList.end(); ++iter){

	(*iter)->timeStepString = timeStepString;

	if(timeStep == opt->timeStep && (*iter)->dumpInterval == 0){

	    //if we haven't gotten an interpolator for this pngWriter yet, lets generate it
	    if((*iter)->interpolatorFlag == false){
		//Call the function to generate it...
		generateImagePlane(*iter); 	
		(*iter)->interpolatorFlag = true;
	    }

	    //Do the write...
	    writePlaneImageForVariable(*iter);	    
	}


	if(timeStep%(*iter)->dumpInterval==0 && (*iter)->dumpInterval > 0){

	    //if we haven't gotten an interpolator for this pngWriter yet, lets generate it
	    if((*iter)->interpolatorFlag == false){
		//Call the function to generate it...
		generateImagePlane(*iter); 	
		(*iter)->interpolatorFlag = true;
	    }

	    //Do the write
	    writePlaneImageForVariable(*iter);	    
	}
    }

/*
    if(useTiming){
	ft2 = MPI_Wtime();
	IF_RANK0 cout << " > writeImg Timing: " << setw(6)  << endl;
    }	 
*/


}

void CurvilinearCSolver::generateImagePlane(PngWriter *pw){

    int N1 = pw->nx, N2 = pw->ny; 
    double (*pointList)[2] = NULL;

    pointList = new double[N1*N2][2];

    //See which case is the limited box size
    double d1, d2;
    if(pw->bboxFlag){
	d1 = (pw->x_max[0] - pw->x_min[0])/((double)N1 + 1.0);
	d2 = (pw->x_max[1] - pw->x_min[1])/((double)N2 + 1.0);
    }else{
	d1 = (msh->x_max - msh->x_min)/((double)N1 + 1.0);
	d2 = (msh->y_max - msh->y_min)/((double)N2 + 1.0);
    }

    double dx = fmax(d1, d2);

    //Get the base locations
    double base1, base2;
    if(pw->bboxFlag){
	base1 = pw->x_min[0] + dx/2;
	base2 = pw->x_min[1] + dx/2;

	if((pw->x_max[1]-pw->x_min[1])>(pw->x_max[0]-pw->x_min[0])){
	    base1 -= (dx*(N1-1))/2.0-(pw->x_max[0]- pw->x_min[0])/2.0;
	}else{
	    base2 -= (dx*(N2-1))/2.0-(pw->x_max[1]- pw->x_min[1])/2.0;
	}

    }else{
	base1 = msh->x_min + dx/2;
	base2 = msh->y_min + dx/2;
    }

    //Now calculate the positions of the pixels
    for(int jp = 0; jp < N2; jp++){
        for(int ip = 0; ip < N1; ip++){
	    int ii = jp*N1 + ip;
	    pointList[ii][0] = base1 + dx*(double)ip;
	    pointList[ii][1] = base2 + dx*(double)jp;	
	}
    }

    cout << pointList[(N2-1)*N1 + N1-1][0] << " " << pointList[(N2-1)*N1 + N1-1][1] << endl;

    pw->ci = new CurvilinearInterpolator(this, pointList, N1*N2);

    delete[] pointList;
}

void CurvilinearCSolver::writePlaneImageForVariable(PngWriter *pw){

    int N1 = pw->nx;
    int N2 = pw->ny;
    double *var = pw->fieldPtr;

    double *ff_ci = new double[N1*N2];
    double *ff    = new double[N1*N2];
 
    for(int ip = 0; ip < N1*N2; ip++){
	ff_ci[ip] = -1000000.0;
	ff[ip]    = -1000000.0;
    }

    pw->ci->interpolateData(var, ff_ci);
 
    for(int ip = 0; ip < pw->ci->pointFoundCount; ip++){
	ff[pw->ci->pointIndex[ip]] = ff_ci[ip];
    }
 
    double dataMin =  100000000.0;
    double dataMax = -100000000.0;


    if(pw->valFlag == true){
	dataMin = pw->valMin;
	dataMax = pw->valMax;
    }else{
	//get the max/min value in the plane...
	for(int ip = 0; ip < N1*N2; ip++){
	    if(ff[ip] > -100000.0){
	        dataMin = fmin(dataMin, ff[ip]);
	        dataMax = fmax(dataMax, ff[ip]);
	    }
	}
    }

    //Scale pixel value to local min and max
    int *r = new int[N1*N2];
    int *g = new int[N1*N2];
    int *b = new int[N1*N2];
    for(int ip = 0; ip < N1*N2; ip++){

	double phitemp = (ff[ip] - dataMin)/(dataMax - dataMin); 

	if((dataMax - dataMin) < 1E-6){
	    phitemp = 0.0;
	}

	if(phitemp > 1.0){
	    phitemp = 1.0;
	}

	if(phitemp < 0.0){
	   phitemp = 0.0;
	}

	if(pw->cm == PngWriter::BWR){
	    pw->getPARAVIEW_BWR(phitemp, r[ip], g[ip], b[ip]);
	}else if(pw->cm == PngWriter::RAINBOW){
	    pw->getRainbowColormap(phitemp, r[ip], g[ip], b[ip]);
	}else if(pw->cm == PngWriter::GREYSCALE){ 
	    r[ip] = (int)(phitemp*255.0);
	    g[ip] = (int)(phitemp*255.0);
	    b[ip] = (int)(phitemp*255.0);
	}	  

    }

    for(int jp = 0; jp < N2; jp++){
	for(int ip = 0; ip < N1; ip++){
	    int ii = jp*N1 + ip;
	    //Grayscale image...
	    if(ff[ii] > -100000.0){
		pw->set(ip, jp, r[ii], g[ii], b[ii]);
	    }else{
		pw->set(ip, jp, 73, 175, 205);
	    }
	}
    }

    string imageName = pw->varName;
    imageName.append(".");
    imageName.append(pw->timeStepString);
    imageName.append(".png");
    pw->write(imageName.c_str());

    delete[] r;
    delete[] g;
    delete[] b;
    delete[] ff;

}

bool CurvilinearCSolver::checkForAndDeleteKillFile(string killFileName){

    bool fileExists;

        if(FILE *file = fopen(killFileName.c_str(),"r")){
            fclose(file);
	    fileExists = true;
        }else{
	    fileExists = false;
        }

        if(fileExists){
	    if(remove(killFileName.c_str()) != 0){
	        cout << " > ERROR DELETING KILLFILE! " << endl;
	    }
        }

    return fileExists;
}

void CurvilinearCSolver::checkEnd(){

    if(checkForAndDeleteKillFile("killsolver")){

	    cout << "===================" << endl;
	    cout << " FOUND KILL SOLVER " << endl;
	    cout << "===================" << endl;
	
	endFlag = true;
    }

    if(time >= ts->maxTime){

	    cout << "=================" << endl;
	    cout << " HIT END OF TIME " << endl;
	    cout << "=================" << endl;

	endFlag = true;
    }

    if(timeStep >= ts->maxTimeStep){
	
	    cout << "=================" << endl;
	    cout << " HIT END OF TIME " << endl;
	    cout << "=================" << endl;

	endFlag = true;

    } 

    if(endFlag){
	dumpSolution();
    }

}

void CurvilinearCSolver::reportAll(){

   cout << "REPORT ALL" << endl;

   getRange(dU1, "dU1", Nx, Ny);
   getRange(dU2, "dU2", Nx, Ny);
   cout << " " << endl;
   getRange(dV1, "dV1", Nx, Ny);
   getRange(dV2, "dV2", Nx, Ny);
   cout << " " << endl;
   getRange(dT1, "dT1", Nx, Ny);
   getRange(dT2, "dT2", Nx, Ny);
   cout << " " << endl;
   getRange(Tau11, "Tau11", Nx, Ny);
   getRange(Tau22, "Tau22", Nx, Ny);
   getRange(Tau12, "Tau12", Nx, Ny);
   cout << " " << endl;
   getRange(cont_1, "cont_1", Nx, Ny);
   getRange(cont_2, "cont_2", Nx, Ny);
   cout << " " << endl;
   getRange(mom1_1, "mom1_1", Nx, Ny);
   getRange(mom1_2, "mom1_2", Nx, Ny);
   cout << " " << endl;
   getRange(mom2_1, "mom2_1", Nx, Ny);
   getRange(mom2_2, "mom2_2", Nx, Ny);
   cout << " " << endl;
   getRange(engy_1, "engy_1", Nx, Ny);
   getRange(engy_2, "engy_2", Nx, Ny);
   cout << " " << endl;
   getRange(rho1, "rho1", Nx, Ny);
   getRange(rhok, "rhok", Nx, Ny);
   getRange(rhok2, "rhok2", Nx, Ny);
   getRange(rho2, "rho2", Nx, Ny);
   cout << " " << endl;
   getRange(rhoU1, "rhoU1", Nx, Ny);
   getRange(rhoUk, "rhoUk", Nx, Ny);
   getRange(rhoUk2, "rhoUk2", Nx, Ny);
   getRange(rhoU2, "rhoU2", Nx, Ny);
   cout << " " << endl;
   getRange(rhoV1, "rhoV1", Nx, Ny);
   getRange(rhoVk, "rhoVk", Nx, Ny);
   getRange(rhoVk2, "rhoVk2", Nx, Ny);
   getRange(rhoV2, "rhoV2", Nx, Ny);
   cout << " " << endl;
   getRange(rhoE1, "rhoE1", Nx, Ny);
   getRange(rhoEk, "rhoEk", Nx, Ny);
   getRange(rhoEk2, "rhoEk2", Nx, Ny);
   getRange(rhoE2, "rhoE2", Nx, Ny);
   cout << " " << endl;
   getRange(p, "p", Nx, Ny);
   getRange(U, "U", Nx, Ny);
   getRange(V, "V", Nx, Ny);
   getRange(Ucurv, "Ucurv", Nx, Ny);
   getRange(Vcurv, "Vcurv", Nx, Ny);
   getRange(T, "T", Nx, Ny);
   getRange(mu, "mu", Nx, Ny);
   getRange(sos, "sos", Nx, Ny);
   cout << " " << endl;

}


