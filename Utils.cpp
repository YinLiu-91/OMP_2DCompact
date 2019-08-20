#include "Utils.hpp"


void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int size)
{
        memcpy(work, b, size*sizeof(double));
        memcpy(x, d, size*sizeof(double));

        for(int ip = 1; ip < size; ip++){
            double m = a[ip]/work[ip-1];
            work[ip] = work[ip] - m*c[ip-1];
            x[ip] = x[ip] - m*x[ip-1];
        }

        x[size-1] /= work[size-1];
        for(int ip = size-2; ip >= 0; ip--){
            x[ip] = (x[ip] - c[ip]*x[ip+1])/work[ip];
        }

}

void cyclic(double *a, double *b, double *c, double alpha, double beta, double *r, int n, double *x)
{
        unsigned long i;
        double fact,gamma,bb[n],u[n],z[n],work[n];

        if (n <= 2) cout << "n too small in cyclic" << endl;

        gamma = -b[0];
        for (i=0;i<n;i++) bb[i]=b[i];
        bb[0]=b[0]-gamma;
        bb[n-1]=b[n-1]-alpha*beta/gamma;

        solveTri(a,bb,c,r,x,work,n);

        for (i=0;i<n;i++) u[i]=0.0;
        u[0]=gamma;
        u[n-1]=alpha;

        solveTri(a,bb,c,u,z,work,n);

        fact=(x[0]+beta*x[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);

        for (i=0;i<n;i++) x[i] -= fact*z[i];

}

void solvePenta(double e[], double a[], double d[], double c[], double f[], double b[], double x[], double B[], double C[], double D[], int size){

    //e is the left-most band (2:size-1), but is defined as e(0) to e(size-3)
    //a is the second to left-most band (1:size-1), but is defined as a(0) to a(size-2)
    //d is the diagonal (0:size-1)
    //c is the second to right-most band (0:size-2), but the last element isn't used
    //f is the right-most band (0:size-3), but the last two elements aren't used
    //b is the RHS 
    //x is the solution
    // B, C, D are all just pre-allocated work arrays for the method

    memcpy(x, a, size*sizeof(double));
    memcpy(B, b, size*sizeof(double));
    memcpy(C, c, size*sizeof(double));
    memcpy(D, d, size*sizeof(double));
    double xmult;
    for(int ip = 1; ip < size-1; ip++){
	xmult   = x[ip-1]/D[ip-1];
	D[ip]   = D[ip] - xmult*C[ip-1];
	C[ip]   = C[ip] - xmult*f[ip-1];
	B[ip]   = B[ip] - xmult*B[ip-1];
	xmult   = e[ip-1]/D[ip-1];
	x[ip]   = x[ip] - xmult*C[ip-1];
	D[ip+1] = D[ip+1] - xmult*f[ip-1];
	B[ip+1] = B[ip+1] - xmult*B[ip-1];
    }
    xmult = x[size-2]/D[size-2];
    D[size-1] =  D[size-1] - xmult*C[size-2];
    x[size-1] = (B[size-1] - xmult*B[size-2])/D[size-1];
    x[size-2] = (B[size-2] - C[size-2]*x[size-1])/D[size-2];


    for(int ip = size-3; ip >= 0; ip--){
	x[ip] = (B[ip] - f[ip]*x[ip+2] - C[ip]*x[ip+1])/D[ip];
    }
}

void cyclicPenta(double a[], double b[], double c[], double d[], double e[], double f[], double cp[6], double x[], int N){

//Solves for x[0:(N-1)] in a pentadiagonal matrix with nonzero entries in the corners
//
//  x     x     x   0   0   ... 0   cp[0] cp[1]
//  x     x     x   x   0   ... 0   0     cp[2]
//  .     .     .   .   .   .   .   .     .
//  cp[4] 0     0   ... 0   0   x   x     x 
//  cp[5] cp[6] 0   ... 0   x   x   x     x

//NOTE: simple change that should improve performance some would be to preallocate these arrays and pass them to the function instead
//of allocating them every time the function is called

	double u[4][N];
	double v[4][N];
	double z[4][N];
	double work1[N], work2[N], work3[N], y[N];

	//These probably don't need to be preallocated for speed because they're small
	double h[4][4], p[4][4];

	for(int ip = 0; ip < 4; ip++){
	    for(int jp = 0; jp < N; jp++){
		u[ip][jp] = 0.0;
		v[ip][jp] = 0.0;
		z[ip][jp] = 0.0;
	    }
	}

	u[0][0]   = 1.0;
	u[1][1]   = 1.0;
	u[2][N-2] = 1.0;
	u[3][N-1] = 1.0;

	v[0][N-2] = cp[0];
	v[0][N-1] = cp[1];
	v[1][N-1] = cp[2];
	v[2][0]   = cp[3];
	v[3][0]   = cp[4];
	v[3][1]   = cp[5];

	solvePenta(a, b, c, d, e, u[0], z[0], work1, work2, work3, N);
	solvePenta(a, b, c, d, e, u[1], z[1], work1, work2, work3, N);
	solvePenta(a, b, c, d, e, u[2], z[2], work1, work2, work3, N);
	solvePenta(a, b, c, d, e, u[3], z[3], work1, work2, work3, N);
	solvePenta(a, b, c, d, e, f,    y,    work1, work2, work3, N);

 	for(int ip = 0; ip < 4; ip++){
	    for(int jp = 0; jp < 4; jp++){
		double sum = 0.0;
		for(int kp = 0; kp < N; kp++){
		    sum += v[ip][kp]*z[jp][kp];
		}
		p[ip][jp] = sum;
	    }
	}
	for(int ip = 0; ip < 4; ip++){
	    p[ip][ip] += 1.0;
	}

	//need to do inverse of matrix...
	double pr[16], hr[16];
	for(int ip = 0; ip < 4; ip++){
	    for(int jp = 0; jp < 4; jp++){
		pr[ip*4 + jp] = p[ip][jp]; // is this backwards need to test
	    }
	}
	invertRowMajor(pr, hr);


	for(int ip = 0; ip < 4; ip++){
	    for(int jp = 0; jp < 4; jp++){
		h[ip][jp] = hr[ip*4 + jp];
	    }
	}

	//luinv(p, 4, 4, indx, h);
	//
	double r[4];
	for(int ip = 0; ip < 4; ip++){
	    r[ip] = 0.0;
	    for(int jp = 0; jp < N; jp++){
		r[ip] += v[ip][jp]*y[jp];
	    } 
	}

	double s[4];
	for(int ip = 0; ip < 4; ip++){
	    s[ip] = 0.0;
	    for(int jp = 0; jp < 4; jp++){
		s[ip] += h[ip][jp]*r[jp];
	    }
	}




	for(int ip = 0; ip < N; ip++){
	    double sum = 0.0;
	    for(int jp = 0; jp < 4; jp++){
		sum += z[jp][ip]*s[jp];
	    }
	    x[ip] = y[ip] - sum;
	}


}	

double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2)
{
	    return m[4*r0+c0] * (m[4*r1+c1] * m[4*r2+c2] - m[4*r2+c1] * m[4*r1+c2]) -
		   m[4*r0+c1] * (m[4*r1+c0] * m[4*r2+c2] - m[4*r2+c0] * m[4*r1+c2]) +
		   m[4*r0+c2] * (m[4*r1+c0] * m[4*r2+c1] - m[4*r2+c0] * m[4*r1+c1]);
}
 
 
void adjoint(double m[16], double adjOut[16])
{
    adjOut[ 0] =  MINOR(m,1,2,3,1,2,3); adjOut[ 1] = -MINOR(m,0,2,3,1,2,3); adjOut[ 2] =  MINOR(m,0,1,3,1,2,3); adjOut[ 3] = -MINOR(m,0,1,2,1,2,3);
    adjOut[ 4] = -MINOR(m,1,2,3,0,2,3); adjOut[ 5] =  MINOR(m,0,2,3,0,2,3); adjOut[ 6] = -MINOR(m,0,1,3,0,2,3); adjOut[ 7] =  MINOR(m,0,1,2,0,2,3);
    adjOut[ 8] =  MINOR(m,1,2,3,0,1,3); adjOut[ 9] = -MINOR(m,0,2,3,0,1,3); adjOut[10] =  MINOR(m,0,1,3,0,1,3); adjOut[11] = -MINOR(m,0,1,2,0,1,3);
    adjOut[12] = -MINOR(m,1,2,3,0,1,2); adjOut[13] =  MINOR(m,0,2,3,0,1,2); adjOut[14] = -MINOR(m,0,1,3,0,1,2); adjOut[15] =  MINOR(m,0,1,2,0,1,2);
}
 
double det(double m[16])
{
    return m[0] * MINOR(m, 1, 2, 3, 1, 2, 3) -
	   m[1] * MINOR(m, 1, 2, 3, 0, 2, 3) +
	   m[2] * MINOR(m, 1, 2, 3, 0, 1, 3) -
	   m[3] * MINOR(m, 1, 2, 3, 0, 1, 2);
}
 
 
void invertRowMajor(double m[16], double invOut[16])
{
    adjoint(m, invOut);    
    double inv_det = 1.0f / det(m);
    for(int i = 0; i < 16; ++i)
        invOut[i] = invOut[i] * inv_det;
}


void transposeMatrix(double *in, int Nx, int Ny, double *out){

    for(int i = 0; i < Ny; ++i)
        for(int j = 0; j < Nx; ++j)
            out[i*Nx + j] =  in[j*Ny + i];

}

void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block){
    for (int i = 0; i < n; i += block) {
        for(int j = 0; j < n; ++j) {
            for(int b = 0; b < block && i + b < n; ++b) {
                out[j*n + i + b] = in[(i + b)*n + j];
            }
        }
    }
}

void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize){
    int i, j, row, col;

    #pragma omp parallel for private(i,j, row, col) collapse(2)
    for ( i = 0; i < n; i += blocksize) {
        for ( j = 0; j < p; j += blocksize) {
            for (row = i; row < i + blocksize && row < n; row++) {
                for (col = j; col < j + blocksize && col < p; col++) {
                    out[row*p + col] = in[col*n + row];
                }
            }
        }
    }
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void getRange(double *phi, string Var, int Nx, int Ny){

    double dataMin = 1000000.0;
    double dataMax = -2000000.0;

    #pragma omp parallel for reduction(max:dataMax) reduction( min:dataMin)
    for(int ip = 0; ip < Nx*Ny; ip++){
	dataMin = fmin(phi[ip], dataMin);
	dataMax = fmax(phi[ip], dataMax);
    }


    cout << "   Range of " << Var << ": " << dataMin << ":" << dataMax << endl;

}

void getRangeValue(double *phi, int Nx, int Ny, double &dataMin, double &dataMax){

    dataMin = 1000000.0;
    for(int ip = 0; ip < Nx*Ny; ip++)
	dataMin = min(phi[ip], dataMin);

    dataMax = -1000000.0;
    for(int ip = 0; ip < Nx*Ny; ip++)
	dataMax = max(phi[ip], dataMax);

}

bool isPointInHexa(double p[3], double vertex[8][3]){


    //Do bounding box check first
    double x_max = -10000000.0;
    double x_min =  10000000.0;
    double y_max = -10000000.0;
    double y_min =  10000000.0;
    double z_max = -10000000.0;
    double z_min =  10000000.0;

    for(int ip = 0; ip < 8; ip++){

	x_max = max(x_max, vertex[ip][0]);
	y_max = max(y_max, vertex[ip][1]);
	z_max = max(z_max, vertex[ip][2]);

	x_min = min(x_min, vertex[ip][0]);
	y_min = min(y_min, vertex[ip][1]);
	z_min = min(z_min, vertex[ip][2]);

    }

    bool isInBB = false;
    if(p[0] <= x_max && p[0] >= x_min){
	if(p[1] <= y_max && p[1] >= y_min){
	    if(p[2] <= z_max && p[2] >= z_min){
		isInBB = true;
	    }
	}
    }

    isInBB = true;
    //We are in the bounding box, continue with volume check
    if(isInBB){

	double volWithCentroid = getHexaVolume(vertex);
	double volWithPoint    = getHexaVolumeWithPoint(vertex, p);

	double eps = 1e-12;

	if(fabs(volWithPoint - volWithCentroid) < eps){
	    return true;
	}else{	
	    return false;
	}


    }else{
	return false;
    }

}

double getHexaVolume(double vertex[8][3]){

    //get hexa center...only guaranteed for convex polyhedra?
    double h_center[3] = {0,0,0};
    for(int ip = 0; ip < 8; ip++){
        FOR_I3 h_center[i] += vertex[ip][i];
    }
    FOR_I3 h_center[i] /= 8.0;

    //THESE ARE IN A SPECIFIC ORDER!!!
    double *A = vertex[0]; //{0.0, 0.0, 0.0} 
    double *B = vertex[2]; //{0.0, 1.0, 0.0}
    double *C = vertex[4]; //{1.0, 0.0, 0.0}
    double *D = vertex[6]; //{1.0, 1.0, 0.0}
    double *E = vertex[1]; //{0.0, 0.0, 1.0}
    double *F = vertex[3]; //{0.0, 1.0, 1.0}
    double *G = vertex[5]; //{1.0, 0.0, 1.0}
    double *H = vertex[7]; //{1.0, 1.0, 1.0}

    double volume_hexa = 0.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, B, C, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, D, C, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(G, F, E, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, G, H, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(E, B, A, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, E, F, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, D, G, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(H, G, D, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, C, G, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, G, E, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, D, B, h_center))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(D, F, H, h_center))/6.0;

    return volume_hexa;
};

double getHexaVolumeWithPoint(double vertex[8][3], double P[3]){

    //THESE ARE IN A SPECIFIC ORDER!!!
    double *A = vertex[0]; //{0.0, 0.0, 0.0} 
    double *B = vertex[2]; //{0.0, 1.0, 0.0}
    double *C = vertex[4]; //{1.0, 0.0, 0.0}
    double *D = vertex[6]; //{1.0, 1.0, 0.0}
    double *E = vertex[1]; //{0.0, 0.0, 1.0}
    double *F = vertex[3]; //{0.0, 1.0, 1.0}
    double *G = vertex[5]; //{1.0, 0.0, 1.0}
    double *H = vertex[7]; //{1.0, 1.0, 1.0}

    double volume_hexa = 0.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, B, C, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, D, C, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(G, F, E, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, G, H, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(E, B, A, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(B, E, F, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, D, G, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(H, G, D, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(A, C, G, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(C, G, E, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(F, D, B, P))/6.0;
    volume_hexa += fabs(SIGNED_TET_VOLUME_6(D, F, H, P))/6.0;



    return volume_hexa;
};


void getDataFromIndex(double *f_in, int *index, int N, double *f_out){

    for(int ip = 0; ip < N; ip++){
	f_out[ip] = f_in[index[ip]];
    }

};
