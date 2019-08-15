#ifndef _CALGEBRAICSINGLEBLOCKMESHH_
#define _CALGEBRAICSINGLEBLOCKMESHH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "AbstractSingleBlockMesh.hpp"

class AlgebraicSingleBlockMesh:public AbstractSingleBlockMesh{

    public:

	AlgebraicSingleBlockMesh(AbstractCSolver *cs, Domain *dom){

	    this->cs = cs;
	    this->d = dom;
	    this->derivX = cs->derivX;
	    this->derivY = cs->derivY;

	    max_xi  = d->Lx;
	    max_eta = d->Ly;

	    Nx = d->Nx;
	    Ny = d->Ny;

   	    cout << endl << " > Running Algebraic single block mesh initialization for curvlinear grid solver! " << endl;
	    
	    //Allocate the memory we need for the mesh...
	    allocateForMesh();

	
	    //Generate the mesh...
	    getMesh();
	    

	    //Read some info out to the console
	    getRange(x, "X", Nx, Ny);
	    getRange(y, "Y", Nx, Ny);
	     cout << "   Range of Xi: "  << d->x[0] <<  ":" << d->x[Nx-1] << endl;
	     cout << "   Range of Eta: " << d->y[0] <<  ":" << d->y[Ny-1] << endl;

	    //Handle the periodic boundary conditions coming in...
	    handlePeriodicStuff();

	    //Do the real workhorse of this class...
	    solveForJacobians();

	    //Initialize the ADT...
	    //initMeshADT();

	}

	void getMesh();

	void generateNozzleGrid(double x_in[3], double x_out[3]); 

	void generateCylinderGrid(double xi_in[3], int cylRIndex, int cylRMax, int cylThetaIndex, int cylThetaMax, double x_out[3]);
};


//Where the magic happens...
void AlgebraicSingleBlockMesh::getMesh(){

	    //Get the global mins/maxes
	    double x_min = 1e10, y_min = 1e10;
	    double x_max =-1e10, y_max =-1e10;
	    bool readInGrid = cs->opt->fromRestart || cs->opt->onlyGridFromRestart;
  
	    if(readInGrid){

/*
	        MPI_File fh;
		MPI_Offset disp, filesize;

		string filename = cs->opt->filename;
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		disp = 0;

		double *x_temp,
		       *y_temp,
		       *z_temp;

		c2d->allocY(x_temp);
		c2d->allocY(y_temp);
		c2d->allocY(z_temp);

	 	//Going to displace by the header with the Nx, Ny, Nz sizes
		disp += 40;

		c2d->readVar(fh, disp, 1, x_temp);
		c2d->readVar(fh, disp, 1, y_temp);
		c2d->readVar(fh, disp, 1, z_temp);

		MPI_File_close(&fh);

		FOR_Z_YPEN{
		    FOR_Y_YPEN{
			FOR_X_YPEN{
			    int ip = GETMAJIND_YPEN;
			    int jp = GETIND_YPEN;

			    x[ip] = x_temp[jp];
			    y[ip] = y_temp[jp];
			    z[ip] = z_temp[jp];

			    x_min_local[0] = fmin(x_min_local[0], x[ip]);
			    x_min_local[1] = fmin(x_min_local[1], y[ip]);
			    x_min_local[2] = fmin(x_min_local[2], z[ip]);
	
			    x_max_local[0] = fmax(x_max_local[0], x[ip]);
			    x_max_local[1] = fmax(x_max_local[1], y[ip]);
			    x_max_local[2] = fmax(x_max_local[2], z[ip]);
			}
		    }
		}

		c2d->deallocXYZ(x_temp);
		c2d->deallocXYZ(y_temp);
		c2d->deallocXYZ(z_temp);
*/

	    }else{
	         //Generate the mesh algebraically...

		    #pragma omp parallel for reduction(min:x_min, y_min) reduction(max:y_max,x_max)
		    FOR_Y{
		        FOR_X{
			    int ip = GET2DINDEX_XY;
		
			    double xi  = d->x[i];
			    double eta = d->y[j];
		
			    double xi_in[2] = {xi, eta};
			    double x_out[2];


			    //generateCylinderGrid(xi_in, jj, Ny, ii, Nx, x_out);
			    //generateNozzleGrid(xi_in, x_out);

			    x[ip] = 2.0*M_PI*xi_in[0];
			    y[ip] = 2.0*M_PI*xi_in[1];

//			    x[ip] = xi_in[0];
//			    y[ip] = xi_in[1];
//			    z[ip] = xi_in[2];

			    //Since we're already in this loop, calculate the local max and mins
			    x_min = fmin(x_min, x[ip]);
			    y_min = fmin(y_min, y[ip]);
	
			    x_max = fmax(x_max, x[ip]);
			    y_max = fmax(x_max, y[ip]);
 
		        }
		    }
	    }

	    //Reduce to get the global bounding box
	    //MPI_Allreduce(x_min_local, x_min, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	    //MPI_Allreduce(x_max_local, x_max, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	
};

void AlgebraicSingleBlockMesh::generateNozzleGrid(double xi_in[3], double x_out[3]){

	double xi = xi_in[0], eta = xi_in[1], zta = xi_in[2];
	double x, y, z;

	x = 5.0*xi;

	double delta_s = 2.0;
	double eta_stretch = 0.5*(1.0 + tanh(delta_s*(eta - 0.5))/tanh(delta_s*0.5));

	if(x >= 1.0 && x < 2.0){
	   double noz_x = x - 1;
	   double y_upper = 0.7 + 0.3*pow(cos(noz_x*M_PI), 4.0);
	   double y_lower = 0.3 - 0.3*pow(cos(noz_x*M_PI), 4.0);

	   y = y_lower + (y_upper - y_lower)*eta_stretch;

	}else{

	   y = eta_stretch;

	}

	z = 0.5*zta;

	x_out[0] = x;
	x_out[1] = y;
	x_out[2] = z;

}

void AlgebraicSingleBlockMesh::generateCylinderGrid(double xi_in[3], int cylRIndex, int cylRMax, int cylThetaIndex, int cylThetaMax, double x_out[3]){

	double xi = xi_in[0], eta = xi_in[1], zta = xi_in[2];
	double x, y, z;


	//radial distribution
	double first_off = 0.0008;
	double R = 0.5;

	double inner_growth_rate = 1.04;
	int Nlayers = 50;

	//sloppy programming 
	int RMAX = cylRMax;
	double r[RMAX];

	for(int ip = 0; ip < Nlayers; ip++){
	    if(ip == 0){
		r[ip] = R;
	    }else if(ip == 1){
		r[ip] = r[ip-1] + first_off;	
	    }else{
		r[ip] = r[ip-1] + inner_growth_rate*(r[ip-1] - r[ip-2]);
	    }

	}

	double outer_growth_rate = 1.02;
	for(int ip = Nlayers; ip < RMAX; ip++){
	    r[ip] = r[ip-1] + outer_growth_rate*(r[ip-1] - r[ip-2]);
	}

	//Theta distribution
	double thetabase = (double)cylThetaIndex/(double)cylThetaMax;
	double deltastretch = 2.0;
	double thetastretch = 0.5*(1.0 + tanh(deltastretch*(thetabase - 0.5))/tanh(deltastretch/2.0));


	double local_theta = 2.0*M_PI*thetastretch;

	//Get x and y out of the r & theta coordinates

	x = r[cylRIndex]*cos(local_theta);
	y = r[cylRIndex]*sin(local_theta);

	//{
	//    cout << "x, y, theta, r = " << x << " "  << y << " " << local_theta << " " << r[cylRIndex] << endl;
	//}

	//take z to be uniform distribution of length pi
	z = zta*M_PI;

	x_out[0] = x;
	x_out[1] = y;
	x_out[2] = z;


}



#endif
