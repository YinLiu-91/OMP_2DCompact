#ifndef _MACROSH_
#define _MACROSH_

#define NUMTHREADSNEST 1

#define FOR_I3 for(int i = 0; i < 3; i++)
#define FOR_J3 for(int j = 0; j < 3; j++)
#define FOR_K3 for(int k = 0; k < 3; k++)

#define FOR_I2 for(int i = 0; i < 2; i++)
#define FOR_J2 for(int j = 0; j < 2; j++)

#define FOR_X for(int i = 0; i < Nx; i++)
#define FOR_Y for(int j = 0; j < Ny; j++)
#define FOR_XY for(int ip = 0; ip < Nx*Ny; ip++)

#define GET2DINDEX_XY j*Nx + i
#define GET2DINDEX_YX i*Ny + j

#define FOR_X0 for(int j = 0; j < Ny; j++){ int i = 0; int ip = GET2DINDEX_XY;
#define FOR_X1 for(int j = 0; j < Ny; j++){ int i = Nx-1; int ip = GET2DINDEX_XY;

#define FOR_Y0 for(int i = 0; i < Nx; i++){ int j = 0; int ip = GET2DINDEX_XY;
#define FOR_Y1 for(int i = 0; i < Nx; i++){ int j = Ny-1; int ip = GET2DINDEX_XY;

#define END_FORX0 }
#define END_FORX1 }
#define END_FORY0 }
#define END_FORY1 }

#define SIGNED_TET_VOLUME_6(A,B,C,D) (((B)[0]-(A)[0])*(((C)[1]-(A)[1])*((D)[2]-(A)[2])-((C)[2]-(A)[2])*((D)[1]-(A)[1])) + \
                                      ((B)[1]-(A)[1])*(((C)[2]-(A)[2])*((D)[0]-(A)[0])-((C)[0]-(A)[0])*((D)[2]-(A)[2])) + \
                                      ((B)[2]-(A)[2])*(((C)[0]-(A)[0])*((D)[1]-(A)[1])-((C)[1]-(A)[1])*((D)[0]-(A)[0])))


#endif
