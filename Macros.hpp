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

#define GET_Xp1 j*Nx + i + 1
#define GET_Xp2 j*Nx + i + 2
#define GET_Xp3 j*Nx + i + 3
#define GET_Xp4 j*Nx + i + 4
#define GET_Xp5 j*Nx + i + 5
#define GET_Xp6 j*Nx + i + 6
#define GET_Xp7 j*Nx + i + 7
#define GET_Xp8 j*Nx + i + 8
#define GET_Xp9 j*Nx + i + 9
#define GET_Xp10 j*Nx + i + 10

#define GET_Xm1 j*Nx + i - 1
#define GET_Xm2 j*Nx + i - 2
#define GET_Xm3 j*Nx + i - 3
#define GET_Xm4 j*Nx + i - 4
#define GET_Xm5 j*Nx + i - 5
#define GET_Xm6 j*Nx + i - 6
#define GET_Xm7 j*Nx + i - 7
#define GET_Xm8 j*Nx + i - 8
#define GET_Xm9 j*Nx + i - 9
#define GET_Xm10 j*Nx + i - 10

#define GET_Yp1 (j+1)*Nx + i
#define GET_Yp2 (j+2)*Nx + i
#define GET_Yp3 (j+3)*Nx + i
#define GET_Yp4 (j+4)*Nx + i
#define GET_Yp5 (j+5)*Nx + i
#define GET_Yp6 (j+6)*Nx + i
#define GET_Yp7 (j+7)*Nx + i
#define GET_Yp8 (j+8)*Nx + i
#define GET_Yp9 (j+9)*Nx + i
#define GET_Yp10 (j+10)*Nx + i

#define GET_Ym1 (j-1)*Nx + i
#define GET_Ym2 (j-2)*Nx + i
#define GET_Ym3 (j-3)*Nx + i
#define GET_Ym4 (j-4)*Nx + i
#define GET_Ym5 (j-5)*Nx + i
#define GET_Ym6 (j-6)*Nx + i
#define GET_Ym7 (j-7)*Nx + i
#define GET_Ym8 (j-8)*Nx + i
#define GET_Ym9 (j-9)*Nx + i
#define GET_Ym10 (j-10)*Nx + i

#define FILL_GET_Xp {GET_Xp1, GET_Xp2, GET_Xp3, GET_Xp4, GET_Xp5, GET_Xp6, GET_Xp7, GET_Xp8, GET_Xp9, GET_Xp10}
#define FILL_GET_Xm {GET_Xm1, GET_Xm2, GET_Xm3, GET_Xm4, GET_Xm5, GET_Xm6, GET_Xm7, GET_Xm8, GET_Xm9, GET_Xm10}

#define FILL_GET_Yp {GET_Yp1, GET_Yp2, GET_Yp3, GET_Yp4, GET_Yp5, GET_Yp6, GET_Yp7, GET_Yp8, GET_Yp9, GET_Yp10}
#define FILL_GET_Ym {GET_Ym1, GET_Ym2, GET_Ym3, GET_Ym4, GET_Ym5, GET_Ym6, GET_Ym7, GET_Ym8, GET_Ym9, GET_Ym10}

#define SIGNED_TET_VOLUME_6(A,B,C,D) (((B)[0]-(A)[0])*(((C)[1]-(A)[1])*((D)[2]-(A)[2])-((C)[2]-(A)[2])*((D)[1]-(A)[1])) + \
                                      ((B)[1]-(A)[1])*(((C)[2]-(A)[2])*((D)[0]-(A)[0])-((C)[0]-(A)[0])*((D)[2]-(A)[2])) + \
                                      ((B)[2]-(A)[2])*(((C)[0]-(A)[0])*((D)[1]-(A)[1])-((C)[1]-(A)[1])*((D)[0]-(A)[0])))


#endif
