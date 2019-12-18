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
