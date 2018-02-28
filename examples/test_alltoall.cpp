
#include <iostream>
#include "LATfield2.hpp"

using namespace LATfield2;



int main(int argc, char **argv)
{
    int n,m;
    int BoxSize = 192;
    int halo = 1;
    int khalo =0;
    int dim = 3;
    int comp = 1;
    double sigma2=0.5;
    double res =0.5;


    for (int i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'n':
				n = atoi(argv[++i]);
				break;
			case 'm':
				m =  atoi(argv[++i]);
				break;
            case 'b':
                BoxSize = atoi(argv[++i]);
                break;
		}
	}
  //parallel.setNodeGeometry(3,4);
	parallel.initialize(n,m);

  Lattice lat;
  lat.initialize(3,BoxSize,halo);

  Lattice latKreal;
  latKreal.initializeRealFFT(lat, khalo);

  Site x(lat);

  Field<Real> phi;
  phi.initialize(lat,comp);

  Field<Imag> phiK;
  phiK.initialize(latKreal,comp);

  PlanFFT<Imag> planReal(&phi,&phiK);

  planReal.execute(FFT_FORWARD);
  planReal.execute(FFT_BACKWARD);

}
