/*! file poissonSolver.cpp
    Created by David Daverio.

    A simple example of LATfield2d usage. This exec solve the poisson equation in fourier space.


 */



#include <iostream>
#include "LATfield2.hpp"
#include <mpi.h>

using namespace LATfield2;



int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);
    int n,m;
    int BoxSize = 64;
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

	parallel.initialize(MPI_COMM_WORLD,n,m);


    Lattice lat;
    lat.initialize(3,BoxSize,halo);


    //Real to complex fourier transform

    Lattice latKreal;
    latKreal.initializeRealFFT(lat, khalo);

    Site x(lat);
    rKSite kReal(latKreal);

    Field<Real> phi;
    phi.initialize(lat,comp);

    Field<Imag> phiK;
    phiK.initialize(latKreal,comp);

    PlanFFT<Imag> planReal(&phi,&phiK);


    for(kReal.first();kReal.test();kReal.next())
    {
      phiK(kReal).real()=0.0;
      phiK(kReal).imag()=0.0;
    }

    planReal.execute(FFT_FORWARD);
    planReal.execute(FFT_BACKWARD);

    
    MPI_Finalize();

}
