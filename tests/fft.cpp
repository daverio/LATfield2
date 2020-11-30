#define BOOST_TEST_MODULE Test FFT
#define FFT3D
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/unit_test.hpp>
#include <random>

#include <fftw3.h>
#include "LATfield2.hpp"

using namespace LATfield2;

BOOST_AUTO_TEST_CASE(simple_fft)
{
    std::default_random_engine rng;
    std::uniform_real_distribution<double> U(0.0,1.0);

    const int n=2,m=2;
    const int BoxSize = 64 ;
    const int halo = 1;
    const int dim = 3;
    const int comp = 1;

    parallel.initialize(n,m);

    Lattice lat(dim,BoxSize,halo);
    Site x(lat);

    Field<Real> phi(lat,comp);
    phi.alloc();

    for (x.first(); x.test(); x.next())
    {
        phi(x) = 1.0;
    }
    
    //Lattice latK(lat,halo);
    //rKSite k(latK);
    //Field<Imag> phiK;
    //phiK.initialize(latK,comp);
    //PlanFFT<Imag> planPhi(&phi,&phiK);

//  parallel.initialize(n,m);
//  //parallel.PleaseNeverFinalizeMPI();
//  
//  Lattice lat;
//  lat.initialize(dim,BoxSize,halo);
//  Site x(lat);
//  Field<Real> phi;
//  phi.initialize(lat,comp);
//  
//  for(x.first();x.test();x.next())
//  {
//      phi(x) = U(rng);
//  }

//    Lattice latKreal;
//    latKreal.initializeRealFFT(lat, halo);
//    rKSite kReal(latKreal);
//    Field<Imag> phiK;
//    phiK.initialize(latKreal,comp);
//
//    PlanFFT<Imag> planReal(&phi,&phiK);
//    
//    for(kReal.first();kReal.test();kReal.next())
//    {
//      phiK(kReal).real()=0.0;
//      phiK(kReal).imag()=0.0;
//    }

//   planReal.execute(FFT_FORWARD);
}

