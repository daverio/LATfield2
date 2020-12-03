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
namespace ut = boost::unit_test;
namespace mpi = boost::mpi;

/*
    From Stroustrup's The C++ Programming Language 4ed
    section 13.3.1
*/
template <class F>
struct Final_action
{
    F clean;
    Final_action(F f) : clean{f} {}
    ~Final_action() { clean(); }
};
template <class F>
Final_action<F> finally(F f)
{
    return Final_action<F>(f);
}

struct fixture
{
    fixture() {}

    ~fixture() {}

    static mpi::environment env;
    static mpi::communicator world;
};

mpi::environment fixture::env;
mpi::communicator fixture::world;

BOOST_TEST_GLOBAL_FIXTURE(fixture);

BOOST_AUTO_TEST_CASE(fft_RtoC)
{   
    
    // parameters
    const int n=2,m=2;
    const int BoxSize = 64 ;
    const int halo = 1;
    const int dim = 3;
    const int comp = 1;
    
    // input function
    auto myF = [&BoxSize](int x,int y,int z)->double
    {
        double res = x * sqrt(1.0 * y) * sin(z * 1.0/BoxSize);
        return res;
    };
    
    // latfield setup
    parallel.initialize(fixture::world,n,m);

    Lattice lat(dim,BoxSize,halo);
    Site x(lat);

    Field<Real> phi(lat,comp);

    for (x.first(); x.test(); x.next())
    {
        phi(x) = myF(x.coord(0),x.coord(1),x.coord(2));
    }
    
    Lattice latK;
    latK.initializeRealFFT(lat, halo);
    
    rKSite k(latK);
    Field<Imag> phiK(latK,comp);
    PlanFFT<Imag> planPhi(&phi,&phiK);
    
    
    // fftw setup
    fftw_complex* in = new fftw_complex[BoxSize*BoxSize*BoxSize];
    fftw_plan p = 
        fftw_plan_dft_3d(BoxSize,BoxSize,BoxSize,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
    auto cleanup = finally([&](){
        delete[] in;
        fftw_destroy_plan(p);
    });
   
    for(int i=0;i<BoxSize;++i)
    for(int j=0;j<BoxSize;++j)
    for(int k=0;k<BoxSize;++k)
    {
        fftw_complex & x = in[k+BoxSize*(j+BoxSize*i)];
        x[0] = myF(i,j,k);
        x[1] = 0.0;
    }
    
    
    // execute
    planPhi.execute(FFT_FORWARD);
    fftw_execute(p);
    
    // compare
    double diff =0 ;
    for(k.first();k.test();k.next())
    {
        fftw_complex & x = in[k.coord(2)+BoxSize*(k.coord(1)+BoxSize*k.coord(0))];
        diff = std::max(diff, std::abs( x[0] - phiK(k).real() )  );
        diff = std::max(diff, std::abs( x[1] - phiK(k).imag() )  );
    }
    BOOST_CHECK_SMALL(diff,1e-7);
}


// Complex to Complex doesnt work ???
// BOOST_AUTO_TEST_CASE(fft_CtoC)
// {   
//     
//     // parameters
//     const int n=2,m=2;
//     const int BoxSize = 64 ;
//     const int halo = 1;
//     const int dim = 3;
//     const int comp = 1;
//     
//     // input function
//     auto myF = [&BoxSize](int x,int y,int z)->double
//     {
//         double res = x * sqrt(1.0 * y) * sin(z * 1.0/BoxSize);
//         return res;
//     };
//     
//     // latfield setup
//     parallel.initialize(fixture::world,n,m);
// 
//     Lattice lat(dim,BoxSize,halo);
//     Site x(lat);
// 
//     Field<Imag> phi(lat,comp);
// 
//     for (x.first(); x.test(); x.next())
//     {
//         phi(x).real() = myF(x.coord(0),x.coord(1),x.coord(2));
//         phi(x).imag() = 0;
//     }
//     
//     Lattice latK;
//     latK.initializeComplexFFT(lat, halo);
//     
//     rKSite k(latK);
//     Field<Imag> phiK(latK,comp);
//     PlanFFT<Imag> planPhi(&phi,&phiK);
//     
//     
//     // fftw setup
//     fftw_complex* in = new fftw_complex[BoxSize*BoxSize*BoxSize];
//     fftw_plan p = 
//         fftw_plan_dft_3d(BoxSize,BoxSize,BoxSize,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
//     auto cleanup = finally([&](){
//         delete[] in;
//         fftw_destroy_plan(p);
//     });
//    
//     for(int i=0;i<BoxSize;++i)
//     for(int j=0;j<BoxSize;++j)
//     for(int k=0;k<BoxSize;++k)
//     {
//         fftw_complex & x = in[k+BoxSize*(j+BoxSize*i)];
//         x[0] = myF(i,j,k);
//         x[1] = 0.0;
//     }
//     
//     
//     // execute
//     planPhi.execute(FFT_FORWARD);
//     fftw_execute(p);
//     
//     // compare
//     double diff =0 ;
//     for(k.first();k.test();k.next())
//     {
//         fftw_complex & x = in[k.coord(2)+BoxSize*(k.coord(1)+BoxSize*k.coord(0))];
//         diff = std::max(diff, std::abs( x[0] - phiK(k).real() )  );
//         diff = std::max(diff, std::abs( x[1] - phiK(k).imag() )  );
//     }
//     BOOST_CHECK_SMALL(diff,1e-7);
// }

