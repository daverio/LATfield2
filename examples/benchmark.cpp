/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */


#include "LATfield2.hpp"
using namespace LATfield2;


int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m;

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
		}
	}

    parallel.initialize(n,m);

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    int latSize = 128;
    int halo = 2;
    Lattice lat(dim,latSize,halo);

    Field<Real> phi(lat);
    Field<Real> rho(lat,3);
    Field<Real> lap(lat);
    Field<Real> mcomp(lat,10);

    LFvector<Real> norm(lat);

    double dx = 1.0/(double)latSize;
    //-----------------------   end   ------------------------

    //--------------   Operations on Fields   ----------------
    Site x(lat);
    Site y(lat);


    //field implementation

    //1 comp
    for(x.first();x.test();x.nextValue())
    {
      phi.value(x) = sin( 2.0 * M_PI * ( (double)x.coord(2) - (double)x.coord(1) ) * dx );
    }

    // 10 comp
    for(x.first();x.test();x.nextValue())
    {
      for(int i=0;i<10;i++)mcomp.value(x,i) = cos(2.0 * M_PI * x.coord(0)*dx);
    }

    //norm of large vector
    for(x.first();x.test();x.next())
    {
      norm = mcomp(x,0)* mcomp(x,0);
      for(int i=1;i<10;i++)norm += mcomp(x,i)* mcomp(x,i);
      norm = vsqrt(norm);
    }

    //3 point stencil derivatives
    for(x.first();x.test();x.next())
    {
      for(int i = 0; i<3;i++)rho(x,i) = (phi(x+i)-phi(x-i))/(2.0 * dx);
    }

    //5 point stencil derivatives
    for(x.first();x.test();x.next())
    {
      for(int i = 0; i<3;i++)rho(x,i) = (phi(x-i-i) - phi(x-i)*8.0 + phi(x+i)*8.0 - phi(x+i+i))/(12.0*dx);
    }

    //3 point stencil Laplacian
    for(x.first();x.test();x.next())
    {
      lap(x) = (phi(x+0)+phi(x-0)+phi(x+1)+phi(x-1)+phi(x+2)+phi(x-2)-6.0*phi(x))/(dx*dx);
    }

    //5 point stencil Laplacian
    for(x.first();x.test();x.next())
    {
      lap(x) = (phi(x+0)*16.0 -phi(x+0+0) + phi(x-0)*16.0 - phi(x-0-0)
                    + phi(x+1)*16.0 -phi(x+1+1) + phi(x-1)*16.0 - phi(x-1-1)
                    + phi(x+2)*16.0 -phi(x+2+0) + phi(x-2)*16.0 - phi(x-2-2)
                    - phi(x)*90) /(dx*dx*12.0);
    }



    //updateHalo

    //1 comp
    phi.updateHalo();

    //3 comp
    rho.updateHalo();

    //10 comp
    mcomp.updateHalo();











}
