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
    int latSize = 512;
    int halo = 1;
    Lattice lat(dim,latSize,halo);


    Field<Real> rho(lat,3);
    Field<Real> phi(lat);
    Field<Real> beta(lat,3);

    //-----------------------   end   ------------------------

    //--------------   Operations on Fields   ----------------
    Site x(lat);


    for(x.first();x.test();x.nextValue())
    {
      for(int i=0;i<3;i++)rho.value(x,i) = x.coord(i);
    }

    rho.updateHalo();

    for(x.first();x.test();x.nextValue())
    {

      if(rho.value(x+0,1) != x.coord(1)) cout<<"error"<<endl;
      if(rho.value(x+0,2) != x.coord(2)) cout<<"error"<<endl;
      if(rho.value(x-0,1) != x.coord(1)) cout<<"error"<<endl;
      if(rho.value(x-0,2) != x.coord(2)) cout<<"error"<<endl;

    }

    rho.updateHalo();


    for(int i = 0;i<1;i++)
      {
    for(x.first();x.test();x.next())
    {
      phi(x)=rho(x,0)+rho(x,1)+rho(x,2);
      for(int i = 0;i<3;i++)beta(x,i)=phi(x+i)-phi(x-i);
    }
      }


    for(x.first();x.test();x.nextValue())
      {

	if(rho.value(x+0,1) != x.coord(1)) cout<<"error"<<endl;
	if(rho.value(x+0,2) != x.coord(2)) cout<<"error"<<endl;
	if(rho.value(x-0,1) != x.coord(1)) cout<<"error"<<endl;
	if(rho.value(x-0,2) != x.coord(2)) cout<<"error"<<endl;

      }

    //    rho.saveHDF5("rho.h5","rho");
    //phi.saveHDF5("phi.h5");


    //rho.loadHDF5("rho.h5","rho");
    //phi.loadHDF5("phi.h5");


    for(x.first();x.test();x.nextValue())
    {
      //for(int i=0;i<3;i++)if(rho.value(x,i) != x.coord(i))cout<<"error"<<endl;
      //if(phi.value(x)!=(rho.value(x+0,0)+rho.value(x+1,1)+rho.value(x+2,2)))cout<<"error"<<endl;
    }

    cout<<"done"<<endl;
    //--------------------------------------------------------
}
