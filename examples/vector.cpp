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
    int latSize = 16;
    int halo = 1;
    Lattice lat(dim,latSize,halo, latSize/4);


    Field<Real> rho(lat,3);

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


    for(int p = 0; p<parallel.size();p++)
    {
      MPI_Barrier(parallel.lat_world_comm());
      if(p==parallel.rank())
      {
        cout<<"======================"<<endl;
        cout<<"process: "<<p<<endl;

        for(x.first();x.test();x.next())
        {
          //cout<< x <<endl;
        }

      }
      MPI_Barrier(parallel.lat_world_comm());
    }


    //--------------------------------------------------------
}
