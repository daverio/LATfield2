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

    parallel.PleaseNeverFinalizeMPI();

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    int latSize[3] = {8,8,8};
    int halo = 2;
    Lattice lat(dim,latSize,halo);

    //-----------   Declaration of the Fields   --------------
    Field<Real> rho(lat);
    Site x(lat);

    for(x.first();x.test();x.next())rho(x)=0;
    rho.updateHalo();

    forallboundary_start(lat,xb)
    rho(xb)=1;
    forallboundary_stop

    for(int i=0; i<lat.sitesLocalGross();i++)if(rho(i)==1)
    {
      x.setIndex(i);
      if(x.coord(0)<lat.size(0) && x.coord(0)>=0 && x.coord(1)<lat.size(1) && x.coord(1)>=0 && x.coord(2)<lat.size(2) && x.coord(2)>=0)
      {
        //cout<<rho(i)<<endl;
       cout<< "wrong boundary points:" << x.coord(0)<<" , "<< x.coord(1)<<" , "<< x.coord(2)<<endl;
      }
    }

    for(int i=0; i<lat.sitesLocalGross();i++)if(rho(i)==0)
    {
      x.setIndex(i);
      if(x.coord(0)>=lat.size(0) || x.coord(0)<0 || x.coord(1)>=lat.size(1) || x.coord(1)<0 || x.coord(2)>=lat.size(2) || x.coord(2)<0)
      {
        cout<< "boundary points not taken into acount:" << x.coord(0)<<" , "<< x.coord(1)<<" , "<< x.coord(2)<<endl;
      }
    }



    MPI_Finalize();





    //--------------------------------------------------------
}
