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

    Lattice lat(3,512,2);
    Field<double> phi;
    COUT<<"phi is init (should not): " <<phi.IsInitialized()<<endl;
    phi.initialize(lat,3);
    COUT<<"phi is init (should): " <<phi.IsInitialized()<<endl;
    phi.alloc();




    MultiField<double> * mphi;

    MultiGrid mg_engine;

    mg_engine.initialize(&lat,18,4);
    mg_engine.intitialize_Field(&phi,mphi);

    Site x(lat);

    COUT<<"mphi[0] is init: " <<mphi[0].IsInitialized()<<endl;

    for(x.first();x.test();x.next())
    {
      mphi[0](x,0)=1;
      mphi[0](x,1)=1;
      mphi[0](x,2)=1;
    }

    int error[3]={0,0,0};
    for(x.first();x.test();x.next())
    {
      if(phi(x,0)!=1)error[0]++;
      if(phi(x,1)!=1)error[1]++;
      if(phi(x,2)!=1)error[2]++;
    }
    parallel.sum(error[0]);
    parallel.sum(error[1]);
    parallel.sum(error[2]);


    COUT<< error[0] <<" "<< error[1] <<" "<< error[2] <<endl;






    //--------------------------------------------------------
}
