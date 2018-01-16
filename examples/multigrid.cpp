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


    Lattice lat(3,128,1);
    Field<double> phi;
    COUT<<"phi is init (should not): " <<phi.IsInitialized()<<endl;
    phi.initialize(lat,3);
    COUT<<"phi is init (should): " <<phi.IsInitialized()<<endl;
    phi.alloc();


    MultiField<double> * mphi;

    MultiGrid mg_engine;
    mg_engine.initialize(&lat,18,8);

    mg_engine.intitialize_Field(&phi,mphi);

    Site x(lat);

    COUT<<"mphi[0] is init: " <<mphi[0].IsInitialized()<<endl;

    for(x.first();x.test();x.next())
    {
      mphi[0](x,0)=x.coord(0);
      mphi[0](x,1)=x.coord(1);
      mphi[0](x,2)=x.coord(2);
    }

    COUT<< "number of level: "<< mg_engine.nl()<<endl;
    COUT<< "number of parallel layer: "<< mg_engine.npl()<<endl;


    for(int i=0;i<mg_engine.nl()-1;i++)
    {
      //COUT<<"Saving level:"<<i <<endl;
      //mphi[i].saveHDF5("outputs/level_"+int2string(i,99)+".h5");
      COUT<<"restriction starting on level:"<<i <<endl;

      mg_engine.restrict(mphi,i);

      //Site x(mphi[i+1].lattice());
    }

    for(int i=0;i<mg_engine.nl();i++)
    {
      if(parallel.layer(mg_engine.player(i)).isPartLayer())
      {
        Site xl(mphi[i].lattice());
        long error[3]={0,0,0};

        for(xl.first();xl.test();xl.next())
        {
          if(mphi[i](xl,0)!=xl.coord(0)*pow(2,i) )error[0]++;
          if(mphi[i](xl,1)!=xl.coord(1)*pow(2,i) )error[1]++;
          if(mphi[i](xl,2)!=xl.coord(2)*pow(2,i) )error[2]++;
          //if(i==3)if( mphi[i](xl,2)!=xl.coord(2)*pow(2,i) ||
          //   mphi[i](xl,1)!=xl.coord(1)*pow(2,i) ||
          //   mphi[i](xl,0)!=xl.coord(0)*pow(2,i) )
          //   cout << xl << mphi[i](xl,0) << " "<< mphi[i](xl,1)<< " "<< mphi[i](xl,2)<<endl;

        }

        parallel.layer(mg_engine.player(i)).sum(error,3);

        if(parallel.layer(mg_engine.player(i)).isRoot())cout<<"error restiction level["<<i<<"] :"<< error[0] <<" "<< error[1] <<" "<< error[2] <<endl;
      }
    }


    for(int i=mg_engine.nl()-1;i>0;i--)
    {
      //COUT<<"Saving level:"<<i <<endl;
      //mphi[i].saveHDF5("outputs/level_"+int2string(i,99)+".h5");
      COUT<<"prologonation starting on level:"<<i <<endl;


      if(mg_engine.isPartLayer(i))
      {
        mphi[i].updateHalo();
        forallboundary_start(mphi[i].lattice(),xg)
          //cout<<i<<" "<<xg<<endl;
          for(int c = 0; c< mphi[i].components();c++)mphi[i](xg,c) = xg.coord(c)*pow(2,i);
        forallboundary_stop


      }


      mg_engine.prolong(mphi,i);


    }


    for(int i=0;i<mg_engine.nl();i++)
    {
      if(parallel.layer(mg_engine.player(i)).isPartLayer())
      {
        Site xl(mphi[i].lattice());
        long error[3]={0,0,0};

        for(xl.first();xl.test();xl.next())
        {
          if(mphi[i](xl,0)!=xl.coord(0)*pow(2,i) )error[0]++;
          if(mphi[i](xl,1)!=xl.coord(1)*pow(2,i) )error[1]++;
          if(mphi[i](xl,2)!=xl.coord(2)*pow(2,i) )error[2]++;
          //if(i==3)if( mphi[i](xl,2)!=xl.coord(2)*pow(2,i) ||
          //   mphi[i](xl,1)!=xl.coord(1)*pow(2,i) ||
          //   mphi[i](xl,0)!=xl.coord(0)*pow(2,i) )
          //   cout << xl << mphi[i](xl,0) << " "<< mphi[i](xl,1)<< " "<< mphi[i](xl,2)<<endl;

        }
        parallel.layer(mg_engine.player(i)).sum(error,3);
        if(parallel.layer(mg_engine.player(i)).isRoot())cout<<"error prologation level["<<i<<"] :"<< error[0] <<" "<< error[1] <<" "<< error[2] <<endl;
      }
    }

    /*
    long error[3]={0,0,0};

    for(x.first();x.test();x.next())
    {
      if(phi(x,0)!=x.coord(0))error[0]++;
      if(phi(x,1)!=x.coord(1))error[1]++;
      if(phi(x,2)!=x.coord(2))error[2]++;
    }
    parallel.sum(error[0]);
    parallel.sum(error[1]);
    parallel.sum(error[2]);


    COUT<< error[0] <<" "<< error[1] <<" "<< error[2] <<endl;

    */




    //--------------------------------------------------------
    return 0;
}
