/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */

#include <scorep/SCOREP_User.h>

#include "LATfield2.hpp"
using namespace LATfield2;



int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m,latSize;

    SCOREP_USER_REGION_DEFINE(init)
    SCOREP_USER_REGION_DEFINE(loop)
    SCOREP_USER_REGION_DEFINE(der1)
    SCOREP_USER_REGION_DEFINE(der2)
    SCOREP_USER_REGION_DEFINE(der3)
    SCOREP_USER_REGION_DEFINE(der4)


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
        case 's':
  				latSize =  atoi(argv[++i]);
  				break;
		}
	}

    parallel.initialize(n,m);

    //  parallel.PleaseNeverFinalizeMPI();

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    int halo = 2;

  SCOREP_USER_REGION_BEGIN( init, "initialisation",
      SCOREP_USER_REGION_TYPE_COMMON )


    Lattice lat(dim,latSize,halo);
    Field<double> rho(lat,1);
    Field<double> phi(lat,3);
    Field<double> q(lat,3,3,LATfield2::symmetric);
    Field<double> K(lat,3,3,LATfield2::symmetric);

    SCOREP_USER_REGION_END( init )

    double dx=0.2;

    Site x(lat);
for(int l = 0;l<5;l++)
{

  SCOREP_USER_REGION_BEGIN( loop, "forallsites",
      SCOREP_USER_REGION_TYPE_LOOP )
    for(x.first();x.test();x.next())
    {
        SCOREP_USER_REGION_BEGIN( der1, "diRho",
            SCOREP_USER_REGION_TYPE_COMMON )
        for(int i=0;i<3;i++)phi(x,i) =rho(x+i)+rho(x-i)/(2.0*dx);
        SCOREP_USER_REGION_END( der1 )

        SCOREP_USER_REGION_BEGIN( der2, "didiRho",
            SCOREP_USER_REGION_TYPE_COMMON )
        for(int i=0;i<3;i++)q(x,i,i) =(rho(x+i) + rho(x-i) - 2.0 * rho(x))/(dx*dx);
        SCOREP_USER_REGION_END( der2 )

        SCOREP_USER_REGION_BEGIN( der3, "didjRhp",
            SCOREP_USER_REGION_TYPE_COMMON )
        for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)q(x,i,j) = (rho(x+i+j) - rho(x-i+j) - rho(x+i-j) + rho(x-i-j))/(4.0*dx*dx);
        SCOREP_USER_REGION_END( der3 )

        SCOREP_USER_REGION_BEGIN( der4, "qijKij",
            SCOREP_USER_REGION_TYPE_COMMON )
        for(int i=0;i<3;i++)
        {
          phi(x,i) = 0;
          for(int j=0;j<3;j++)rho(x,i) += q(x,i,j)*K(x,i,j);
        }
        SCOREP_USER_REGION_END( der4 )
    }
  SCOREP_USER_REGION_END( loop )
}
/*
    for(x.first();x.test();x.next())
    {
        for(int i=0;i<3;i++)q(x,i,i) = (rho(x+i) + rho(x-i) - 2.0 * rho(x))/(dx*dx);
        for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)q(x,i,j) = rho(x+i+j) - rho(x-i+j) - rho(x+i-j) + rho(x-i-j);
    }

    for(x.first();x.test();x.next())
    {
        rho(x) = 0;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)rho(x) += q(x,i,j)*K(x,i,j);
    }
*/

/*
    for(int i=0;i<32;i++)
    {
      a[i]=i;
      b[i]=0.1 + i;
      c[i]=0.2 + i;
      //COUT<<a[i]<<endl;;
    }

    Site y(lat);

    for(x.first();x.test();x.next())
    {
        for(int i=0;i<lat.size(0);i++ )
        {
          x.setCoord0(i);



          phi.value(x,0)=x.coord(0);
          phi.value(x,1)=x.coord(1)*2;
          phi.value(x,2)=x.coord(2)*3;

        }

    }

    for(x.first();x.test();x.next())
    {

        COUT<<x<<": "<<phi(x,0)[0];
        for(int i=1;i<32;i++)
        {
          COUT<<","<<phi(x,0)[i];
        }
        COUT<<endl;
        COUT<<x<<": "<<phi(x,1)[0];
        for(int i=1;i<32;i++)
        {
          COUT<<","<<phi(x,1)[i];
        }
        COUT<<endl;
        COUT<<x<<": "<<phi(x,2)[0];
        for(int i=1;i<32;i++)
        {
          COUT<<","<<phi(x,2)[i];
        }
        COUT<<endl;
    }

/*
    //-----------   Declaration of the Fields   --------------
    Field<Real> rho(lat,3);


    Site x(lat);


    Site xr(lat);


    for(x.first();x.test();x.next())
    {
	for(int j=0;j<3;j++)rho(x,j)=x.coord(j);
    }
    rho.updateHalo();
*/


    //--------------------------------------------------------
}
