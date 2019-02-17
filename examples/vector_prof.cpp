/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */


#include "LATfield2.hpp"
using namespace LATfield2;

#include <scorep/SCOREP_User.h>

#include "timer.hpp"



#define NTIMER 16
#define T_FADD_COMP 0
#define T_FDER3     1
#define T_FADD_COMP_OLD 2
#define T_FDER3_OLD     3

int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m;
    int latSize = 64;
    int iteration = 1;

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
      case 'i':
        iteration = atoi(argv[++i]);
  			break;
		}
	}

    parallel.initialize(n,m);

    SCOREP_USER_REGION_DEFINE(fadd_comp)
    SCOREP_USER_REGION_DEFINE(fder3)

    int dim = 3;
    //int latSize = 256;
    int halo = 2;
    int vectorSize = 64;//latSize / 2;
    Lattice lat(dim,latSize,halo);
    Lattice lat_part(dim,latSize,0);

    Site x(lat);
    long npts3d = latSize*latSize*latSize;

    Field<Real> rho(lat,3);
    Field<Real> phi(lat);
    Field<Real> beta(lat,3);

    MPI_timer timer(NTIMER);


    for(x.first();x.test();x.next())
    {
      for(int i=0;i<3;i++)rho(x,i) = x.coord(i);
    }

    rho.updateHalo();

    for(int i = 0;i<iteration;i++)
      {
        timer.start(T_FADD_COMP);
        SCOREP_USER_REGION_BEGIN( fadd_comp, "fadd_comp", SCOREP_USER_REGION_TYPE_COMMON )
        for(x.first();x.test();x.next())
        {
          phi(x)=rho(x,0)+rho(x,1)+rho(x,2);
        }
        SCOREP_USER_REGION_END( fadd_comp )
        timer.stop(T_FADD_COMP);

        timer.start(T_FDER3);
        SCOREP_USER_REGION_BEGIN( fder3, "fadd_comp", SCOREP_USER_REGION_TYPE_COMMON )
        for(x.first();x.test();x.next())
        {
          for(int i = 0;i<3;i++)beta(x,i)=phi(x+i)-phi(x-i);
        }
        SCOREP_USER_REGION_END( fder3 )
        timer.stop(T_FDER3);
      }


    for(int i = 0;i<iteration;i++)
      {
        timer.start(T_FADD_COMP_OLD);
        for(x.first();x.test();x.next())
        {
          phi(x)=rho(x,0)+rho(x,1)+rho(x,2);
        }
        timer.stop(T_FADD_COMP_OLD);

        timer.start(T_FDER3_OLD);
        for(x.first();x.test();x.next())
        {
          for(int i = 0;i<3;i++)beta(x,i)=phi(x+i)-phi(x-i);
        }
        timer.stop(T_FDER3_OLD);
      }






    COUT<<"field comp sum time: "<<timer.timer(T_FADD_COMP)<<endl;
    COUT<<"field deriv_3p time: "<<timer.timer(T_FDER3)<<endl;
    COUT<<"field comp sum old time: "<<timer.timer(T_FADD_COMP_OLD)<<endl;
    COUT<<"field deriv_3p old time: "<<timer.timer(T_FDER3_OLD)<<endl;
    cout<<"done"<<endl;
    //--------------------------------------------------------

}
