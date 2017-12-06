/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */


#include "LATfield2.hpp"
using namespace LATfield2;

#include <scorep/SCOREP_User.h>

int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m;
    int ompTasks = 1;
    int pass;
    int boxsize;
    string filename;

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
         case 't':
    				 ompTasks =  atoi(argv[++i]);
    				 break;
         case 'p':
  				   pass =  atoi(argv[++i]);
  				   break;
         case 'f':
             filename = argv[++i];
             break;
         case 'b':
      		   boxsize =  atoi(argv[++i]);
  				   break;
		}
	}

    parallel.initialize(n,m);

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------



    SCOREP_USER_REGION_DEFINE(impl_1comp)
    SCOREP_USER_REGION_DEFINE(impl_10comp)
    SCOREP_USER_REGION_DEFINE(largenorm)
    SCOREP_USER_REGION_DEFINE(deriv_3p)
    SCOREP_USER_REGION_DEFINE(deriv_5p)
    SCOREP_USER_REGION_DEFINE(laplace_3p)
    SCOREP_USER_REGION_DEFINE(laplace_5p)
    SCOREP_USER_REGION_DEFINE(udHalo_1comp)
    SCOREP_USER_REGION_DEFINE(udHalo_3comp)
    SCOREP_USER_REGION_DEFINE(udHalo_10comp)


    fstream file;

    double timer,timer_ref,timer_start, timer_all;

    double timer_impl_1comp;
    double timer_impl_10comp;
    double timer_largenorm;
    double timer_deriv_3p;
    double timer_deriv_5p;
    double timer_laplace_3p;
    double timer_laplace_5p;
    double timer_udHalo_1comp;
    double timer_udHalo_3comp;
    double timer_udHalo_10comp;

timer_impl_1comp = 0.0;
timer_impl_10comp = 0.0;
timer_largenorm = 0.0;
timer_deriv_3p = 0.0;
timer_deriv_5p = 0.0;
timer_laplace_3p = 0.0;
timer_laplace_5p = 0.0;
timer_udHalo_1comp = 0.0;
timer_udHalo_3comp = 0.0;
timer_udHalo_10comp = 0.0;


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    int halo = 2;
    Lattice lat(dim,boxsize,halo);

    Field<Real> phi(lat);
    Field<Real> rho(lat,3);
    Field<Real> lap(lat);
    Field<Real> mcomp(lat,10);

    double norm;

    double dx = 1.0/(double)boxsize;
    //-----------------------   end   ------------------------

    //--------------   Operations on Fields   ----------------
    Site x(lat);
    Site y(lat);


    timer_start = MPI_Wtime();

    for(int p=0;p<pass;p++)
    {
      //field implementation

      //1 comp
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( impl_1comp, "impl_1comp",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        phi(x) = cos( 2.0 * M_PI * (double)x.coord(0) * dx );
      }
      SCOREP_USER_REGION_END( impl_1comp )
      timer = MPI_Wtime() - timer_ref;
      timer_impl_1comp += timer;
      // 10 comp
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( impl_10comp, "impl_10comp",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        for(int i=0;i<10;i++)mcomp(x,i) = (double)i*cos(2.0 * M_PI * (double)x.coord(0)*dx);
      }
      SCOREP_USER_REGION_END( impl_10comp )
      timer = MPI_Wtime() - timer_ref;
      timer_impl_10comp += timer;

      //norm of large vector
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( largenorm, "largenorm",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        norm = mcomp(x,0)* mcomp(x,0);
        for(int i=1;i<10;i++)norm += mcomp(x,i)* mcomp(x,i);
        norm = sqrt(norm);
      }
      SCOREP_USER_REGION_END( largenorm )
      timer = MPI_Wtime() - timer_ref;
      timer_largenorm += timer;

      //3 point stencil derivatives
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( deriv_3p, "deriv_3p",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        for(int i = 0; i<3;i++)rho(x,i) = (phi(x+i)-phi(x-i))/(2.0 * dx);
      }
      SCOREP_USER_REGION_END( deriv_3p )
      timer = MPI_Wtime() - timer_ref;
      timer_deriv_3p += timer;

      //5 point stencil derivatives
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( deriv_5p, "deriv_5p",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        for(int i = 0; i<3;i++)rho(x,i) = (phi(x-i-i) - phi(x-i)*8.0 + phi(x+i)*8.0 - phi(x+i+i))/(12.0*dx);
      }
      SCOREP_USER_REGION_END( deriv_5p )
      timer = MPI_Wtime() - timer_ref;
      timer_deriv_5p += timer;

      //3 point stencil Laplacian
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( laplace_3p, "laplace_3p",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        lap(x) = (phi(x+0)+phi(x-0)+phi(x+1)+phi(x-1)+phi(x+2)+phi(x-2)-6.0*phi(x))/(dx*dx);
      }
      SCOREP_USER_REGION_END( laplace_3p )
      timer = MPI_Wtime() - timer_ref;
      timer_laplace_3p += timer;

      //5 point stencil Laplacian
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( laplace_5p, "laplace_5p",
      SCOREP_USER_REGION_TYPE_COMMON )
      for(x.first();x.test();x.next())
      {
        lap(x) = (phi(x+0)*16.0 -phi(x+0+0) + phi(x-0)*16.0 - phi(x-0-0)
                      + phi(x+1)*16.0 -phi(x+1+1) + phi(x-1)*16.0 - phi(x-1-1)
                      + phi(x+2)*16.0 -phi(x+2+2) + phi(x-2)*16.0 - phi(x-2-2)
                      - phi(x)*90) /(dx*dx*12.0);
      }
      SCOREP_USER_REGION_END( laplace_5p )
      timer = MPI_Wtime() - timer_ref;
      timer_laplace_5p += timer;

      //updateHalo

      //1 comp
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( udHalo_1comp, "udHalo_1comp",
      SCOREP_USER_REGION_TYPE_COMMON )
      phi.updateHalo();
      SCOREP_USER_REGION_END( udHalo_1comp )
      timer = MPI_Wtime() - timer_ref;
      timer_udHalo_1comp += timer;

      //3 comp
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( udHalo_3comp, "udHalo_3comp",
      SCOREP_USER_REGION_TYPE_COMMON )
      rho.updateHalo();
      SCOREP_USER_REGION_END( udHalo_3comp )
      timer = MPI_Wtime() - timer_ref;
      timer_udHalo_3comp += timer;

      //10 comp
      timer_ref =MPI_Wtime();
      SCOREP_USER_REGION_BEGIN( udHalo_10comp, "udHalo_10comp",
      SCOREP_USER_REGION_TYPE_COMMON )
      mcomp.updateHalo();
      SCOREP_USER_REGION_END( udHalo_10comp )
      timer = MPI_Wtime() - timer_ref;
      timer_udHalo_10comp += timer;

    }

    timer_all = MPI_Wtime() - timer_start;


    parallel.max(timer_impl_1comp);
    parallel.max(timer_impl_10comp);
    parallel.max(timer_largenorm);
    parallel.max(timer_deriv_3p);
    parallel.max(timer_deriv_5p);
    parallel.max(timer_laplace_3p);
    parallel.max(timer_laplace_5p);
    parallel.max(timer_udHalo_1comp);
    parallel.max(timer_udHalo_3comp);
    parallel.max(timer_udHalo_10comp);
    parallel.max(timer_all);


    if(parallel.isRoot())
    {
      file.open( filename.c_str() , std::fstream::out | std::fstream::app); //does not need to trunk as only enter here is the file does not exist
      if(!file.is_open())
      {
        cout << "cannot open file: " << filename << ", exiting" << endl;
        exit(2);
      }
      else
      {
        file << n << "," << m << "," << ompTasks << "," << boxsize << "," << pass;
        file << "," << timer_impl_1comp;
        file << "," << timer_impl_10comp;
        file << "," << timer_largenorm;
        file << "," << timer_deriv_3p;
        file << "," << timer_deriv_5p;
        file << "," << timer_laplace_3p;
        file << "," << timer_laplace_5p;
        file << "," << timer_udHalo_1comp;
        file << "," << timer_udHalo_3comp;
        file << "," << timer_udHalo_10comp;
        file << endl;
      }
      file.close();
    }





}
