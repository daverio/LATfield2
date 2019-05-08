/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */


#include "LATfield2.hpp"
using namespace LATfield2;




#define NTIMER 3
#define T_PROJ 0
#define T_FOP  1
#define T_FOPVECTOR  2

int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m;
    int latSize = 64;
    int loop = 10;
    int ratio = 2;

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
    		loop =  atoi(argv[++i]);
    		break;
      case 'r':
      	ratio =  atoi(argv[++i]);
      	break;
		}
	}

    parallel.initialize(n,m);

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    //int latSize = 256;
    int halo = 2;
    int vectorSize = latSize/ratio;
    Lattice lat(dim,latSize,halo,vectorSize);
    Lattice lat_part(dim,latSize,0,vectorSize);

    Site x(lat);
    long npts3d = latSize*latSize*latSize;

    Field<Real> rho(lat,3);
    Field<Real> phi(lat);
    Field<Real> beta(lat,3);

    MPI_timer timer(NTIMER);



    Real boxSize[3] = {(Real)latSize,(Real)latSize,(Real)latSize};
    part_simple_info particles_global_info;
    part_simple_dataType particles_dataType;

    particles_global_info.mass=8;
    particles_global_info.relativistic=false;
    set_parts_typename(&particles_global_info,"part_simple");

    Particles<part_simple,part_simple_info,part_simple_dataType> parts;
    parts.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);

    Site xp(lat_part);



    part_simple part;

    for(xp.first();xp.test();xp.next())
    {
        part.ID=0;
        part.pos[0]= ((Real)xp.coord(0)+0.5) * parts.res();
        part.pos[1]= ((Real)xp.coord(1)+0.5) * parts.res();
        part.pos[2]= ((Real)xp.coord(2)+0.5) * parts.res();
        part.vel[0]=1.0;
        part.vel[1]=1.0;
        part.vel[2]=1.0;
        //part.mass=0.22;
        parts.addParticle_global(part);
    }

    timer.start(T_PROJ);
    projection_init(&phi);
    scalarProjectionCIC_project(&parts,&phi);
    scalarProjectionCIC_comm(&phi);
    timer.stop(T_PROJ);

    for(x.first();x.test();x.next())
    {
      if(phi(x) != 8)cout<<" projection error: phi"<<x <<"= "<<phi(x)<<endl;
    }

    for(x.first();x.test();x.next())
    {
      for(int i=0;i<3;i++)rho(x,i) = x.coord(i);
    }

    rho.updateHalo();

    for(int i = 0;i<loop;i++)
      {


        timer.start(T_FOP);
        for(x.first();x.test();x.next())
        {
          phi(x)=rho(x,0)+rho(x,1)+rho(x,2);
          for(int c = 0;c<3;c++)beta(x,c)=phi(x+c)-phi(x-c);
        }

        for(x.first();x.test();x.next())
        {
            for(int c = 0;c<3;c++)beta(x,c)=phi(x+c)+phi(x-c);
        }
        timer.stop(T_FOP);


        timer.start(T_FOPVECTOR);
        for(x.first();x.test();x.nextVector())
        {
          //COUT<<"+++++++++++++++++++++"<<endl;
          phi.vect(x)=rho.vect(x,0)+rho.vect(x,1)+rho.vect(x,2);
          for(int c = 0;c<3;c++)beta.vect(x,c)=phi.vect(x+c)-phi.vect(x-c);
        }

        for(x.first();x.test();x.nextVector())
        {
          for(int c = 0;c<3;c++)
          beta.vect(x,c)=phi.vect(x+c)+phi.vect(x-c);
        }
        timer.stop(T_FOPVECTOR);

      }


    for(x.first();x.test();x.next())
    {
       for(int i=0;i<3;i++)if(rho(x,i) != x.coord(i))cout<<"error rho"<<endl;
       if(phi(x)!=(rho(x,0)+rho(x,1)+rho(x,2)))cout<<"error phi"<<endl;
    }


    COUT<<"projection time: "<<timer.aveTimer(T_PROJ)<<endl;
    COUT<<"field op time: "<<timer.aveTimer(T_FOP)<<endl;
    COUT<<"field op vector time: "<<timer.aveTimer(T_FOPVECTOR)<<endl;

    //--------------------------------------------------------

}
