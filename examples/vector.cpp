/*! file gettingStarted.cpp
    Created by David Daverio.

    A simple example of LATfield2 usage.
 */


#include "LATfield2.hpp"
using namespace LATfield2;

#include "timer.hpp"



#define NTIMER 16
#define T_PROJ 0
#define T_FOP  1

int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m;
    int latSize = 64;

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

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    //int latSize = 256;
    int halo = 2;
    int vectorSize = latSize / 2;
    Lattice lat(dim,latSize,halo,vectorSize);
    Lattice lat_part(dim,latSize,0,vectorSize);

    Site x(lat);
    long npts3d = latSize*latSize*latSize;

    Field<Real> rho(lat,3);
    Field<Real> phi(lat);
    Field<Real> beta(lat,3);

    MPI_timer timer(NTIMER);

    /*

    double boxSize[3] = {latSize,latSize,latSize};
    part_simple_info particles_global_info;
    part_simple_dataType particles_dataType;

    particles_global_info.mass=8;
    particles_global_info.relativistic=false;
    set_parts_typename(&particles_global_info,"part_simple");

    Particles<part_simple,part_simple_info,part_simple_dataType> parts;
    parts.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);

    Site xp(lat_part);



    part_simple part;

    for(xp.first();xp.test();xp.nextValue())
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

    for(x.first();x.test();x.nextValue())
    {
      if(phi.value(x) != 8)cout<<" projection error: phi"<<x <<"= "<<phi.value(x)<<endl;
    }

    //phi.saveHDF5("phi.h5");


    //-----------------------   end   ------------------------

    //--------------   Operations on Fields   ----------------


    */

    for(x.first();x.test();x.nextValue())
    {
      for(int i=0;i<3;i++)rho.value(x,i) = x.coord(i);
    }

    //for(x.first();x.test();x.next())
    //{
    //  phi(x) = 1;
    //}

    //Real sum = 0;

    //for(x.first();x.test();x.next())
    //{
    //  sum += phi(x);
    //}

    //parallel.sum(sum);

    //COUT<<"npts3d: "<<npts3d<< " ; sum: "<<sum<<endl;




    rho.updateHalo();

/*
    for(x.first();x.test();x.nextValue())
    {

      y=x.move(0,-1);
      if(y.coord(0)==-1){
        if(rho.value(y,0)!=latSize-1) cout<< "error down 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error down 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error down 0, comp 2" << endl;
      }

      y=x.move(0,-2);
      if(y.coord(0)==-2){
        if(rho.value(y,0)!=latSize-2) cout<< "error down2 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error down2 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error down2 0, comp 2" << endl;
      }

      y=x.move(0,1);
      if(y.coord(0)==latSize){
        if(rho.value(y,0)!=0) cout<< "error up 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error up 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error up 0, comp 2" << endl;
      }

      y=x.move(0,2);
      if(y.coord(0)==latSize+1){
        if(rho.value(y,0)!=1) cout<< "error up2 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error up2 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error up2 0, comp 2" << endl;
      }


      y=x.move(1,-1);
      if(y.coord(1)==-1){
        if(rho.value(y,1)!=latSize-1) cout<<y<< "error down 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error down 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }
      y=x.move(1,-2);
      if(y.coord(1)==-2){
        if(rho.value(y,1)!=latSize-2) cout<<y<< "error down2 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down2 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error down2 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x.move(1,1);
      if(y.coord(1)==latSize){
        if(rho.value(y,1)!=0) cout<<y<< "error up 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error up 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }
      y=x.move(1,2);
      if(y.coord(1)==latSize+1){
        if(rho.value(y,1)!=1) cout<<y<< "error up2 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up2 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error up2 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x.move(2,-1);
      if(y.coord(2)==-1){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error down 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=latSize-1) cout<<y<< "error down 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x.move(2,-2);
      if(y.coord(2)==-2){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down2 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error down2 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=latSize-2) cout<<y<< "error down2 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x.move(2,1);
      if(y.coord(2)==latSize){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error up 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=0) cout<<y<< "error up 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x.move(2,2);
      if(y.coord(2)==latSize+1){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up2 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error up2 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=1) cout<<y<< "error up2 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }




    }
*/
    rho.updateHalo();


    for(int i = 0;i<1;i++)
      {
        timer.start(T_FOP);
        for(x.first();x.test();x.next())
        {
          phi(x)=rho(x+0,0)+rho(x+1,1)+rho(x+2,2);
          for(int i = 0;i<3;i++)beta(x,i)=phi(x+i)-phi(x-i);
        }
        timer.stop(T_FOP);
      }




    //rho.saveHDF5("rho.h5","rho");
    //phi.saveHDF5("phi.h5");


    //rho.loadHDF5("rho.h5","rho");
    //phi.loadHDF5("phi.h5");

/*
    for(x.first();x.test();x.nextValue())
    {
       for(int i=0;i<3;i++)if(rho.value(x,i) != x.coord(i))cout<<"error"<<endl;
       if(phi.value(x)!=(rho.value(x+0,0)+rho.value(x+1,1)+rho.value(x+2,2)))cout<<"error"<<endl;
    }
*/

    COUT<<"projection time: "<<timer.timer(T_PROJ)<<endl;
    COUT<<"field op time: "<<timer.timer(T_FOP)<<endl;

    cout<<"done"<<endl;
    //--------------------------------------------------------

}
