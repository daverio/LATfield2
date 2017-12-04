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
    int latSize = 18;
    int halo = 2;
    Lattice lat(dim,latSize,halo);


    Field<Real> rho(lat,3);
    Field<Real> phi(lat);
    Field<Real> beta(lat,3);

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        //cout << " "<<ID<< endl;
        printf("MPI rank %d, openMP task %d \n", parallel.rank(),ID);
    }


    //-----------------------   end   ------------------------

    //--------------   Operations on Fields   ----------------
    Site x(lat);
    Site y(lat);


    for(x.first();x.test();x.next())
    {
      for(int i=0;i<3;i++)rho(x,i) = i;
    }


    for(x.first();x.test();x.next())
    {
      for(int i=0;i<3;i++)cout<<"comp "<< i <<" : " <<x<<rho(x,i)<<endl;

    }
    //if(x.setCoord(0,0,0))cout<<rho(x,0)<<endl;
    /*
    for(x.first();x.test();x.next())
    {
      for(int i=0;i<3;i++)cout<<"comp "<< i <<" : " <<x<<rho(x,i)<<endl;
    }
    */
    rho.updateHalo();
    /*

    for(x.first();x.test();x.nextValue())
    {

      y=x-0;
      if(y.coord(0)==-1){
        if(rho.value(y,0)!=latSize-1) cout<< "error down 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error down 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error down 0, comp 2" << endl;
      }

      y=x-0-0;
      if(y.coord(0)==-2){
        if(rho.value(y,0)!=latSize-2) cout<< "error down2 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error down2 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error down2 0, comp 2" << endl;
      }

      y=x+0;
      if(y.coord(0)==latSize){
        if(rho.value(y,0)!=0) cout<< "error up 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error up 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error up 0, comp 2" << endl;
      }

      y=x+0+0;
      if(y.coord(0)==latSize+1){
        if(rho.value(y,0)!=1) cout<< "error up2 0, comp 0" << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<< "error up2 0, comp 1" << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<< "error up2 0, comp 2" << endl;
      }


      y=x-1;
      if(y.coord(1)==-1){
        if(rho.value(y,1)!=latSize-1) cout<<y<< "error down 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error down 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }
      y=x-1-1;
      if(y.coord(1)==-2){
        if(rho.value(y,1)!=latSize-2) cout<<y<< "error down2 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down2 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error down2 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x+1;
      if(y.coord(1)==latSize){
        if(rho.value(y,1)!=0) cout<<y<< "error up 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error up 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }
      y=x+1+1;
      if(y.coord(1)==latSize+1){
        if(rho.value(y,1)!=1) cout<<y<< "error up2 1, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up2 1, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=y.coord(2)) cout<<y<< "error up2 1, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x-2;
      if(y.coord(2)==-1){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error down 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=latSize-1) cout<<y<< "error down 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x-2-2;
      if(y.coord(2)==-2){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error down2 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error down2 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=latSize-2) cout<<y<< "error down2 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x+2;
      if(y.coord(2)==latSize){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error up 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=0) cout<<y<< "error up 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }

      y=x+2+2;
      if(y.coord(2)==latSize+1){
        if(rho.value(y,0)!=y.coord(0)) cout<<y<< "error up2 2, comp 0: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2) << endl;
        if(rho.value(y,1)!=y.coord(1)) cout<<y<< "error up2 2, comp 1: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
        if(rho.value(y,2)!=1) cout<<y<< "error up2 2, comp 2: "<< rho.value(y,0) << "," << rho.value(y,1) << "," << rho.value(y,2)  << endl;
      }




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

    rho.saveHDF5("rho.h5","rho");
    phi.saveHDF5("phi.h5");


    rho.loadHDF5("rho.h5","rho");
    phi.loadHDF5("phi.h5");


    for(x.first();x.test();x.nextValue())
    {
      //for(int i=0;i<3;i++)if(rho.value(x,i) != x.coord(i))cout<<"error"<<endl;
      //if(phi.value(x)!=(rho.value(x+0,0)+rho.value(x+1,1)+rho.value(x+2,2)))cout<<"error"<<endl;
    }

    cout<<"done"<<endl;
    //--------------------------------------------------------
    */
}
