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

    //  parallel.PleaseNeverFinalizeMPI();

    COUT << "Parallel grid size: ("<<parallel.grid_size()[0]<<","<<parallel.grid_size()[1]<<"). "<<endl;
    //-----------------------   end   ------------------------


    //------------   Declaration of a Lattice   --------------
    int dim = 3;
    int latSize[3] = {128,128,128};
    int halo = 3;
    Lattice lat(dim,latSize,halo);

    //-----------   Declaration of the Fields   --------------
    Field<Real> rho(lat,3);
    Site x(lat);
    Site xr(lat);


    for(x.first();x.test();x.next())
    {
	for(int j=0;j<3;j++)rho(x,j)=x.coord(j);
    }
    rho.updateHalo();

    forallboundary_start(lat,xb)
	for(int j=0;j<3;j++)rho(xb,j)=xb.coord(j);
    forallboundary_stop

    int error_count=0;

    for(int p=0;p<parallel.size();p++)
    {
	MPI_Barrier(MPI::COMM_WORLD);
	if(parallel.rank()==0)
	{
	    for(x.first();x.test();x.next())
	    {
		for(int k=0;k<3;k++)
		{
		    xr=x-k;
		    if((x.coord(k)-1)!=rho(xr,k))
		    {
			cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
			error_count++;
		    }
		    xr=x+k;
		    if((x.coord(k)+1)!=rho(xr,k))
		    {
			cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
			error_count++;
		    }
		    xr=x-k-k;
                    if((x.coord(k)-2)!=rho(xr,k))
                    {
                        cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
                        error_count++;
                    }
                    xr=x+k+k;
                    if((x.coord(k)+2)!=rho(xr,k))
                    {
                        cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
                        error_count++;
                    }
		    xr=x-k-k-k;
                    if((x.coord(k)-3)!=rho(xr,k))
                    {
                        cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
                        error_count++;
                    }
                    xr=x+k+k+k;
                    if((x.coord(k)+3)!=rho(xr,k))
                    {
                        cout<<"issues on site"<<x<<"direction "<< k <<" ; "<< xr <<"field: ("<<rho(xr,0)<<","<<rho(xr,1)<<","<<rho(xr,2)<<")"<<endl;
                        error_count++;
                    }
                    /*
		    xr=x-k-k;
		    if((x.coord(k)-2)==rho(xr,k))cout<<"issues on site"<<x<<", direction "<< k <<" : "<< xr <<endl;
		    xr=x+k+k;
		    if((x.coord(k)+2)==rho(xr,k))cout<<"issues on site"<<x<<", direction "<< k <<" : "<< xr <<endl;
		    xr=x-k-k-k;
		    if((x.coord(k)-3)==rho(xr,k))cout<<"issues on site"<<x<<", direction "<< k <<" : "<< xr <<endl;
		    xr=x+k+k+k;
		    if((x.coord(k)+3)==rho(xr,k))cout<<"issues on site"<<x<<", direction "<< k <<" : "<< xr <<endl;
		    */
		}
	    }
	}
	MPI_Barrier(MPI::COMM_WORLD);
    }
    cout<<"process: "<<parallel.rank()<<", numer of error: "<<error_count<<endl;
/*
    for(int i=0; i<lat.sitesLocalGross();i++)if(rho(i)!=0)
    {
      x.setIndex(i);
      if(x.coord(0)<lat.size(0) && x.coord(0)>=0 && x.coord(1)<lat.size(1) && x.coord(1)>=0 && x.coord(2)<lat.size(2) && x.coord(2)>=0)
      {
        //cout<<rho(i)<<endl;
	  cout<< "wrong boundary points:" << x.coord(0)<<" , "<< x.coord(1)<<" , "<< x.coord(2)<<endl;
      }
    }

    for(int i=0; i<lat.sitesLocalGross();i++)if(rho(i)!=1)
    {
      x.setIndex(i);
      if(x.coord(0)>=lat.size(0) || x.coord(0)<0 || x.coord(1)>=lat.size(1) || x.coord(1)<0 || x.coord(2)>=lat.size(2) || x.coord(2)<0)
      {
	  cout<< "boundary points not taken into acount:" << x.coord(0)<<" , "<< x.coord(1)<<" , "<< x.coord(2)<<endl;
      }
    }
*/


//    MPI_Finalize();





    //--------------------------------------------------------
}
