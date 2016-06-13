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
    int latSize[3] = {512,512,512};
    int halo = 2;
    Lattice lat(dim,latSize,halo);
/*
    LFvector<double> a(&lat);
    LFvector<double> b(&lat);
    LFvector<double> c(&lat);

    for(int i=0;i<32;i++)
    {
      a[i]=i;
      b[i]=i;
      c[i]=2;
      //COUT<<a[i]<<endl;;
    }


    c += 2.0+(a*b);

    for(int i=0;i<32;i++)
    {
      COUT<<"a["<<i<<"]: "<<a[i]<<" , b["<<i<<"]: "<<b[i]<<" , c["<<i<<"]: "<<c[i]<<endl;;
    }
*/

    Field<double> rho(lat,1);
    Field<double> phi(lat,3);
    Field<double> q(lat,3,3,LATfield2::symmetric);
    Field<double> K(lat,3,3,LATfield2::symmetric);
    double dx=0.2;

    Site x(lat);

    for(x.first();x.test();x.next())
    {
        for(int i=0;i<3;i++)phi(x,i) = (rho(x+i)-rho(x-i))/(2.0*dx);
    }

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
