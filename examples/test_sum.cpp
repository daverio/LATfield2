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

    int a_int  = 1;
    long a_long = 1;
    float a_float = 1;
    double a_double = 1;

    parallel.sum(a_int);
    parallel.sum(a_double);
    parallel.sum(a_float);
    parallel.sum(a_long);
    cout<<parallel.rank()<<" , parallel.sum (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.sum (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.sum (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.sum (a_long) : "<<a_long<<endl;
    cout<<flush;
    parallel.barrier();
    COUT<<"============================================="<<endl;
    COUT<<"============================================="<<endl;


    a_int  = 1;
    a_long = 1;
    a_float = 1;
    a_double = 1;
    parallel.sum_to(a_int);
    parallel.sum_to(a_double);
    parallel.sum_to(a_float);
    parallel.sum_to(a_long);
    cout<<parallel.rank()<<" , parallel.sum_to 0 (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 0 (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 0 (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 0 (a_long) : "<<a_long<<endl;
    cout<<flush;
    parallel.barrier();
    COUT<<"============================================="<<endl;
    COUT<<"============================================="<<endl;


    a_int  = 1;
    a_long = 1;
    a_float = 1;
    a_double = 1;
    parallel.sum_to(a_int,1);
    parallel.sum_to(a_double,1);
    parallel.sum_to(a_float,1);
    parallel.sum_to(a_long,1);
    cout<<parallel.rank()<<" , parallel.sum_to 1 (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 1 (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 1 (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.sum_to 1 (a_long) : "<<a_long<<endl;
    cout<<flush;
    parallel.barrier();
    COUT<<"============================================="<<endl;
    COUT<<"============================================="<<endl;

    a_int  = parallel.rank();
    a_long = parallel.rank();
    a_float = parallel.rank();
    a_double = parallel.rank();
    parallel.min(a_int);
    parallel.min(a_double);
    parallel.min(a_float);
    parallel.min(a_long);
    cout<<parallel.rank()<<" , parallel.min (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.min (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.min (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.min (a_long) : "<<a_long<<endl;
    cout<<flush;
    parallel.barrier();
    COUT<<"============================================="<<endl;
    COUT<<"============================================="<<endl;

    a_int  = parallel.rank();
    a_long = parallel.rank();
    a_float = parallel.rank();
    a_double = parallel.rank();
    parallel.min_to(a_int);
    parallel.min_to(a_double);
    parallel.min_to(a_float);
    parallel.min_to(a_long);
    cout<<parallel.rank()<<" , parallel.min_to 0 (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 0 (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 0 (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 0 (a_long) : "<<a_long<<endl;
    cout<<flush;
    parallel.barrier();
    COUT<<"============================================="<<endl;
    COUT<<"============================================="<<endl;

    a_int  = parallel.rank();
    a_long = parallel.rank();
    a_float = parallel.rank();
    a_double = parallel.rank();
    parallel.min_to(a_int,1);
    parallel.min_to(a_double,1);
    parallel.min_to(a_float,1);
    parallel.min_to(a_long,1);
    cout<<parallel.rank()<<" , parallel.min_to 1 (a_int) : "<<a_int<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 1 (a_double) : "<<a_double<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 1 (a_float) : "<<a_float<<endl;
    cout<<parallel.rank()<<" , parallel.min_to 1 (a_long) : "<<a_long<<endl;
    cout<<flush;
    


    //--------------------------------------------------------
}
