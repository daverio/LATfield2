#include <iostream>
#include <iomanip>
#include "LATfield2.hpp"

using namespace LATfield2;

#define TOLERANCE 1.0e-10

double chop(const double val, const double tol)
{
	if (val > tol || val < -tol) return val;
	else return 0.;
}

int main(int argc, char **argv)
{
    int n,m;
    int BoxSize = 64;
    int halo = 1;
    int khalo =0;
    int dim = 3;
    int comp = 1;
    int i,j,l,rnk;
    double val;
    int count = 0;


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
            case 'b':
                BoxSize = atoi(argv[++i]);
                break;
		}
	}

    parallel.initialize(n,m);

    Lattice lat(dim,BoxSize,halo);
    Lattice latK;
    latK.initializeRealFFT(lat, khalo);

    Site x(lat);
    rKSite k(latK);

    Field<Real> phi;
    phi.initialize(lat,comp);
    Field<Imag> phiK;
    phiK.initialize(latK,comp);
    PlanFFT<Imag> planPhi(&phi,&phiK);

	double cpu_time_start;
	double cpu_time_total = 0.;
	double fft_time_start;
	double fft_time_total = 0.;


#ifndef LOOP_COUNT
#define LOOP_COUNT 4
#endif 

    for (i = 0; i < LOOP_COUNT; i++)
    {
         for (j = 0; j < LOOP_COUNT; j++)
         {
             for (l = 0; l < LOOP_COUNT; l++)
             {
				cpu_time_start = MPI_Wtime();
		for (k.first(); k.test(); k.next())
		{
			phiK(k) = Imag(cos(2. * M_PI * i * k.coord(0) / (double) BoxSize) * cos(2. * M_PI * j * k.coord(1) / (double) BoxSize) * cos(2. * M_PI * l * k.coord(2) / (double) BoxSize), 0.);
		}
				cpu_time_total += MPI_Wtime() - cpu_time_start;

#ifdef FULL_OUTPUT
		COUT << endl << "iteration (" << i << ", " << j << ", " << l << ")" << endl << endl;
#endif
		fft_time_start = MPI_Wtime();
		planPhi.execute(FFT_BACKWARD);
		fft_time_total += MPI_Wtime() - fft_time_start;

#ifdef FULL_OUTPUT
		for (rnk = 0; rnk < parallel.size(); rnk++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (parallel.rank() == rnk)
			{
				cout << " rank = " << rnk << endl;
				for (x.first(); x.test(); x.next())
				{
					if (x.coord(0) == 0)
						cout << " " << setfill('0') << setw(3) << x.coord(1) << "-" << setfill('0') << setw(3) << x.coord(2) << " ";
					
					cout << setprecision(2) << chop(phi(x)/lat.sites(), 1.0e-12) << " ";

					if (x.coord(0) == lat.size(0)-1)
						cout << endl;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
#else
		for (x.first(); x.test(); x.next())
		{
			val = 0.;
			if (x.coord(0) == i || x.coord(0) == BoxSize-i)
			{
				if (x.coord(1) == j || x.coord(1) == BoxSize-j)
				{
					if (x.coord(2) == l || x.coord(2) == BoxSize-l)
					{
						if (i == 0 && j == 0 && l == 0)
						{
							val = 1.0;				
						}
						else if ((i == 0 && j == 0) || (i == 0 && l == 0) || (j == 0 && l == 0))
						{
							val = 0.5;
						}
						else if (i == 0 || j == 0 || l == 0)
						{
							val = 0.25;
						}
						else
						{
							val = 0.125;
						}
					}
				}
			}
			if (fabs(phi(x)/lat.sites() - val) > TOLERANCE)
			{
#ifndef NO_OUTPUT
				cout << " proc#" << parallel.rank() << " : " << setfill('0') << setw(3) << x.coord(0) << "-" << setfill('0') << setw(3) << x.coord(1) << "-" << setfill('0') << setw(3) << x.coord(2) << "  " << phi(x)/lat.sites() << " exceeding tolerance!" << endl;
#endif
				count++;
			}
		}
#endif
             }
         }
    }

	parallel.max(cpu_time_total);
	parallel.max(fft_time_total);

	COUT << " timing information: field setup " << cpu_time_total << " sec, FFTs " << fft_time_total << " sec" << endl;


    parallel.sum(count);

    exit(count);
}

