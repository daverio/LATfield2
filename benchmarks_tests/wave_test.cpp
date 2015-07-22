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
    double val_re, val_im;
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

#ifdef FFT3D_ACC
    double fft_time_noh2d = 0.;
#endif

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
		for (x.first(); x.test(); x.next())
		{
			phi(x) = cos(2. * M_PI * i * x.coord(0) / (double) BoxSize) * cos(2. * M_PI * j * x.coord(1) / (double) BoxSize) * cos(2. * M_PI * l * x.coord(2) / (double) BoxSize);
		}
				cpu_time_total += MPI_Wtime() - cpu_time_start;

#ifdef FULL_OUTPUT
		COUT << endl << "iteration (" << i << ", " << j << ", " << l << ")" << endl << endl;
#endif
		fft_time_start = MPI_Wtime();
		planPhi.execute(FFT_FORWARD);
		fft_time_total += MPI_Wtime() - fft_time_start;
#ifdef FFT3D_ACC
		fft_time_noh2d += planPhi.timing;
#endif

#ifdef FULL_OUTPUT
		for (rnk = 0; rnk < parallel.size(); rnk++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (parallel.rank() == rnk)
			{
				cout << " rank = " << rnk << endl;
				for (k.first(); k.test(); k.next())
				{
					if (k.coord(0) == 0)
						cout << " " << setfill('0') << setw(3) << k.coord(1) << "-" << setfill('0') << setw(3) << k.coord(2) << " ";
					
					cout << setprecision(2) << chop(phiK(k).real()/lat.sites(), 1.0e-12) << "+" << setprecision(2) << chop(phiK(k).imag()/lat.sites(), 1.0e-12) << "i ";

					if (k.coord(0) == latK.size(0)-1)
						cout << endl;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
#else
		for (k.first(); k.test(); k.next())
		{
			val_re = 0.;
			val_im = 0.;
			if (k.coord(0) == i)
			{
				if (k.coord(1) == j || k.coord(1) == BoxSize-j)
				{
					if (k.coord(2) == l || k.coord(2) == BoxSize-l)
					{
						if (i == 0 && j == 0 && l == 0)
						{
							val_re = 1.0;				
						}
						else if ((i == 0 && j == 0) || (i == 0 && l == 0) || (j == 0 && l == 0))
						{
							val_re = 0.5;
						}
						else if (i == 0 || j == 0 || l == 0)
						{
							val_re = 0.25;
						}
						else
						{
							val_re = 0.125;
						}
					}
				}
			}
			if (fabs(phiK(k).real()/lat.sites() - val_re) > TOLERANCE || fabs(phiK(k).imag()/lat.sites() - val_im) > TOLERANCE)
			{
#ifndef NO_OUTPUT
				cout << " proc#" << parallel.rank() << " : " << setfill('0') << setw(3) << k.coord(0) << "-" << setfill('0') << setw(3) << k.coord(1) << "-" << setfill('0') << setw(3) << k.coord(2) << "  " << phiK(k).real() << "+" << phiK(k).imag() << "i  exceeding tolerance!" << endl;
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
#ifdef FFT3D_ACC
	parallel.max(fft_time_noh2d);
#endif

	COUT << " timing information: field setup " << cpu_time_total << " sec, FFTs " << fft_time_total << " sec" << endl;

#ifdef FFT3D_ACC
	COUT << " FFTs w/o h2d " << fft_time_noh2d << " sec" << endl;
#endif

    parallel.sum(count);

    exit(count);
}

