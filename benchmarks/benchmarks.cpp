/*! file benchmarks.cpp
    Created by David Daverio.
 
    LATfield2d Benchmarks.
 
 */

#include <unistd.h>

#include <iostream>
#include "LATfield2d.hpp"

using namespace LATfield2d;



int main(int argc, char **argv)
{
    int n,m;
    int BoxSize=64;
    int runs=3;
    double maxTime=300;
    string str_filename;
    
    
    
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
				BoxSize =  atoi(argv[++i]);
				break;
            case 'r':
                runs = atoi(argv[++i]);
                break;
            case 'o':
                str_filename = argv[++i];
                break;
		}
	}

	
    parallel.initialize(n,m);
    
    
    int halo = 1;
    int khalo =0;
    int dim = 3;
    int comp[4] = {1,2,3,6};
    double sigma2=0.5;
    int performedRuns;
    int halos[4]={1,2,3,4};
    
    
    ioserver_file file;
    ofstream textfile;
    
    
    double timerFFTreal[2][4]={{0,0,0,0},{0,0,0,0}};
    double timerFFTImag[2][4]={{0,0,0,0},{0,0,0,0}};
    
    double timerFillGaussian[4]={0,0,0,0};
    double timerUpDateHalo[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    double timerNeighBorOpp[4]={0,0,0,0};
    double timerSaveHDF5[4] = {0,0,0,0};
    double timerSend2Server = 0;
    
    
    double timerRef;
    double timerRef2;
    
    Lattice lat(dim,BoxSize,halo);
    Lattice latKReal,latKImag;
    
    latKReal.initializeRealFFT(lat, khalo);
    latKImag.initializeComplexFFT(lat, khalo);
    
    Site x(lat);
    rKSite kReal(latKReal);
    cKSite kImag(latKImag);
    
    
    //fourier field
    Field<Real> phiReal;
    Field<Imag> phiImag;
    
    Field<Imag> phiKReal;
    Field<Imag> phiKImag;
    
    PlanFFT<Imag> planReal;
    PlanFFT<Imag> planImag;
    
    //real space field
    Field<Real> lapRho;
    Field<Real> rho;
    
    
    COUT << "FFT Benchmark"<<endl;
    for(int i = 0;i<4;i++)
    {
        
        //COUT<< "working with "<<comp[i]<<" components."<<endl;
        
        phiReal.initialize(lat,comp[i]);
        phiKReal.initialize(latKReal,comp[i]);
        
        phiImag.initialize(lat,comp[i]);
        phiKImag.initialize(latKImag,comp[i]);
        
        
        planImag.initialize(&phiImag,&phiKImag);
        planReal.initialize(&phiReal,&phiKReal);
        
        for(x.first();x.test();x.next())
        {
            double x2 = pow(0.5l + x.coord(0) - lat.size(0)/2,2);
            x2 += pow(0.5l + x.coord(1) - lat.size(1)/2,2);
            x2 += pow(0.5l + x.coord(2) - lat.size(2)/2,2);
            double r=1.0 * exp(-x2/sigma2);
            for(int c=0;c<comp[i];c++){
                phiReal(x,c)= r;
                phiImag(x,c)=Imag(r,r);
            }
        }
        
        performedRuns =0;
        timerRef2 = MPI_Wtime();
        for(int nb=0;nb<runs;nb++)
        {
            
            timerRef = MPI_Wtime();
            planReal.execute(FFT_FORWARD);
            timerFFTreal[0][i]+=MPI_Wtime() - timerRef;
            
            timerRef = MPI_Wtime();
            planImag.execute(FFT_FORWARD);
            timerFFTImag[0][i]+=MPI_Wtime() - timerRef;
            
            
            timerRef = MPI_Wtime();
            planReal.execute(FFT_BACKWARD);
            timerFFTreal[1][i]+=MPI_Wtime() - timerRef;
            
            
            timerRef = MPI_Wtime();
            planImag.execute(FFT_BACKWARD);
            timerFFTImag[1][i]+=MPI_Wtime() - timerRef;
            
            performedRuns++;
            if((MPI_Wtime() - timerRef2)>maxTime)nb=runs;
            
        }
        
        timerFFTreal[0][i]/=performedRuns;
        timerFFTImag[0][i]/=performedRuns;
        timerFFTreal[1][i]/=performedRuns;
        timerFFTImag[1][i]/=performedRuns;
        
        /*
        COUT << performedRuns << endl;
        COUT << timerFFTreal[0][i]<<endl;
        COUT << timerFFTImag[0][i]<<endl;
        COUT << timerFFTreal[1][i]<<endl;
        COUT << timerFFTImag[1][i]<<endl;*/
        
        phiReal.dealloc();
        phiKReal.dealloc();
        phiImag.dealloc();
        phiKImag.dealloc();
    }
    //cout << parallel.rank()<<" okokok"<<endl;
    
    COUT<<"Opp benchmark"<<endl;
    
    for(int i = 0;i<4;i++)
    {
        rho.initialize(lat,comp[i]);
        rho.alloc();
        lapRho.initialize(lat,comp[i]);
        lapRho.alloc();
        
        
        
        performedRuns =0;
        timerRef2 = MPI_Wtime();
        for(int nb=0;nb<runs;nb++)
        {
            timerRef=MPI_Wtime();
            for(x.first();x.test();x.next())
            {
                double x2 = pow(0.5l + x.coord(0) - lat.size(0)/2,2);
                x2 += pow(0.5l + x.coord(1) - lat.size(1)/2,2);
                x2 += pow(0.5l + x.coord(2) - lat.size(2)/2,2);
                for(int c=0;c<comp[i];c++)rho(x,c)= 1.0 * exp(-x2/sigma2);
            }
            timerFillGaussian[i] += MPI_Wtime()-timerRef;
            
            timerRef=MPI_Wtime();
            for(int c=0;c<comp[i];c++){for(int l=0;l<3;l++)lapRho(x,c) += rho(x+l,c) - 2 * rho(x,c) + rho(x-l,c);}
            timerNeighBorOpp[i]+= MPI_Wtime()-timerRef;
            
            performedRuns++;
            if((MPI_Wtime() - timerRef2)>maxTime)nb=runs;
        }
        
        timerFillGaussian[i]/=performedRuns;
        timerNeighBorOpp[i]/=performedRuns;
        
        
        performedRuns =0;
        timerRef2 = MPI_Wtime();
        for(int nb=0;nb<runs;nb++)
        {
            timerRef=MPI_Wtime();
            rho.saveHDF5("./hdf5/test.h5");
            timerSaveHDF5[i] += MPI_Wtime()-timerRef;
            
            
            performedRuns++;
            if((MPI_Wtime() - timerRef2)>maxTime)nb=runs;
        }
        timerSaveHDF5[i]/=performedRuns;
        
        
        rho.dealloc();
        lapRho.dealloc();
    }
    
    //updatehalo benchmark
    
    COUT<<"Halo"<<endl;
    
    for(int h=0;h<4;h++)
    {
        lat.initialize(3,BoxSize,halos[h]);
        Site X(lat);
        
        for(int i = 0;i<4;i++)
        {
            
            rho.initialize(lat,comp[i]);
            rho.alloc();
            
            
            for(int c=0;c<comp[i];c++)rho(x,c)= 1.0; 
            
            timerRef2 = MPI_Wtime();
            for(int nb=0;nb<runs;nb++)
            {
                timerRef = MPI_Wtime();
                rho.updateHalo();
                timerUpDateHalo[h][i]+= MPI_Wtime()-timerRef;
                
                
                performedRuns++;
                if((MPI_Wtime() - timerRef2)>maxTime)nb=runs;
            }
            timerUpDateHalo[h][i]/=performedRuns;
            //COUT<< timerUpDateHalo[h][i] << endl;
            rho.dealloc();
            
        }
        
        
    }
    
    
    if(parallel.isRoot())
    {
        textfile.open(str_filename.c_str(),ios::out | ios::app); 
        textfile<< n<<","<<m<< ","<<BoxSize<< "," << io_size << ","<< io_groupe_size<<",";
        textfile<< timerFFTreal[0][0]<< ","<< timerFFTreal[0][1]<< ","<< timerFFTreal[0][2]<< ","<< timerFFTreal[0][3]<< ",";
        textfile<< timerFFTreal[1][0]<< ","<< timerFFTreal[1][1]<< ","<< timerFFTreal[1][2]<< ","<< timerFFTreal[1][3]<< ",";
        textfile<< timerFFTImag[0][0]<< ","<< timerFFTImag[0][1]<< ","<< timerFFTImag[0][2]<< ","<< timerFFTImag[0][3]<< ",";
        textfile<< timerFFTImag[1][0]<< ","<< timerFFTImag[1][1]<< ","<< timerFFTImag[1][2]<< ","<< timerFFTImag[1][3]<< ",";
        textfile<< timerFillGaussian[0]<< ","<< timerFillGaussian[1]<< ","<< timerFillGaussian[2]<< ","<< timerFillGaussian[3]<< ",";
        textfile<< timerNeighBorOpp[0]<< ","<< timerNeighBorOpp[1]<< ","<< timerNeighBorOpp[2]<< ","<< timerNeighBorOpp[3]<< ",";
        textfile<< timerUpDateHalo[0][1]<< ","<< timerUpDateHalo[0][1]<< ","<< timerUpDateHalo[0][2]<< ","<< timerUpDateHalo[0][3]<< ",";
        textfile<< timerUpDateHalo[1][1]<< ","<< timerUpDateHalo[1][1]<< ","<< timerUpDateHalo[1][2]<< ","<< timerUpDateHalo[1][3]<< ",";
        textfile<< timerUpDateHalo[2][1]<< ","<< timerUpDateHalo[2][1]<< ","<< timerUpDateHalo[2][2]<< ","<< timerUpDateHalo[2][3]<< ",";
        textfile<< timerUpDateHalo[3][1]<< ","<< timerUpDateHalo[3][1]<< ","<< timerUpDateHalo[3][2]<< ","<< timerUpDateHalo[3][3]<< ",";
        textfile<< timerSaveHDF5[0]<< ","<< timerSaveHDF5[1]<< ","<< timerSaveHDF5[2]<< ","<< timerSaveHDF5[3]<< endl;
        textfile.close();
    }



}
   
