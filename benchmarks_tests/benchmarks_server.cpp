/*! file benchmarks.cpp
    Created by David Daverio.
 
    LATfield2d Benchmarks.
 
 */

#include <unistd.h>

#include <iostream>
#include "LATfield2.hpp"

using namespace LATfield2;



int main(int argc, char **argv)
{
    int n,m;
    int io_size;
    int io_groupe_size;
    int runs=10;
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
            case 'i':
                io_size =  atoi(argv[++i]);
                break;
            case 'g':
                io_groupe_size = atoi(argv[++i]);
                break;
            case 'r':
                runs = atoi(argv[++i]);
                break;
            case 'o':
                str_filename = argv[++i];
                break;
		}
	}

	
    parallel.initialize(n,m,io_size,io_groupe_size);
    
    if(parallel.isIO())ioserver.start();
    else
    {
               
        
        double * buffer;
        long size = 32l*1024l*1024l*sizeof(double);
        long total_ubin = 32l*sizeof(double);
        parallel.sum(total_ubin); // total in Mb
        total_ubin*=runs;
        
        buffer = (double*)malloc(size);
        
        
        Lattice lat(3,1024,0);
        Field<double> phi(lat,3);
        
        hsize_t size_dset[lat.dim()];
        hsize_t offset[lat.dim()];
        
        
        for(int i=0;i<lat.dim();i++)size_dset[i] = lat.sizeLocal(i);
        
        offset[0]=0;
        offset[1]=lat.coordSkip()[1];
        offset[2]=lat.coordSkip()[0];
        
        
        
        
        int nparts=1000000;
        part_simple pcls[nparts];
        part_simple_dataType pcls_types;
        
        
        
        
        ioserver_file ubin_file,uh5_file,sh5_file;
        
        double timerRef,timerSend2Server_ubin, timerSend2Server_uh5, timerSend2Server_sh5, timerSend2Server_fieldSlice;
        
        
        
        
        
        while(!ioserver.openOstream())usleep(50);
        
        ubin_file = ioserver.openFile("./server/TestServer_ubin",UNSTRUCTURED_BIN_FILE);
        uh5_file = ioserver.openFile("./server/TestServer_uh5",UNSTRUCTURED_H5_FILE,pcls_types.part_memType,pcls_types.part_fileType);
        sh5_file = ioserver.openFile("./server/TestServer_sh5",STRUCTURED_H5_FILE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,&lat,3,1);
        phi.saveHDF5_server_open("testPhi",10,64);
        
        
        
        
        timerRef = MPI_Wtime();
        for(int i=0;i<runs;i++)ioserver.sendData(ubin_file,(char*)buffer,size);
        timerSend2Server_ubin = MPI_Wtime()-timerRef;
        
        timerRef = MPI_Wtime();
        for(int i=0;i<runs;i++)ioserver.sendData(uh5_file,(char*)pcls,size*sizeof(part_simple));
        timerSend2Server_uh5 = MPI_Wtime()-timerRef;
        
        timerRef = MPI_Wtime();
        ioserver.sendData(sh5_file,(char*)phi.data(),size_dset,offset);
        timerSend2Server_sh5 = MPI_Wtime()-timerRef;
        
        timerRef = MPI_Wtime();
        phi.saveHDF5_server_write(4);
        timerSend2Server_fieldSlice = MPI_Wtime()-timerRef;
        
        
        ioserver.closeFile(ubin_file);
        ioserver.closeFile(uh5_file);
        ioserver.closeFile(sh5_file);
        
        ioserver.closeOstream();
        
        ofstream textfile;
        
        if(parallel.isRoot())
        {
        textfile.open(str_filename.c_str(),ios::out | ios::app); 
        textfile<< n<<","<<m<< ","<< io_size<< ","<<io_groupe_size<< ","<<timerSend2Server_ubin<< ","<<total_ubin<<endl;
        textfile.close();
        }
        ioserver.stop();
    }
    
    
    
    
}
   
