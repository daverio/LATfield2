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
    
    if(parallel.isIO())IO_Server.start();
    else
    {
               
        
        double * buffer;
        long size = 4l*1024l*1024l*sizeof(double);
        long total = 4l*sizeof(double);
        parallel.sum(total); // total in Mb
        total*=runs;
        
        buffer = (double*)malloc(size);
        
        
        
        while(IO_Server.openOstream()==OSTREAM_FAIL)usleep(50);
        
        timerRef = MPI_Wtime();
        file = IO_Server.createFile("./server/TestServer");
        for(int i=0;i<runs;i++)IO_Server.writeBuffer(file,(char*)buffer,size);
        IO_Server.closeFile(file);
        timerSend2Server = MPI_Wtime()-timerRef;
        
        COUT<< timerSend2Server << endl;
        
        IO_Server.closeOstream();
        
        
        if(parallel.isRoot())
        {
        textfile.open(str_filename.c_str(),ios::out | ios::app); 
        textfile<< n<<","<<m<< ","<< io_size<< ","<<io_groupe_size<< ","<<timerSend2Server<< ","<<total<<endl;
        textfile.close();
        }
        IO_Server.stop();
    }
    
    
    
    
}
   
