/*! file IOserver.cpp
    Created by David Daverio.
 
    A simple example of LATfield2d usage. This exec solve the poisson equation in fourier space.
    
 
 */



#include <iostream>
#include "LATfield2.hpp"

using namespace LATfield2;



int main(int argc, char **argv)
{
    int n,m;
    int io_size;
    int io_groupe_size;
    
    
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
		}
	}
	
	parallel.initialize(n,m,io_size,io_groupe_size);
    
    /*!
     Part of the code which will be executed by the IOserver process. In fact this cores should only start the server
     */
    if(parallel.isIO()) IO_Server.start();
    else
    {
        /*!
         Part of the code which will be executed by the compute process.
        
         For this example each processor will write 100 times the same line within one file with the following path: "./testfile"
         */
        
        /*!
         Create line: I am the "rank_world" MPI proces. My rank is "rank" in the compute group. I have the position ("N","M") in the process grid. 
         */
        
        string filename = "./testfile";
        ioserver_file file;
        
        string sentence;
        sentence = "I am the " + int2string(parallel.world_rank(),99999) + " MPI process. My rank is ";
        sentence += int2string(parallel.rank(),99999) + " in the compute group. I have the position (";
        sentence += int2string(parallel.grid_rank()[0],99999) + "," + int2string(parallel.grid_rank()[1],99999);
        sentence += ") in the processes grid." ;
        
        char * sendbuffer;
        sendbuffer = (char*)malloc(sentence.size()+1);
        for(int i=0;i<sentence.size();i++)sendbuffer[i]=sentence[i];
        sendbuffer[sentence.size()]='\n'; 
        
        
        
        while(IO_Server.openOstream()==OSTREAM_FAIL)usleep(50);
        
        file = IO_Server.createFile(filename);
        IO_Server.writeBuffer(file,sendbuffer,sentence.size()+1);
        IO_Server.closeFile(file);
        
        IO_Server.closeOstream();
        
        IO_Server.stop();
    }
    
    
    
    
}