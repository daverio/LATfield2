#ifndef LATFIELD2_IOSERVER_HPP
#define LATFIELD2_IOSERVER_HPP

/*! \file LATfield2_IO_server.hpp
 \brief LATfield2_IO_server.hpp contains the class IOserver definition.
 \author David Daverio
 */ 


typedef int ioserver_file

struct ioserver_file_description{
    int sctructure
    int maxSize
    
    //structured
    int dim //must be at least 2 and max 4
    int Size[4];
    hid_t datatype;
    
    
};


class IOserver {

public:
    
    
void start();
void stop();
    
int openOstream();

void closeOstream();
    
ioserver_file openFile(string filename);
    
void closeFile(ioserver_file fileID);    
    
}


#endif
