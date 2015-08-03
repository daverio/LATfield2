#ifndef LATFIELD2_IOSERVER_HPP
#define LATFIELD2_IOSERVER_HPP

/*! \file LATfield2_IO_server.hpp
 \brief LATfield2_IO_server.hpp contains the class IOserver definition.
 \author David Daverio
 */ 

#include "mpi.h"
#include "int2string.hpp"

#define SERVER_STATE_TAG 1

#define SERVER_CONTROL_TAG 2

#define CONTROL_OPEN_OSTREAM 1
#define CONTROL_STOP 2

#define IO_FILE_CONTROL_TAG 3

#define CONTROL_CREATE_FILE 1
#define CONTROL_CLOSE_OSTREAM 2

#define IO_FILE_CONTROL_FILEID_TAG 4
#define IO_FILE_CONTROL_FILENAME_TAG 5


#define OSTREAM_SUCCESS 1
#define OSTREAM_FAIL 0
#define FILE_FAIL -32000
#define MAX_FILE_NUMBER 5
#define IO_FILE_CONTROL_CLOSE_TAG 15

#define FILETYPE_UNSTRUCTURED 1

//byte given;
#define IO_BUFFERS_TOTAL_SIZE 873741824 

typedef int ioserver_file;

//! A structure to describe a file for the I/O server (dedicated MPI processes for writing to disks)
struct file_struct{
    //! path to the file
    string filename;
    //! data array of the file
    char * data;
    //! size of the local part of the file
    long long size; 
    //! type of the file (currently only FILETYPE_UNSTRUCTURED)
    int type;
};



/*! \class IOserver  
 
 \brief A class to handle the I/O using MPI process reserved for IO purpose on which the files are defined
 
 This server is in beta stage, but as such a functionality is very useful, it has been added to the stable part of LATfield2. An example of the usage of this class is given in the IOserver example. User should never instanciate an IOserver object. The IOserver objet (IO_Server) is instanciate within the library header
 
 
 */
class IOserver {
    
public:
      
    ~IOserver();
    
    
  
    
    /*! \brief Server method (only called by server nodes)
     
     Method which is called to start the server.
     */
    void start(); 

    
    

    
    /*! \brief Client method (only called by compute nodes)
     
     Method which is called to stop the server.
     
     */
    void stop();  //client
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to open an Ostream. Meaning a stream from the compute to the server processes.
     
     \return OSTREAM_SUCCESS if the stream is open.
     \return OSTREAM_FAIL  if the stream cannot be open.
     */
    int openOstream();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close the current Ostream. After the stream is closed, the server will start to write the files it have in memory.
     
    */
    void closeOstream();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to create a new file, it return the fileID.
     
     \param  filename: name of the file (including the path...)
     \return  fileID.
     */
    ioserver_file createFile(string filename);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close a new file: fileID.
     
     \param ioserver_file fileID: file to close.
     */
    void closeFile(ioserver_file fileID);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to write to a file.
     !!! Beta, this method work only if fileID have been created and not closed!!!
     
     \param  fileID: file where to write data.
     \param  buffer: pointer to the buffer to add to the file fileID.
     \param  size: size of "buffer", in byte. 
     */
    void writeBuffer(ioserver_file fileID,char * buffer, int size);
  
    
    //2 steps initialisation :
    
    /*!
     Initialize the I/O server, this method is called by parallel.initialize(...). Should never be used!!!
     
     */
    void initialize(int proc_size0,int proc_size1, int IOserver_size, int IO_node_size); //called by avery cores.... initialize global variable
    
    
private:
    
    bool serverOn_flag;
    bool serverReady_flag;
    bool ostreamFile_flag;
    
    MPI_Group world_group_;
    
    MPI_Group IO_Group_;
    MPI_Group computeGroup_;
    MPI_Comm IO_Comm_;
    MPI_Comm computeComm_;
    
    MPI_Group syncLineGroup_; // root IO and compute 
    MPI_Comm  syncLineComm_;
    
    
    MPI_Group masterClientGroup_;
    MPI_Comm  masterClientComm_;
    
    MPI_Group IO_NodeGroup_;
    MPI_Comm IO_NodeComm_; 
    
    
    int IO_Rank_;
    int computeRank_;
    int syncLineRank_;
    int IO_NodeRank_;

    file_struct * files;
    
    
    int IO_ClientSize_;
    int IO_NodeSize_;
    int IO_Node_;
    
    MPI_Request sendRequest;
    
protected:
    
    char * dataBuffer;
    

};

 IOserver::~IOserver()
{
    
}


void IOserver::initialize(int proc_size0,int proc_size1, int IOserver_size, int IO_node_size)
{
    int rang[3];
    int totalMPIsize;
    int itemp;
    
    MPI_Group groupTemp1,groupTemp2;
   
    MPI_Comm_group(MPI_COMM_WORLD,&world_group_);
    MPI_Group_size(world_group_,&totalMPIsize);
    
    
    if((proc_size0*proc_size1) % IOserver_size!=0 || IOserver_size % IO_node_size!=0)
    {
        //cout<<"IOserver wrong number of process"<<endl;
        exit(-44);
    }
    
    
    rang[0]=0;
    rang[1]=proc_size0*proc_size1-1;
    rang[2]=1;
    MPI_Group_range_incl(world_group_,1,&rang,&computeGroup_);
    MPI_Comm_create(MPI_COMM_WORLD,computeGroup_ , &computeComm_);    
    
    MPI_Group_rank(computeGroup_, &computeRank_);
    
    
    rang[0]=proc_size0*proc_size1;
    rang[1]=proc_size0*proc_size1 + IOserver_size - 1 ;
    rang[2]=1;
    MPI_Group_range_incl(world_group_,1,&rang,&IO_Group_);
    MPI_Comm_create(MPI_COMM_WORLD,IO_Group_ , &IO_Comm_);
    
    MPI_Group_rank(IO_Group_, &IO_Rank_);
    
    
    
    rang[0]=proc_size0*proc_size1;
    rang[1]=0;
    MPI_Group_incl(world_group_,2,&rang[0],&syncLineGroup_);
    MPI_Comm_create(MPI_COMM_WORLD,syncLineGroup_ , &syncLineComm_);
    
    MPI_Group_rank(syncLineGroup_, &syncLineRank_);
    
    
    IO_ClientSize_=proc_size0*proc_size1/IOserver_size;
    IO_NodeSize_=IO_node_size;
    
    if(computeRank_!=MPI_UNDEFINED)itemp = floor((float)computeRank_/(float)IO_ClientSize_) * IO_ClientSize_;
    else itemp = IO_Rank_ * IO_ClientSize_;
    //if(computeRank_!=MPI_UNDEFINED)cout<< "compute core: "<< computeRank_ <<" , "<<  itemp<<endl;
    //if(IO_Rank_!=MPI_UNDEFINED)cout<< "IO core: "<< IO_Rank_ <<" , "<<  itemp<<endl;
    
    rang[0] = itemp;
    rang[1] = itemp + IO_ClientSize_ -1;
    rang[2]=1;
    MPI_Group_range_incl(world_group_,1,&rang,&groupTemp2);
    
    if(computeRank_!=MPI_UNDEFINED)itemp = proc_size0*proc_size1 + floor((float)computeRank_/(float)IO_ClientSize_);
    else itemp = proc_size0*proc_size1 + IO_Rank_;
    //if(computeRank_!=MPI_UNDEFINED)cout<< "compute core: "<< computeRank_ <<" , "<<  itemp<<endl;
    //if(IO_Rank_!=MPI_UNDEFINED)cout<< "IO core: "<< IO_Rank_ <<" , "<<  itemp<<endl;

    MPI_Group_incl(world_group_,1,&itemp,&groupTemp1);
    
    MPI_Group_union(groupTemp1,groupTemp2,&masterClientGroup_);
    MPI_Comm_create(MPI_COMM_WORLD,masterClientGroup_ , &masterClientComm_);
    
    
    //if(computeRank_!=MPI_UNDEFINED)cout<< "compute core: "<< computeRank_ <<" , "<<  itemp<<endl;
    //if(IO_Rank_!=MPI_UNDEFINED)cout<< "IO core: "<< IO_Rank_ <<" , "<<  itemp<<endl;
    
    
    
    if(IO_Rank_!=MPI_UNDEFINED)
    {
                
        itemp= floor((float)IO_Rank_/ (float)IO_node_size) * IO_node_size;
        
        rang[0] = itemp;
        rang[1] = itemp + IO_node_size -1;
        rang[2]=1;
        MPI_Group_range_incl(IO_Group_,1,&rang,&IO_NodeGroup_);
        MPI_Comm_create(IO_Comm_,IO_NodeGroup_ , &IO_NodeComm_);
        MPI_Group_rank(IO_NodeGroup_, &IO_NodeRank_);
        
        files = new file_struct[MAX_FILE_NUMBER];
        dataBuffer = (char*)malloc(IO_BUFFERS_TOTAL_SIZE);
        IO_Node_=floor((float)IO_Rank_/ (float)IO_node_size) ;
        
    }
    
    
    
    sendRequest = MPI_REQUEST_NULL;
    
}
int IOserver::openOstream()
{
    MPI_Status status;
    int flag=1;
    int count=1;
    
    MPI_Barrier(computeComm_);
    
    if(computeRank_==0)
    {
        bool getState=true;
        while(getState)
        {
            
            MPI_Iprobe(0,SERVER_STATE_TAG,syncLineComm_,&flag,&status);
            if(flag==true)
            {
                //cout<<"count : "<<count<<endl;
                count++;
                MPI_Recv(&serverReady_flag,1,MPI::BOOL,0,SERVER_STATE_TAG,syncLineComm_,&status);
                //cout<<"flag : "<< serverReady_flag<<endl;
            }
            else  getState=false;
        }
        //cout << "state : "<< serverReady_flag<<endl;
        
        if(serverReady_flag){int send=CONTROL_OPEN_OSTREAM ; MPI_Send(&send,1,MPI::INT,0,SERVER_CONTROL_TAG,syncLineComm_);}
    }
    
    
    
    MPI_Bcast(&serverReady_flag,1,MPI::BOOL,0,computeComm_);
    
    if(serverReady_flag) return OSTREAM_SUCCESS;
    else return OSTREAM_FAIL;
    
    
}
void IOserver::closeOstream()
{
    if(computeRank_==0)
    {
        int send=CONTROL_CLOSE_OSTREAM;
        MPI_Send(&send,1,MPI::INT,0,IO_FILE_CONTROL_TAG,syncLineComm_);
    }
    MPI_Barrier(computeComm_);
}
int IOserver::createFile(string filename)
{
    MPI_Status status;
    int fileID;
    MPI_Barrier(computeComm_);
    
    if(computeRank_==0)
    {
        int send=CONTROL_CREATE_FILE;
        MPI_Send(&send,1,MPI::INT,0,IO_FILE_CONTROL_TAG,syncLineComm_);
        
        MPI_Recv(&fileID,1,MPI::INT,0,IO_FILE_CONTROL_FILEID_TAG,syncLineComm_,&status);
        if(fileID!=FILE_FAIL)
        {
            int len=filename.size()+1;
            char  temp_filename[len];
            for(int i = 0;i<len;i++)temp_filename[i]=filename[i];
            
            MPI_Ssend(temp_filename,len,MPI::CHAR,0,IO_FILE_CONTROL_FILENAME_TAG,syncLineComm_);
            //cout<<len<<endl;
        }
        
    }
    
    MPI_Bcast(&fileID,1,MPI::INT,0,computeComm_);
    
    return fileID;
    
    
    
    
}
void IOserver::closeFile(int fileID)
{
    int send =fileID;
    MPI_Send(&send,1,MPI::INT,0,IO_FILE_CONTROL_CLOSE_TAG,masterClientComm_);
    MPI_Barrier(computeComm_);
}
void IOserver::writeBuffer(int fileID, char * buffer, int size)
{
    MPI_Send(buffer,size,MPI::CHAR,0,fileID,masterClientComm_);
}
void IOserver::stop()
{
    if(computeRank_==0)
    {
        int send=CONTROL_STOP;
        MPI_Send(&send,1,MPI::INT,0,SERVER_CONTROL_TAG,syncLineComm_);
    }
    MPI_Barrier(computeComm_);
}
void IOserver::start()
{
    //starting IO server (Sinc line)
    
    MPI_Status status,status2;
    MPI_Request send_statut_request;
    int flag=1;
    int control;
    
    int count=0;
    int fileID;
    int fileCounter;
    int newFilenameLength;
    int len;
    bool getingDataFlag;
    int getingData[IO_ClientSize_];
    
    
    
    
    serverOn_flag=true;
    
    while(serverOn_flag)
    {
        //cout<<count<<endl;
        
        
        ostreamFile_flag=false;
        
        if(IO_Rank_==0)
        {
            //cout<<count<<endl;
            //count++;
            
            //send non-blocking for status server free. (Sync line)
            
            serverReady_flag=true;
            MPI_Isend(&serverReady_flag,1,MPI::BOOL,1,SERVER_STATE_TAG,syncLineComm_,&send_statut_request);
        
            //test if an ostream need to be open or if the server need to be stoped (Sync line) blocking rec
      
            MPI_Recv(&control,1,MPI::INT,1,SERVER_CONTROL_TAG,syncLineComm_,&status);            
        }
        
        MPI_Bcast(&control,1,MPI::INT,0,IO_Comm_);
        //MPI_Bcast(&ostreamFile_flag,1,MPI::BOOL,0,IO_Comm_);
        
        if(control==CONTROL_STOP)
        {
            //if(IO_Rank_==0)cout<<"stop server"<<endl;
            serverOn_flag=false;
        }
        if(control==CONTROL_OPEN_OSTREAM)
        {
            //if(IO_Rank_==0)cout<<"ostream"<<endl;
            ostreamFile_flag=true;
            serverReady_flag=false;
            if(IO_Rank_==0)MPI_Isend(&serverReady_flag,1,MPI::BOOL,1,SERVER_STATE_TAG,syncLineComm_,&send_statut_request);
            fileCounter=0;
            
            
            while(ostreamFile_flag)
            {
                //cout<<"rank: "<<IO_Rank_<<", ostream open"<<endl;
                
                //test if there is a new file (Sync line)
                
                
                
                if(IO_Rank_==0)MPI_Recv(&control,1,MPI::INT,1,IO_FILE_CONTROL_TAG,syncLineComm_,&status);
                
                MPI_Bcast(&control,1,MPI::INT,0,IO_Comm_);
                
                if(control==CONTROL_CLOSE_OSTREAM)
                {
                    //if(IO_Rank_==0)cout<<"stop ostream"<<endl;
                    ostreamFile_flag=false;
                }
                
                if(control==CONTROL_CREATE_FILE)
                {
                    //create new file ID
                    if(fileCounter>= MAX_FILE_NUMBER)
                    {
                        if(IO_Rank_==0)cout<<"Already to much file"<<endl;
                        fileID=FILE_FAIL;
                    }
                    else
                    {
                        fileID=fileCounter;
                        
                        //if(IO_Rank_==0)cout<<"IO_server creating new file : "<<fileID<<endl;
                    }
                    
                    if(IO_Rank_==0)MPI_Ssend(&fileID,1,MPI::INT,1,IO_FILE_CONTROL_FILEID_TAG,syncLineComm_);
                    if(fileID  != FILE_FAIL)
                    {
                        if(IO_Rank_==0)
                        {
                            MPI_Probe(1,IO_FILE_CONTROL_FILENAME_TAG,syncLineComm_,&status);
                            MPI_Get_count(&status,MPI::CHAR,&newFilenameLength);
                        }
                        MPI_Bcast(&newFilenameLength,1,MPI::INT,0,IO_Comm_);
                        char * filename;
                        filename = (char*)malloc(newFilenameLength*sizeof(char));
                        if(IO_Rank_==0)MPI_Recv(filename,newFilenameLength,MPI::CHAR,1,IO_FILE_CONTROL_FILENAME_TAG,syncLineComm_,&status2);
                        MPI_Bcast(filename,newFilenameLength,MPI::CHAR,0,IO_Comm_);
                        files[fileID].filename=filename;
                        free(filename);
                        files[fileID].type=FILETYPE_UNSTRUCTURED; // no structured file implemented!!!!
                        //cout<< files[fileID].filename <<endl;
                        
                        if(fileID==0)files[fileID].data=dataBuffer;
                        else files[fileID].data=&(files[fileID-1].data[files[fileID-1].size]);
                        
                        files[fileID].size=0;
                        
                        fileCounter++;
                        
                        
                        if(files[fileID].type==FILETYPE_UNSTRUCTURED)
                        {
                            for(int i=0;i<IO_ClientSize_;i++)getingData[i]=1;
                            getingDataFlag=true;
                            
                            
                            while(getingDataFlag)
                            {
                                MPI_Iprobe(MPI_ANY_SOURCE,fileID, masterClientComm_,&flag,&status);
                                if(flag==true)
                                {
                                    char * send;
                                    send=&(files[fileID].data[files[fileID].size]);
                                    int size;
                                    MPI_Get_count(&status,MPI::CHAR,&size);
                                    MPI_Recv(send,size,MPI::CHAR,status.MPI_SOURCE,fileID,masterClientComm_,&status2);
                                    files[fileID].size+=size;
                                }
                                else
                                {
                                    MPI_Iprobe(MPI_ANY_SOURCE,IO_FILE_CONTROL_CLOSE_TAG,masterClientComm_ ,&flag,&status);
                                    if(flag==true)
                                    {
                                        //cout<<"file closed"<<endl;  
                                        int tempFileID;
                                        int total=0;
                                        int source = status.MPI_SOURCE;
                                        //cout<<source<<endl;
                                        MPI_Recv(&tempFileID,1,MPI::INT,source,IO_FILE_CONTROL_CLOSE_TAG,masterClientComm_,&status2);
                                        //cout<<source<<endl;
                                        getingData[source-1]=0;
                                        for(int i=0;i<IO_ClientSize_;i++)total+=getingData[i];
                                        if(total==0)getingDataFlag=false;
                                        //cout<<total<<endl;
                                    }
                                    
                                }
                            }
                            //cout<<"file closed"<<endl;
                            MPI_Barrier(IO_Comm_);
                        }
                        MPI_Barrier(IO_Comm_);
                        
                        
                        
                    }
                    
                }
                
            }//close stream and write

            
            
            if(fileCounter!=0)
            {
                for(int i=0;i < fileCounter;i++)
                {
                    //cout<<"writing file: "<< i<<endl;
                    
                    MPI_File ofile;
                    
                    long file_Offset=0;
                    long temp_long=0;
                    string str_fname;
                    str_fname = files[i].filename + int2string(IO_Node_,999)+".dat";
                    char * fname = &(str_fname[0]);
                    
                    for(int k=0;k<(IO_NodeSize_-1); k++)
                    {
                        if(IO_NodeRank_==k)
                        {
                            temp_long = file_Offset + files[i].size;
                            MPI_Send(&temp_long,1,MPI_LONG, k+1 , 0, IO_NodeComm_ );
                        }
                        if(IO_NodeRank_==k+1)
                        {
                            MPI_Recv( &file_Offset, 1, MPI_LONG, k, 0, IO_NodeComm_, &status);
                        }
                    }
                    
                    if(files[i].size!=0)
                    {
                        MPI_File_open(IO_NodeComm_,fname,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&ofile);
                        MPI_File_set_view(ofile,file_Offset,MPI_CHAR,MPI_CHAR,(char*)"native",MPI_INFO_NULL);
                        MPI_File_write_all(ofile,files[i].data,files[i].size,MPI_CHAR,&status);
                        MPI_File_close(&ofile);
                    } 
                }
            }
            fileCounter=0;
            
            MPI_Barrier(IO_Comm_);
        
        
        }
        
        
                

        
    }
}

#endif
