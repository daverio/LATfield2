#ifndef LATFIELD2_IOSERVER_HPP
#define LATFIELD2_IOSERVER_HPP

/*! \file LATfield2_IO_server.hpp
 \brief LATfield2_IO_server.hpp contains the class IOserver definition.
 \author David Daverio
 */

#include "int2string.hpp"
#include <cstring>
#include "mpi.h"


#define SERVER_READY          0
#define SERVER_BUSY           1

#define OSTREAM_SUCCESS 0
#define OSTREAM_FAIL    1


#define CONTROL_STOP          0
#define CONTROL_OPEN_OSTREAM  1
#define CONTROL_CLOSE_OSTREAM 2




#define SERVER_STATE_TAG   0
#define SERVER_CONTROL_TAG 1


//tags 0 to 9999 reserved for files_ID

#define FILE_OPEN_TAG      10000
#define FILE_CLOSE_TAG     10001
#define CLOSE_OSTREAM_TAG  10002
#define GET_DATA_TAG       10003
#define GET_ATTR_TAG       10004
#define GET_DSET_TAG       10005
#define GET_HEADER_TAG     10006
#define OFFSET_ATTR_TAG    10050
#define OFFSET_DSET_TAG    20050

#define UNSTRUCTURED_BIN_FILE  0
#define UNSTRUCTURED_H5_FILE   1
#define STRUCTURED_H5_FILE     2

#define DEFAULT_SEND  0
#define BUFFERED_SEND 1




#define IO_CORE_BUFFER_SIZE  1073741824


#define MPI_BSEND_BUFFER_SIZE 40960



/*! \struct ioserver_file
    \brief Structure describing a output server file (used by compute processes)
 */
struct ioserver_file
{
    int ID;
    bool is_open;
    int type;
    int sizeof_mem_dtype; //used only by h5 files
    int dim; //used for h5 structured file
    int components; //used for h5 structured file
    int array_size; //used for h5 structured file
    ioserver_file(){
        ID=-1;
        type=-1;
        is_open = false;
        sizeof_mem_dtype = -1;
        dim =-1;
        components=-1;
        array_size=-1;
    }
};


#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct h5_attr
{
    string name;
    hid_t dtype;
    hsize_t size;
    char * attr;
};
struct h5_dset
{
    string name;
    hid_t dtype;
    hsize_t dim;
    hsize_t * size;
    char * dset;
};
struct unstruct_message
{
    char * data;
    long size;
    int core;
};
bool operator<(struct unstruct_message m,struct unstruct_message m1){return (m.core<m1.core);}

struct unstruct_bin_file
{
    char * filename;
    int fnSize;
    int ID;
    int isClosed;
    int sizeof_mem_dtype;
    int * isClosed_client;
    int totalSize;
    int headerSize;
    char * header;
    list<unstruct_message> msg;
};
bool operator<(struct unstruct_bin_file f,struct unstruct_bin_file f1){return (f.ID<f1.ID);}

struct unstruct_h5_file : unstruct_bin_file
{
    hid_t mem_dtype;
    hid_t file_dtype;
    list<h5_attr> attr;
    list<h5_dset> dset;
};
bool operator<(struct unstruct_h5_file f,struct unstruct_h5_file f1){return (f.ID<f1.ID);}

struct struct_message
{
    char * data;
    hsize_t * size;
    hsize_t * offset;
};
struct struct_h5_file
{
    char * filename;
    int fnSize;
    int ID;
    int isClosed;
    int * isClosed_client;
    int sizeof_mem_dtype;
    
    hid_t mem_dtype;
    hid_t file_dtype;
    int dim;
    hsize_t * size;
    hsize_t offset;
    int components;
    int array_size;
    list<struct_message> msg;
    list<h5_attr> attr;
    list<h5_dset> dset;
};
bool operator<(struct struct_h5_file f,struct struct_h5_file f1){return (f.ID<f1.ID);}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */



/*! \class IOserver
 
 \brief \brief A class to handle the outputs using MPI process reserved for output purpose.
 
 The output server is a group of processes reserved for output purpose. It's usage is exposed in details by the Ioserver example.
 As for the parallel object, users should never instanciate an IOserver object as the IOserver objet (ioserver) is instanciate within the library header.

 */

class IOserver {

public:
    ~IOserver();
    /*!
     Initialize the I/O server, this method is called by parallel.initialize(...). Should never be called. When using the server, the method parallel.initalize() is used to initialize both the parallel and ioserver objects.
     
     */
    void initialize(int proc_size0,int proc_size1,
                    int IOserver_size, int IO_node_size);
    
    /*! \brief Server method (only called by server nodes)
     
     Method which is called to start the server.
     */
    void start();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method which is called to stop the server.
     
     */
    void stop();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to open an Ostream. Meaning a stream from the compute to the server processes.
     
     \return OSTREAM_SUCCESS if the stream is open.
     \return OSTREAM_FAIL  if the stream cannot be open.
     */
    int  openOstream();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close the current Ostream. After the stream is closed, the server will start to write the files it have in memory.
     
     */
    void closeOstream();
    
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to create a new file, it return the fileID.
     
     \param  string filename: name of the file (including the path...)
     \param  int type: type of the file. Can take the value: UNSTRUCTURED_BIN_FILE, UNSTRUCTURED_H5_FILE or STRUCTURED_H5_FILE.
     \param  hid_t mem_datatype: used only for HDF5 files, hdf5 datatype used for memory.
     \param  hid_t file_datatype: used only for HDF5 files, hdf5 datatype used for memory.
     \param  LATfield2::Lattice * lat: used only for structured HDF5 files, pointer to a Lattice object. The lattice describe the dataset used for the file, and its distribution in the compute process grid.
     \param  int components: used only for structured HDF5 files, number of elements per lattice per lattice site.
     \param  int array_size: used only for structured HDF5 files, each elements/components can be an array of a given type, size of this array.
     \return  oiserver_file
     */
    ioserver_file openFile(string file_name,
                           int type,
                           hid_t mem_datatype = -1,
                           hid_t file_datatype = -1,
                           LATfield2::Lattice * lat = NULL,
                           int components = 1,
                           int array_size = 1);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close a given file.
     
     \param ioserver_file file: file to close.
     */
    void closeFile(ioserver_file file);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method send data to unstructured files (UNSTRUCTURED_BIN_FILE and UNSTRUCTURED_H5_FILE).
     
     \param ioserver_file file: destination file
     \param char * message: pointer to the array containing the data to be send.
     \param long size: size of the message array, in byte.
     \param int sendType: method use for the MPI send. Default is: DEFAULT_SEND, using MPI_Send. BUFFERED_SEND will use MPI_Bsend
     */
    void sendData(ioserver_file file, char * message,long size,int sendType);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method send data to an UNSTRUCTURED_H5_FILE files.
     
     \param ioserver_file file: destination file
     \param char * message: pointer to the array containing the data to be send. Stored continuously, in row major order.
     \param hsize_t * size: pointer to an array contining the size of each demension of the sub-lattice which is sent.
     \param hsize_t * offset: offset of the sub-lattice in respect to the entire lattice (global coordinate of the lowest site of the message.).
     \param int sendType: method use for the MPI send. Default is: DEFAULT_SEND, using MPI_Send. BUFFERED_SEND will use MPI_Bsend
     */
    void sendData(ioserver_file file, char * message,hsize_t * size, hsize_t * offset,int sendType);

    
    
    /*! \brief Client method (only called by compute nodes)
     
     Method send a header to an UNSTRUCTURED_BIN_FILE files. The header is sent only the root compute process of the server nodes. This list of processes can be selected using: if(isClientFileRoot_){}
     Headers are send using MPI_Bsend.
     
     \param ioserver_file file: destination file
     \param char * header: pointer to the array containing the header to send.
     \param int size: size of the header in byte.
     */
    void sendHeader(ioserver_file file,char* header,int size);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method send a attribute to HDF5 files. The attribute is sent only the root compute process of the server nodes. This list of processes can be selected using: if(isClientFileRoot_){}
     Attributes are send using MPI_Bsend.
     
     \param ioserver_file file: destination file
     \param string attr_name: string containing the name of the attribute to write.
     \param char * attr: pointer to the array containing the attribute to send.
     \param int size: size of the attr array. Given in units of H5Tget_size(dtype).
     \param hid_t dtype: HDF5 datatype of the attribute.
     */
    void sendATTR(ioserver_file file,string attr_name,char * attr, int size,hid_t dtype);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method send a dataset to an HDF5 files. The dataset is sent only the root compute process of the server nodes. This list of processes can be selected using: if(isClientFileRoot_){}
     Attributes are send using MPI_Bsend.
     
     \param ioserver_file file: destination file
     \param string dset_name: string containing the name of the dataset to write.
     \param char * dset: pointer to the array containing the dataset to send.
     \param hsize_t dim: number of dismension of the dataset.
     \param hsize_t * size: size of each dimension of the dataset.
     \param hid_t dtypee: HDF5 datatype of the attribute.
     */
    void sendDataset(ioserver_file file,string dset_name, char * dset, hsize_t dim, hsize_t * size, hid_t dtype);
    
    /*! \brief Client method (only called by compute nodes)
     
     \return MPI_Comm compute_file_comm_: This communicator contains all compute processes linked to the same IO node. All compute processes contained in this communicator will write within the same files. This communicator is used to broadcast the ID of a file to each processes which share the same file. It is used by the openFile(...) method of the IOserver class.

     */
    MPI_Comm compute_file_comm(){return compute_file_comm_;};
    /*! \brief Client method (only called by compute nodes)
     
     \return int compute_file_size_: size of the compute_file_comm_ communicator.
     \sa compute_file_comm()
     
     */
    int compute_file_size(){return compute_file_size_;};
    
    /*! \brief Client method (only called by compute nodes)
     
     \return int compute_file_rank_: rank in the compute_file_comm_ communicator.
     \sa compute_file_comm()
     
     */
    int compute_file_rank(){return compute_file_rank_;};
    /*! \brief Client & server method
     
     \return int io_node_number_: number of server nodes.
     */
    int io_node_number(){return io_node_number_;};
    /*! \brief Client & server method
     
     \return int io_node_number_: return the node in which this process is.
     */
    int my_node(){return my_node_;};
    
    /*! \brief Client method (only called by compute nodes)
     
     \return bool isClientFileRoot: return true if the process is the root of the compute_file_comm_.
     */
    bool isClientFileRoot(){return isClientFileRoot_;};
    
private:
    
    
    list<unstruct_bin_file> usb_files_;
    list<unstruct_h5_file> ush_files_;
    list<struct_h5_file> sh_files_;
    
    void ostream();
    void write_files();
    
    MPI_Comm  compute_world_comm_;
    MPI_Group compute_world_group_;
    
    MPI_Comm  io_world_comm_;
    MPI_Group io_world_group_;
    
    MPI_Comm  compute_file_comm_;
    MPI_Group compute_file_group_;
    
    MPI_Comm  io_node_comm_;
    MPI_Group io_node_group_;
    
    
    MPI_Comm  client_comm_;
    MPI_Group client_group_;
    
    MPI_Comm  sync_global_comm_;
    MPI_Group sync_global_group_;
    
    
    int io_node_number_;
    int io_node_size_;
    int my_node_;
    int io_node_rank_;
    int compute_file_size_;
    int compute_file_rank_;
    int client_size_;
    
    int proc_size0_;
    int proc_size1_;
    
    
    bool isIO_;
    bool isRoot_;
    bool isClientRoot_;
    bool isClientFileRoot_;
    bool isIONodeRoot_;
    
    int state_;
    
    int filesNumber_;
    
    char * data_;
    char * wp_data_;
    
    char * sendBuffer_;
    
    
};


IOserver::~IOserver()
{
    int * size;
    *size =MPI_BSEND_BUFFER_SIZE;
    MPI_Buffer_detach(sendBuffer_,size);
    free(sendBuffer_);
    //cout<<size<<endl;
    delete[] data_;
}

void IOserver::initialize(int proc_size0,int proc_size1,
                          int IOserver_size, int IO_node_size)
{
    
    if((proc_size0*proc_size1) % IOserver_size!=0)
    {
        cout<<"IOserver: wrong number of total process"<<endl;
        exit(-44);
    }
    if( IOserver_size % IO_node_size!=0)
    {
        cout<<"IOserver: IOserver_size % IO_node_size!=0"<<endl;
        exit(-44);
    }
    if(proc_size1 % (IOserver_size/IO_node_size) !=0)
    {
        cout<<"IOserver: proc_size1 % (IOserver_size/IO_node_size)!=0"<<endl;
        exit(-44);
    }
    if(IOserver_size % proc_size1 !=0)
    {
        cout<<"IOserver: IOserver_size % proc_size1 !=0"<<endl;
        exit(-44);
    }
    if(proc_size0 % (IOserver_size / proc_size1) !=0)
    {
        cout<<"IOserver: proc_size0 % (IOserver_size / proc_size1) !=0"<<endl;
        exit(-44);

    }
    if(IO_node_size % (IOserver_size / proc_size1) !=0)
    {
        cout<<"IOserver: IO_node_size % (IOserver_size / proc_size1) !=0"<<endl;
        exit(-44);
        
    }
    
    data_ = new char[IO_CORE_BUFFER_SIZE];
    sendBuffer_ = (char*)malloc(MPI_BSEND_BUFFER_SIZE);
    MPI_Buffer_attach(sendBuffer_,MPI_BSEND_BUFFER_SIZE);
    
    proc_size0_ = proc_size0;
    proc_size1_ = proc_size1;
    
    MPI_Group groupTemp1,groupTemp2;
    
    //int world_rank;
    int rank;
    int compute_rank;
    int io_world_rank;
    int rang[3];
    
    MPI_Group world_group;
    
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);
    
    rang[0]=0;
    rang[1]=(proc_size0*proc_size1)-1;
    rang[2]=1;
    MPI_Group_range_incl(world_group,1,&rang,&compute_world_group_);
    MPI_Comm_create(MPI_COMM_WORLD,compute_world_group_ , &compute_world_comm_);
    
    MPI_Group_rank(compute_world_group_, &compute_rank);
    
    
    rang[0]=proc_size0*proc_size1;
    rang[1]=(proc_size0*proc_size1) + IOserver_size - 1 ;
    rang[2]=1;
    MPI_Group_range_incl(world_group,1,&rang,&io_world_group_);
    MPI_Comm_create(MPI_COMM_WORLD,io_world_group_ , &io_world_comm_);
    
    MPI_Group_rank(io_world_group_, &io_world_rank);
    if(io_world_rank != MPI_UNDEFINED)isIO_=true;
    
    io_node_size_ = IO_node_size;
    io_node_number_ = IOserver_size / IO_node_size;
    compute_file_size_ = proc_size0*proc_size1 / io_node_number_;
    client_size_ = proc_size0*proc_size1 / IOserver_size;
    
    if(io_world_rank != MPI_UNDEFINED)
        my_node_ =  io_world_rank / IO_node_size;
    else
        my_node_ = compute_rank / compute_file_size_;

    rang[0]= (proc_size0*proc_size1) + (my_node_ * io_node_size_) ;
    rang[1]= (proc_size0*proc_size1) + ((my_node_+1) * io_node_size_) -1;
    rang[2]=1;
    MPI_Group_range_incl(world_group,1,&rang,&io_node_group_);
    MPI_Comm_create(MPI_COMM_WORLD,io_node_group_ , &io_node_comm_);
    
    MPI_Group_rank(io_node_group_, &io_node_rank_);
    if(io_node_rank_==0)isIONodeRoot_=true;
    else isIONodeRoot_=false;
    
    rang[0]= my_node_ * compute_file_size_ ;
    rang[1]= ((my_node_+1) * compute_file_size_) -1;
    rang[2]=1;
    MPI_Group_range_incl(world_group,1,&rang,&compute_file_group_);
    MPI_Comm_create(MPI_COMM_WORLD,compute_file_group_ , &compute_file_comm_);
    
    MPI_Group_rank(compute_file_group_, &rank);
    if(rank!=MPI_UNDEFINED)compute_file_rank_=rank;
    else compute_file_rank_ = -1;
    
    if(rank==0)isClientFileRoot_ =true;
    else isClientFileRoot_=false;

    
    
    rang[0]=proc_size0*proc_size1;
    rang[1]=0;
    MPI_Group_incl(world_group,2,&rang[0],&sync_global_group_);
    MPI_Comm_create(MPI_COMM_WORLD,sync_global_group_ , &sync_global_comm_);
    
    MPI_Group_rank(sync_global_group_, &rank);
    if(rank!=MPI_UNDEFINED)isRoot_=true;
    else isRoot_=false;
    
    
    int which_master;
    
    if(io_world_rank != MPI_UNDEFINED) which_master = io_world_rank;
    else which_master = compute_rank / client_size_;
    
    //cout<< compute_rank<< " , "<< io_world_rank<<" , "<< my_node_<< " , " << which_master <<endl;
    
    rang[0]=  which_master * client_size_;
    rang[1]=  ((which_master+1) * client_size_)-1;
    rang[2]= 1;
    MPI_Group_range_incl(world_group,1,&rang,&groupTemp1);
    
    rang[0]=(proc_size0*proc_size1) + which_master;
    MPI_Group_incl(world_group,2,&rang[0],&groupTemp2);
    
    
    MPI_Group_union(groupTemp1,groupTemp2,&client_group_);
    MPI_Comm_create(MPI_COMM_WORLD,client_group_ , &client_comm_);
    
    MPI_Group_rank(client_group_, &rank);
    if(rank==0)isClientRoot_ = true ;
    else isClientRoot_ = false ;
}

void IOserver::stop()
{
    if(isRoot_)
    {
        int send=0;
        MPI_Ssend(&send,1,MPI::INT,0,SERVER_CONTROL_TAG,sync_global_comm_);
        //cout<< "call stop"<<endl;
    }
    
    MPI_Barrier(compute_world_comm_);
}
int IOserver::openOstream()
{
    MPI_Status status;
    int flag=1;
    
    if(isRoot_)
    {
        bool getState=true;
        while(getState)
        {
            MPI_Iprobe(0,SERVER_STATE_TAG,sync_global_comm_,&flag,&status);
            if(flag==true)
            {
                MPI_Recv(&state_,1,MPI::INT,0,SERVER_STATE_TAG,sync_global_comm_,&status);
            }
            else  getState=false;
        }
        
        if(state_ == SERVER_READY)
        {
            int send=CONTROL_OPEN_OSTREAM ;
            MPI_Bsend(&send,1,MPI::INT,0,SERVER_CONTROL_TAG,sync_global_comm_);
            //cout << "call openostream"<<endl;
        }
    }
    
    MPI_Bcast(&state_,1,MPI::INT,0,compute_world_comm_);
    
    
    
    if(state_ == SERVER_READY) return OSTREAM_SUCCESS;
    else return OSTREAM_FAIL;
    
}

void IOserver::closeOstream()
{
    int send = CONTROL_CLOSE_OSTREAM;
    MPI_Send(&send,1,MPI::INT,client_size_,CLOSE_OSTREAM_TAG,client_comm_);
}

ioserver_file IOserver::openFile(string file_name,
                                 int type,
                                 hid_t mem_datatype,
                                 hid_t file_datatype,
                                 LATfield2::Lattice * lat,
                                 int components,
                                 int array_size
                                 )
{

    
    
    //verification:
    
    if(type ==UNSTRUCTURED_H5_FILE)
    {
        if(mem_datatype == -1|| file_datatype == -1 )
        {
            if(isRoot_)cout<< "IOserver::openFile: UNSTRUCTURED_H5_FILE have to have their datatype set!"<<endl;
            if(isRoot_)cout<< "Aborting" << endl;
            MPI_Abort(MPI_COMM_WORLD,0);
        }
    }

    
    if(type ==STRUCTURED_H5_FILE)
    {
            if(mem_datatype == -1|| file_datatype == -1 || lat==NULL)
            {
                if(isRoot_)cout<< "IOserver::openFile: STRUCTURED_H5_FILE have to have their datatype,lattice set!"<<endl;
                if(isRoot_)cout<< "Aborting" << endl;
                MPI_Abort(MPI_COMM_WORLD,0);
            }
    }
    
    
    
    ioserver_file file;
    
    
    MPI_Status status;
    char * send;
    int send_size;
    
    int fileType = type;
    
    int  dim;
    hsize_t * size;
    hsize_t offset;
    int itemp;
    
    string my_filename;
    my_filename = file_name + "_" + int2string(my_node_,999);
    
    int filename_len=my_filename.size();
    
    //cout<<filename_len<<endl;
    if(type==UNSTRUCTURED_BIN_FILE)
    {
        //cout<<"compute: opening unstructured bin file: "<<my_filename<<endl;
        send_size = (2*sizeof(int)) + filename_len;
        send = new char[send_size];
        
        
        
    }
    else if(type == UNSTRUCTURED_H5_FILE)
    {
        //cout<<"compute: opening unstructured bin file: "<<my_filename<<endl;
        size_t mem_datatype_size;
        size_t file_datatype_size;
        char * buf_datatype;
        
        buf_datatype =NULL;
        H5Tencode(mem_datatype,buf_datatype,&mem_datatype_size);
        H5Tencode(file_datatype,buf_datatype,&file_datatype_size);
        
        file.sizeof_mem_dtype = H5Tget_size(mem_datatype);
        
        send_size = (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size;
        send = new char[send_size];
        
        memcpy(send + ( (2*sizeof(int))+ filename_len ) ,(char*)&mem_datatype_size,sizeof(hsize_t));
        memcpy(send + ( (2*sizeof(int))+ filename_len ) + sizeof(hsize_t) ,(char*)&file_datatype_size,sizeof(hsize_t));
        
        
        buf_datatype = new char[mem_datatype_size];
        H5Tencode(mem_datatype,buf_datatype,&mem_datatype_size);
        memcpy(send+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) ),buf_datatype,mem_datatype_size);
        delete[] buf_datatype;
        
        
        buf_datatype = new char[file_datatype_size];
        H5Tencode(file_datatype,buf_datatype,&file_datatype_size);
        memcpy(send+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_datatype_size),buf_datatype,file_datatype_size);
        delete[] buf_datatype;

    }
    else if(type ==STRUCTURED_H5_FILE)
    {
        dim = lat->dim();
        size = new hsize_t[dim];
        for(int i=0;i<dim;i++)size[i]=lat->size(i);
        
        int bproc = my_node_ * proc_size1_ / io_node_number_;
        offset=0;
        for(int i=0;i < bproc;i++) offset += lat->sizeLocalAllProcDim1()[i];
        size[dim-2]=0;
        for(int i=0;i < (proc_size1_ / io_node_number_);i++)size[dim-2]+= lat->sizeLocalAllProcDim1()[bproc+i];
        
        
        size_t mem_datatype_size;
        size_t file_datatype_size;
        char * buf_datatype;
        
        buf_datatype =NULL;
        H5Tencode(mem_datatype,buf_datatype,&mem_datatype_size);
        H5Tencode(file_datatype,buf_datatype,&file_datatype_size);
        
        file.sizeof_mem_dtype = H5Tget_size(mem_datatype);
        
        send_size = (5*sizeof(int)) + filename_len + ((2+dim+1)*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size ;
        send = new char[send_size];
        
        memcpy(send + ( (2*sizeof(int))+ filename_len ) ,(char*)&mem_datatype_size,sizeof(hsize_t));
        memcpy(send + ( (2*sizeof(int))+ filename_len ) + sizeof(hsize_t) ,(char*)&file_datatype_size,sizeof(hsize_t));
        
        
        buf_datatype = new char[mem_datatype_size];
        H5Tencode(mem_datatype,buf_datatype,&mem_datatype_size);
        memcpy(send+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) ),buf_datatype,mem_datatype_size);
        delete[] buf_datatype;
        
        
        buf_datatype = new char[file_datatype_size];
        H5Tencode(file_datatype,buf_datatype,&file_datatype_size);
        memcpy(send+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_datatype_size),buf_datatype,file_datatype_size);
        delete[] buf_datatype;
        
        memcpy(send+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size  ),(char*)&dim,sizeof(int));
        memcpy(send+( (3*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size  ),(char*)size,dim*sizeof(hsize_t));
        memcpy(send+( (3*sizeof(int))+ filename_len + ((2+dim)*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size),(char*)&offset,sizeof(hsize_t));
        itemp =components;
        memcpy(send+( (3*sizeof(int))+ filename_len + ((2+dim+1)*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size  ),(char*)&itemp,sizeof(int));
        itemp =array_size;
        memcpy(send+( (4*sizeof(int))+ filename_len + ((2+dim+1)*sizeof(hsize_t)) + mem_datatype_size + file_datatype_size  ),(char*)&itemp,sizeof(int));
        
        file.dim = dim;
        file.components = components;
        file.array_size = array_size;
        
        delete[] size;
    }
    
    
    memcpy(send,(char*)&fileType,sizeof(int));
    memcpy(send+sizeof(int),(char*)&filename_len,sizeof(int));
    for(int i=0;i<filename_len;i++)send[(2*sizeof(int))+i] = my_filename[i];
    
    int fileInfo[2];
    
    if(isClientRoot_)MPI_Bsend(send,send_size,MPI::CHAR,client_size_,FILE_OPEN_TAG,client_comm_);
    
    
    if(isClientFileRoot_)MPI_Recv(fileInfo,2,MPI::INT,client_size_,FILE_OPEN_TAG,client_comm_,&status);
    
    
    //(&fileID,1,MPI::INT,0,io_node_comm_);
    MPI_Bcast(fileInfo,2,MPI::INT,0,compute_file_comm_);
    
    file.ID=fileInfo[0];
    file.type = fileInfo[1];
    file.is_open = true;
    
    
    if(file.type != type)cout<< "ARGARGARG file created but wrong type"<<endl;
    
    delete[] send;
    
    
    MPI_Barrier(compute_world_comm_);
    
    return file;
    
}
void IOserver::closeFile(ioserver_file file)
{
    if(file.is_open == false)
    {
        cout<< "IOServer: file is not open, cant close it!"<<endl;
        return;
    }
    else
    {
        int send[2];
        send[0]= file.ID;
        send[1]= file.type;
        
        MPI_Bsend(send,2,MPI::INT,client_size_,FILE_CLOSE_TAG,client_comm_);
        file.is_open = false;
    }
}

void IOserver::sendData(ioserver_file file, char * message,long size,int sendType = DEFAULT_SEND)
{
    //verify if the file is unstructured...
    if(file.type!=UNSTRUCTURED_BIN_FILE && file.type!=UNSTRUCTURED_H5_FILE)
    {
        cout<< "not a unstructured file! wrong method to send data, more argument needed..."<<endl;
        return;
    }
    
    long sendinfo[3];
    
    sendinfo[0] = file.ID;
    sendinfo[1] = size;
    sendinfo[2] = file.type;
    
    MPI_Bsend(sendinfo,3,MPI::LONG,client_size_,GET_DATA_TAG,client_comm_);
    
    //cout<< "send data to file id : "<<sendinfo[0]<<", size :"<<sendinfo[1]<<endl;
    
    if(sendType == DEFAULT_SEND)MPI_Send(message,size,MPI::CHAR,client_size_,file.ID,client_comm_);
    else if(sendType == BUFFERED_SEND)MPI_Bsend(message,size,MPI::CHAR,client_size_,file.ID,client_comm_);
    //cout<< "send has been sent data to file id : "<<sendinfo[0]<<endl;
}
void IOserver::sendData(ioserver_file file, char * message,hsize_t * size, hsize_t * offset,int sendType = DEFAULT_SEND)
{
    if(file.type!=STRUCTURED_H5_FILE)
    {
        cout<< "not a structured h5file! wrong method to send data, less argument needed..."<<endl;
        return;
    }
    
    long sendinfo[3 + (file.dim*2)];
    long mSize=1;
    
    for(int i = 0;i<file.dim;i++) mSize *= size[i];
    
    
    sendinfo[0]= file.ID;
    sendinfo[1]= mSize * file.sizeof_mem_dtype * file.components * file.array_size;
    //cout<<sendinfo[1]<<endl;
    sendinfo[2]= file.type;
    for(int i = 0;i<file.dim;i++)sendinfo[3+i] = size[i];
    for(int i = 0;i<file.dim;i++)sendinfo[3+file.dim+i] = offset[i];
    
    
     MPI_Bsend(sendinfo,3 + (file.dim*2),MPI::LONG,client_size_,GET_DATA_TAG,client_comm_);
    
    if(sendType == DEFAULT_SEND)MPI_Send(message,sendinfo[1],MPI::CHAR,client_size_,file.ID,client_comm_);
    else if(sendType == BUFFERED_SEND)MPI_Bsend(message,sendinfo[1],MPI::CHAR,client_size_,file.ID,client_comm_);
}


void IOserver::sendHeader(ioserver_file file,char* header,int size)
{
    if(file.type != UNSTRUCTURED_BIN_FILE)
    {
        cout<<"IOserver::sendHeader: trying to add header to a non binary file..."<<endl;
        cout<<"IOserver::sendHeader: returning without adding the header"<<endl;
        return;
    }
    if(isClientFileRoot_)
    {
        char * send;
        int send_size;
        int ID=file.ID;
        int header_size = size;
        
        send_size = size + (2 * sizeof(int));
        send = new char[send_size];
        
        memcpy(send,(char*)&ID,sizeof(int));
        memcpy(send + sizeof(int),(char*)&header_size,sizeof(int));
        memcpy(send + (2*sizeof(int)),header,header_size);
        
        MPI_Bsend(send,send_size,MPI::CHAR,client_size_,GET_HEADER_TAG,client_comm_);
        delete[] send;
    }
    
}


void IOserver::sendATTR(ioserver_file file,string attr_name, char * attr, int size,hid_t dtype)
{
    if(file.type != UNSTRUCTURED_H5_FILE && file.type!=STRUCTURED_H5_FILE)
    {
        cout<<"IOserver::sendATTR: trying to add attribute to a non hdf5 file..."<<endl;
        cout<<"IOserver::sendATTR: returning without adding the attribute"<<endl;
        return;
    }
    if(isClientFileRoot_)
    {
        //cout<<"adding attribute"<<endl;
        char * send;
        char * buf_datatype;
        size_t buf_datatype_size;
        int ID = file.ID;
        int attr_size = size;
        int send_size;
        int attr_name_size = attr_name.size();
        
        buf_datatype =NULL;
        H5Tencode(dtype,buf_datatype,&buf_datatype_size);
        //cout<<buf_datatype_size<<endl;
        buf_datatype = new char[buf_datatype_size];
        H5Tencode(dtype,buf_datatype,&buf_datatype_size);
        
        send_size = buf_datatype_size + sizeof(size_t) + (3*sizeof(int)) + attr_name_size + (attr_size*H5Tget_size(dtype));
        send = new char[send_size];
        
        //cout<<buf_datatype_size<<endl;
        memcpy(send,(char*)&ID,sizeof(int));
        memcpy(send + sizeof(int),(char*)&buf_datatype_size,sizeof(size_t));
        memcpy(send + sizeof(int) + sizeof(size_t),buf_datatype,buf_datatype_size);
        memcpy(send + sizeof(int) + sizeof(size_t)+ buf_datatype_size ,(char*)&attr_size,sizeof(int));
        memcpy(send + (2*sizeof(int)) + sizeof(size_t)+ buf_datatype_size ,(char*)&attr_name_size,sizeof(int));
        strcpy(send + (3*sizeof(int)) + sizeof(size_t)+ buf_datatype_size,attr_name.c_str());
        memcpy(send + (3*sizeof(int)) + sizeof(size_t)+ buf_datatype_size + attr_name_size,attr,attr_size * H5Tget_size(dtype));
        
        MPI_Bsend(send,send_size,MPI::CHAR,client_size_,GET_ATTR_TAG,client_comm_);
        
        delete[] buf_datatype;
        delete[] send;
         //cout<<"sended attribute"<<endl;
    }
}
void IOserver::sendDataset(ioserver_file file,string dset_name, char * dset, hsize_t dim, hsize_t * size, hid_t dtype)
{
    if(file.type != UNSTRUCTURED_H5_FILE && file.type!=STRUCTURED_H5_FILE)
    {
        cout<<"IOserver::sendATTR: trying to add attribute to a non hdf5 file..."<<endl;
        cout<<"IOserver::sendATTR: returning without adding the attribute"<<endl;
        return;
        }
        if(isClientFileRoot_)
        {
            char * send;
            char * buf_datatype;
            size_t buf_datatype_size;
            int ID = file.ID;
            int dtype_size;
            int send_size;
            int dset_size;
            int dset_name_size = dset_name.size();
            
            hsize_t dimT = dim;
            hsize_t sizeT[dim];
            for(int i = 0 ; i<dim;i++)sizeT[i]=size[i];
            
            dset_size = 1;
            for(int i=0;i<dim;i++)dset_size*=size[i];
            dtype_size = H5Tget_size(dtype);
            dset_size *= dtype_size;
            
            buf_datatype =NULL;
            H5Tencode(dtype,buf_datatype,&buf_datatype_size);
            buf_datatype = new char[buf_datatype_size];
            H5Tencode(dtype,buf_datatype,&buf_datatype_size);
            
            send_size= (3*sizeof(int)) + sizeof(size_t) + buf_datatype_size + ((1+dim)*sizeof(hsize_t)) + dset_name_size + dset_size;
            
            send = new char[send_size];
 
            memcpy(send,(char*)&ID,sizeof(int));
            memcpy(send + sizeof(int),(char*)&buf_datatype_size,sizeof(size_t));
            memcpy(send + sizeof(int) + sizeof(size_t),buf_datatype,buf_datatype_size);
            memcpy(send + sizeof(int) + sizeof(size_t) + buf_datatype_size,(char*)&dimT,sizeof(hsize_t));
            memcpy(send + sizeof(int) + sizeof(size_t) + buf_datatype_size + sizeof(hsize_t),(char*)&sizeT,dim * sizeof(hsize_t));
            memcpy(send + sizeof(int) + sizeof(size_t) + buf_datatype_size + ((1+dim)*sizeof(hsize_t)),(char*)&dset_size,sizeof(int));
            memcpy(send + (2*sizeof(int)) + sizeof(size_t) + buf_datatype_size + ((1+dim)*sizeof(hsize_t)),(char*)&dset_name_size,sizeof(int));
            strcpy(send + (3*sizeof(int)) + sizeof(size_t) + buf_datatype_size + ((1+dim)*sizeof(hsize_t)),dset_name.c_str());
            memcpy(send + (3*sizeof(int)) + sizeof(size_t) + buf_datatype_size + ((1+dim)*sizeof(hsize_t)) + dset_name_size ,dset,dset_size);
            
            MPI_Bsend(send,send_size,MPI::CHAR,client_size_,GET_DSET_TAG,client_comm_);
            //MPI_Bsend(dset,dset_size,MPI::CHAR,client_size_,file.ID + OFFSET_DSET_TAG,client_comm_);
            
            delete[] buf_datatype;
            delete[] send;
            
        }
}

void IOserver::start()
{
    MPI_Status status;
    MPI_Request send_statut_request;
    bool serverOn_flag = true;
    int control;
    
    while(serverOn_flag)
    {
        int control;
        
        state_ = SERVER_READY;
        if(isRoot_)
        {
            //MPI_Bsend(&state_,1,MPI::INT,1,SERVER_STATE_TAG,sync_global_comm_);
            MPI_Isend(&state_,1,MPI::INT,1,SERVER_STATE_TAG,sync_global_comm_,&send_statut_request);
        
            MPI_Recv(&control,1,MPI::INT,1,SERVER_CONTROL_TAG,sync_global_comm_,&status);
        }
        MPI_Bcast(&control,1,MPI::INT,0,io_world_comm_);
        //cout<< "control: " <<control<<endl;
        
        if(control==CONTROL_STOP)
        {
            serverOn_flag=false;
        }
        if(control==CONTROL_OPEN_OSTREAM)
        {
            state_ = SERVER_BUSY;
            if(isRoot_)MPI_Isend(&state_,1,MPI::INT,1,SERVER_STATE_TAG,sync_global_comm_,&send_statut_request);
            
            //cout<< "ostream open"<<endl;
            
            this->ostream();
            
            MPI_Barrier(io_world_comm_);
            
            this->write_files();
            MPI_Barrier(io_world_comm_);
        }

    }
}
void IOserver::ostream()
{
    MPI_Status status,status2;
    MPI_Request send_statut_request;
    
    bool ostream_flag=true;
    int message_flag;
    
    int ostream_close_client_flags[client_size_];
    for(int i = 0;i<client_size_;i++)ostream_close_client_flags[i]=0;
    int ostream_close_flag=0;
    int control;
    
    
    int newFile_ID = 0;
    
    usb_files_.clear();
    ush_files_.clear();
    sh_files_.clear();
    filesNumber_ =0;
    wp_data_ = data_;
    
    bool allfiles_closed=true;
    bool flag_continue= false;
    
    list<unstruct_bin_file>::iterator it_usb;
    list<unstruct_h5_file>::iterator it_ush;
    list<struct_h5_file>::iterator it_sh;
    
    while(ostream_flag)
    {
        flag_continue = false;
        
        
        MPI_Iprobe(0,FILE_OPEN_TAG,client_comm_,&message_flag,&status);
        if(message_flag)
        {
            int msize;
            MPI_Get_count(&status,MPI::CHAR,&msize);
            char getdata[msize];
            int fyleType;
            int filename_len;
            int fileID;

            
            MPI_Barrier(io_world_comm_);
            
            MPI_Recv(getdata,msize,MPI::CHAR,0,FILE_OPEN_TAG,client_comm_,&status);
            
            memcpy((char*)&fyleType,getdata,sizeof(int));
            memcpy( (char*)&filename_len,getdata + sizeof(int) ,sizeof(int));
            
            
            
            fileID = newFile_ID;
            newFile_ID++;
            
            MPI_Bcast(&fileID,1,MPI::INT,0,io_node_comm_);
            //cout<< "generated id "<< fileID <<endl;
            
            
            if(fyleType == UNSTRUCTURED_BIN_FILE)
            {
                unstruct_bin_file file;
                file.filename = new char[filename_len+5];
                for(int i=0;i<filename_len;i++)file.filename[i] = getdata[(2*sizeof(int)) + i];
                file.filename[filename_len] = '.';
                file.filename[filename_len+1] = 'd';
                file.filename[filename_len+2] = 'a';
                file.filename[filename_len+3] = 't';
                file.filename[filename_len+4] = '\0';
                file.fnSize = filename_len+5;
                file.ID = fileID;
                file.isClosed = 0;
                file.isClosed_client = new int[client_size_];
                for(int i=0;i<client_size_;i++)file.isClosed_client[i]=0;
                
                file.headerSize = 0;
                file.totalSize = 0;
                file.msg.clear();
                
                
                
                usb_files_.push_back(file);
                filesNumber_ ++;
                
                string str(file.filename,file.fnSize);
                //cout<< "opening usb file : "<< str << ", id: "<< file.ID<<endl;
                
                
            }
            else if(fyleType == UNSTRUCTURED_H5_FILE)
            {
                unstruct_h5_file file;
                file.filename = new char[filename_len+4];
                for(int i=0;i<filename_len;i++)file.filename[i] = getdata[(2*sizeof(int)) + i];
                file.filename[filename_len] = '.';
                file.filename[filename_len+1] = 'h';
                file.filename[filename_len+2] = '5';
                file.filename[filename_len+3] = '\0';
                file.fnSize = filename_len+4;
                file.ID = fileID;
                file.isClosed = 0;
                file.isClosed_client = new int[client_size_];
                for(int i=0;i<client_size_;i++)file.isClosed_client[i]=0;

                
                hid_t mem_datatype_id,file_datatype_id;
                hsize_t mem_size,file_size;
                
                memcpy((char*)&mem_size,getdata+( (2*sizeof(int))+ filename_len ),sizeof(hsize_t));
                memcpy((char*)&file_size,getdata+( (2*sizeof(int))+ filename_len  + sizeof(hsize_t) ),sizeof(hsize_t));
                
                char * buf_datatype;
                buf_datatype = new char[mem_size];
                memcpy((char*)buf_datatype,getdata+( (2*sizeof(int))+ filename_len  + (2*sizeof(hsize_t)) ),mem_size);
                mem_datatype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                
                file.sizeof_mem_dtype = H5Tget_size(mem_datatype_id);
                
                buf_datatype = new char[file_size];
                memcpy((char*)buf_datatype,getdata+( (2*sizeof(int))+ filename_len  + (2*sizeof(hsize_t)) + mem_size),file_size);
                file_datatype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                
                file.mem_dtype = mem_datatype_id;
                file.file_dtype = file_datatype_id;
                
                file.totalSize = 0;
                file.msg.clear();
                
                ush_files_.push_back(file);
                filesNumber_ ++;
                
                string str(file.filename,file.fnSize);
                //cout<< "opening ush file : "<< str << ", id: "<< file.ID<<endl;
                
            }
            else if(fyleType == STRUCTURED_H5_FILE)
            {
                struct_h5_file file;
                file.filename = new char[filename_len+3];
                for(int i=0;i<filename_len;i++)file.filename[i] = getdata[(2*sizeof(int)) + i];
                file.filename[filename_len] = '.';
                file.filename[filename_len+1] = 'h';
                file.filename[filename_len+2] = '5';
                file.filename[filename_len+3] = '\0';
                file.fnSize = filename_len+4;
                file.ID = fileID;
                file.isClosed = 0;
                file.isClosed_client = new int[client_size_];
                for(int i=0;i<client_size_;i++)file.isClosed_client[i]=0;
                
                hid_t mem_datatype_id,file_datatype_id;
                hsize_t mem_size,file_size;
                
                memcpy((char*)&mem_size,getdata+( (2*sizeof(int))+ filename_len ),sizeof(hsize_t));
                memcpy((char*)&file_size,getdata+( (2*sizeof(int))+ filename_len  + sizeof(hsize_t) ),sizeof(hsize_t));
                
                char * buf_datatype;
                buf_datatype = new char[mem_size];
                memcpy((char*)buf_datatype,getdata+( (2*sizeof(int))+ filename_len  + (2*sizeof(hsize_t)) ),mem_size);
                mem_datatype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                
                file.sizeof_mem_dtype = H5Tget_size(mem_datatype_id);
                
                buf_datatype = new char[file_size];
                memcpy((char*)buf_datatype,getdata+( (2*sizeof(int))+ filename_len  + (2*sizeof(hsize_t)) + mem_size),file_size);
                file_datatype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                
                file.mem_dtype = mem_datatype_id;
                file.file_dtype = file_datatype_id;
                
                int itemp;
                memcpy((char*)&itemp,getdata+( (2*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_size + file_size  ),sizeof(int));
                file.dim = itemp;
                file.size = new hsize_t[itemp];
                hsize_t stemp[itemp];
                memcpy((char*)stemp, getdata + ( (3*sizeof(int))+ filename_len + (2*sizeof(hsize_t)) + mem_size + file_size ),file.dim * sizeof(hsize_t));
                for(int i = 0;i<file.dim;i++)file.size[i]=stemp[i];
                memcpy((char*)stemp, getdata +( (3*sizeof(int))+ filename_len + ((2+file.dim)*sizeof(hsize_t)) + mem_size + file_size ),sizeof(hsize_t));
                file.offset = stemp[0];
                memcpy((char*)&itemp,getdata+( (3*sizeof(int))+ filename_len + ((2+file.dim+1)*sizeof(hsize_t)) + mem_size + file_size),sizeof(int));
                file.components = itemp;
                memcpy((char*)&itemp,getdata+( (4*sizeof(int))+ filename_len + ((2+file.dim+1)*sizeof(hsize_t)) + mem_size + file_size   ),sizeof(int));
                file.array_size = itemp;
                
                file.msg.clear();
                
                sh_files_.push_back(file);
                filesNumber_ ++;
                
                string str(file.filename,file.fnSize);
                //cout<< "opening sh file : "<< str << ", id: "<< file.ID<<endl;

            }
            else
            {
                cout<< "unrecognised file type"<<endl;
                exit(-444);
            }
           
            
            int isend[2];
            isend[0]=fileID;
            isend[1]=fyleType;
            //cout<<fileID<<endl;
            
            if(isIONodeRoot_)MPI_Send(isend,2,MPI::INT,0,FILE_OPEN_TAG,client_comm_);
            
            
            allfiles_closed = false;
            
            
            
            //cout<< "IO: opening file: "<< fyleType <<" , "<< filename_len <<endl;
            //create file
        }
        
        MPI_Iprobe(MPI_ANY_SOURCE,FILE_CLOSE_TAG,client_comm_,&message_flag,&status);
        if(message_flag)
        {
            // close file
            
            int iget[2];
            int ID;
            int client = status.MPI_SOURCE;
            
            
            MPI_Recv(&iget,2,MPI::INT,client,FILE_CLOSE_TAG,client_comm_,&status);
            ID = iget[0];
            
            //cout<< "closing file : "<<ID<< ", " << client<< " , "<<  iget[1]<<endl;
            
            if(iget[1] == UNSTRUCTURED_BIN_FILE)
            {
                //cout<< "closing file 2222222 : "<<ID<< ", " << client<<endl;
                for(it_usb=usb_files_.begin(); it_usb != usb_files_.end(); ++it_usb)
                {
                    if((*it_usb).ID==ID)
                    {
                        int sum=0;
                        (*it_usb).isClosed_client[status.MPI_SOURCE] = 1;
                        for(int i=0;i<client_size_;i++)sum+=(*it_usb).isClosed_client[i];
                        if(sum==client_size_)(*it_usb).isClosed = 1;
                        break;
                    }
                }
                
                
            }
            if(iget[1] == UNSTRUCTURED_H5_FILE)
            {
                //cout<< "closing file 2222222 : "<<ID<< ", " << client<<"  : "<< ush_files_.size() <<endl;
                for(it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)
                {
                    if((*it_ush).ID==ID)
                    {
                        int sum=0;
                        (*it_ush).isClosed_client[status.MPI_SOURCE] = 1;
                        for(int i=0;i<client_size_;i++)sum+=(*it_ush).isClosed_client[i];
                        if(sum==client_size_)(*it_ush).isClosed = 1;
                        break;
                    }
                    
                }
            }
            if(iget[1] == STRUCTURED_H5_FILE)
            {
                for(it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)
                {
                    if((*it_sh).ID==ID)
                    {
                        int sum=0;
                        (*it_sh).isClosed_client[status.MPI_SOURCE] = 1;
                        for(int i=0;i<client_size_;i++)sum+=(*it_sh).isClosed_client[i];
                        if(sum==client_size_)(*it_sh).isClosed = 1;
                        break;
                    }
                }
            }
            
            allfiles_closed = true;
            for (it_usb=usb_files_.begin(); it_usb != usb_files_.end(); ++it_usb)if( (*it_usb).isClosed == 0 )allfiles_closed=false;
            for (it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)if( (*it_ush).isClosed == 0 )allfiles_closed=false;
            for (it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)if( (*it_sh).isClosed == 0 )allfiles_closed=false;
            
        }
        
        
        MPI_Iprobe(MPI_ANY_SOURCE,CLOSE_OSTREAM_TAG,client_comm_,&message_flag,&status);
        if(message_flag)
        {
            //close ostream
            MPI_Recv(&control,1,MPI::INT,status.MPI_SOURCE,CLOSE_OSTREAM_TAG,client_comm_,&status2);
            if(control==CONTROL_CLOSE_OSTREAM)
            {
                ostream_close_client_flags[status.MPI_SOURCE] = 1;
                //cout<< "call to close ostream "<< status.MPI_SOURCE <<endl;
            }
        }
        
        MPI_Iprobe(MPI_ANY_SOURCE,GET_DATA_TAG,client_comm_,&message_flag,&status);
        if(message_flag)
        {
            //get data
            // get info data...
            int infosize;
            int source;
            MPI_Get_count(&status,MPI::LONG,&infosize);
            source = status.MPI_SOURCE;
            long info[infosize];
            
            //cout<< infosize<<endl;//
            MPI_Recv(info,infosize,MPI::LONG,source,GET_DATA_TAG,client_comm_,&status);
            
            if(info[2]==UNSTRUCTURED_BIN_FILE)
            {
                unstruct_message message;
                message.core = source;
                message.data = wp_data_;
                
                //cout<< "receiving usb data file:"<<info[0]<<" , size : "<<info[1]<<endl;
                
                MPI_Recv(message.data,info[1],MPI::CHAR,source,info[0],client_comm_,&status);
                wp_data_ += info[1];
                for(it_usb=usb_files_.begin(); it_usb != usb_files_.end(); ++it_usb)
                {
                    if( (*it_usb).ID==info[0] )
                    {
                        //cout<<"adding data to file: "<< info[0]<<endl;
                        message.size = info[1];
                        (*it_usb).totalSize+=message.size;
                        (*it_usb).msg.push_back(message);
                    }
                }
            
            }
            else if(info[2]==UNSTRUCTURED_H5_FILE)
            {
                unstruct_message message;
                message.core = source;
                message.data = wp_data_;
                
                //cout<< "receiving ush data file:"<<info[0]<<" , size : "<<info[1]<<endl;
                
                MPI_Recv(message.data,info[1],MPI::CHAR,source,info[0],client_comm_,&status);
                wp_data_ += info[1];
                
                for(it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)
                {
                    if( (*it_ush).ID==info[0] )
                    {
                        //cout<<"adding data to file: "<< info[0]<<endl;
                        message.size = info[1]/(*it_ush).sizeof_mem_dtype;
                        (*it_ush).totalSize+=message.size;
                        (*it_ush).msg.push_back(message);
                        break;
                    }
                }
            }
            else if(info[2]==STRUCTURED_H5_FILE)
            {
                //cout<<"structured"<<endl;
                struct_message message;
                message.data = wp_data_;
                MPI_Recv(message.data,info[1],MPI::CHAR,source,info[0],client_comm_,&status);
                wp_data_ += info[1];
                
                for(it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)
                {
                    if( (*it_sh).ID==info[0] )
                    {
                        //cout<<"adding data to file: "<< info[0]<<endl;
                        message.size = new hsize_t[(*it_sh).dim];
                        message.offset = new hsize_t[(*it_sh).dim];
                        for(int i = 0;i <(*it_sh).dim;i++)
                        {
                            message.size[i] = info[3 + i];
                            message.offset[i] = info[3 + (*it_sh).dim + i];
                        }
                        
                        
                        (*it_sh).msg.push_back(message);
                        break;
                    }
                }

            
            
            }
            flag_continue = true;
        }
        
        
        if(io_node_rank_==0)
        {
            MPI_Iprobe(MPI_ANY_SOURCE,GET_HEADER_TAG,client_comm_,&message_flag,&status);
            if(message_flag)
            {
                int infosize;
                int source;
                MPI_Get_count(&status,MPI::CHAR,&infosize);
                source = status.MPI_SOURCE;
                int file_id;
                int header_size;
                char info[infosize];
                bool notfinded = true;
                list<unstruct_bin_file>::iterator it_ubin;
                
                
                
                MPI_Recv(info,infosize,MPI::CHAR,source,GET_HEADER_TAG,client_comm_,&status);
                
                memcpy((char*)&file_id,info,sizeof(int));
                memcpy((char*)&header_size,info+sizeof(int),sizeof(int));
                
                if(usb_files_.size()!=0)
                {
                    for(it_ubin=usb_files_.begin(); it_ubin != usb_files_.end(); ++it_ubin)
                    {
                        if((*it_ubin).ID == file_id)
                        {
                            (*it_ubin).headerSize = header_size;
                            (*it_ubin).header = new char[header_size];
                            memcpy((*it_ubin).header,info+(2*sizeof(int)),header_size);
                            notfinded = false;
                        }
                    }
                }
                

                
                if(notfinded)cout<<"argargarga"<<endl;
                
                
                
                flag_continue = true;
                
                
            }
            MPI_Iprobe(MPI_ANY_SOURCE,GET_ATTR_TAG,client_comm_,&message_flag,&status);
            if(message_flag)
            {
                
                h5_attr attr;
                int infosize;
                int source;
                MPI_Get_count(&status,MPI::CHAR,&infosize);
                source = status.MPI_SOURCE;
                char info[infosize];
                size_t buf_dtype_size;
                int file_id;
                hid_t dtype_id;
                int attr_size;
                int attr_name_size;
                list<unstruct_h5_file>::iterator it_ush;
                list<struct_h5_file>::iterator it_sh;
                bool notfinded = true;
                
                MPI_Recv(info,infosize,MPI::CHAR,source,GET_ATTR_TAG,client_comm_,&status);
                
                memcpy((char*)&file_id,info,sizeof(int));
                
                memcpy((char*)&buf_dtype_size,info + sizeof(int) ,sizeof(size_t));
                char * buf_datatype = new char[buf_dtype_size];
                memcpy(buf_datatype,info + sizeof(int) + sizeof(size_t) ,buf_dtype_size);
                dtype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                
                memcpy((char*)&attr_size,info + sizeof(int) + sizeof(size_t)+ buf_dtype_size ,sizeof(int));
                memcpy((char*)&attr_name_size,info + (2*sizeof(int)) + sizeof(size_t)+ buf_dtype_size ,sizeof(int));
                char * attr_name = new char[attr_name_size];
                memcpy(attr_name,info + (3*sizeof(int)) + sizeof(size_t)+ buf_dtype_size ,attr_name_size);

                attr.name.assign(attr_name,attr_name_size);
                attr.dtype = dtype_id;
                attr.size = attr_size;
                attr.attr = new char[attr_size*H5Tget_size(dtype_id)];
                memcpy(attr.attr,info + (3*sizeof(int)) + sizeof(size_t)+ buf_dtype_size + attr_name_size ,attr_size*H5Tget_size(dtype_id));
                
                
                if(ush_files_.size()!=0)
                {
                    for(it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)
                    {
                        if((*it_ush).ID == file_id)
                        {
                            (*it_ush).attr.push_back(attr);
                            notfinded = false;
                        }
                    }
                }
                
                if(sh_files_.size()!=0)
                {
                    for(it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)
                    {
                        if((*it_sh).ID == file_id)
                        {
                            (*it_sh).attr.push_back(attr);
                            notfinded = false;
                        }
                    }
                }
                
                if(notfinded)cout<<"argargarga"<<endl;
                
                delete[] attr_name;
                
                
                flag_continue = true;
            }
            MPI_Iprobe(MPI_ANY_SOURCE,GET_DSET_TAG,client_comm_,&message_flag,&status);
            if(message_flag)
            {
                h5_dset dset;
                
                int infosize;
                int source;
                MPI_Get_count(&status,MPI::CHAR,&infosize);
                source = status.MPI_SOURCE;
                char info[infosize];
                size_t buf_dtype_size;
                int file_id;
                hid_t dtype_id;
                hsize_t dim;
                int dset_size;
                int dset_name_size;
                list<unstruct_h5_file>::iterator it_ush;
                list<struct_h5_file>::iterator it_sh;
                bool notfinded = true;

                
                MPI_Recv(info,infosize,MPI::CHAR,source,GET_DSET_TAG,client_comm_,&status);
                
                
                memcpy((char*)&file_id,info,sizeof(int));
                
                memcpy((char*)&buf_dtype_size,info + sizeof(int) ,sizeof(size_t));
                char * buf_datatype = new char[buf_dtype_size];
                memcpy(buf_datatype,info + sizeof(int) + sizeof(size_t) ,buf_dtype_size);
                dtype_id =  H5Tdecode(buf_datatype);
                delete[] buf_datatype;
                memcpy((char*)&dim,info + sizeof(int) + sizeof(size_t) + buf_dtype_size,sizeof(hsize_t));
                
                dset.dtype = dtype_id;
                dset.dim=dim;
                dset.size = new hsize_t[dim];
                memcpy((char*)dset.size,info + sizeof(int) + sizeof(size_t) + buf_dtype_size + sizeof(hsize_t),dim * sizeof(hsize_t));
                memcpy((char*)&dset_size,info + sizeof(int) + sizeof(size_t) + buf_dtype_size + ((1+dim)*sizeof(hsize_t)),sizeof(int));
                memcpy((char*)&dset_name_size,info + (2*sizeof(int)) + sizeof(size_t) + buf_dtype_size + ((1+dim)*sizeof(hsize_t)),sizeof(int));
                char * dset_name = new char[dset_name_size];
                memcpy(dset_name,info + (3*sizeof(int)) + sizeof(size_t) + buf_dtype_size + ((1+dim)*sizeof(hsize_t)) ,dset_name_size);
                dset.name.assign(dset_name,dset_name_size);
                dset.dset = new char[dset_size];
                memcpy(dset.dset,info + (3*sizeof(int)) + sizeof(size_t) + buf_dtype_size + ((1+dim)*sizeof(hsize_t)) + dset_name_size ,dset_size);
                
                
                if(ush_files_.size()!=0)
                {
                    for(it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)
                    {
                        if((*it_ush).ID == file_id)
                        {
                            (*it_ush).dset.push_back(dset);
                            notfinded = false;
                        }
                    }
                }
                
                if(sh_files_.size()!=0)
                {
                    for(it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)
                    {
                        if((*it_sh).ID == file_id)
                        {
                            (*it_sh).dset.push_back(dset);
                            notfinded = false;
                        }
                    }
                }
                
                if(notfinded)cout<<"argargarga"<<endl;
                
                delete[] dset_name;
                
                flag_continue = true;
            }
        }
        
        
        if(flag_continue)continue;
        
        
        //look if the ostream can be closed
        
        
        
        

        ostream_close_flag = 0;
        for(int i = 0;i<client_size_;i++)ostream_close_flag+=ostream_close_client_flags[i];
        if(ostream_close_flag==client_size_ && allfiles_closed)
        {
            ostream_flag=false;
            //cout<< "closing ostream"<<endl;
            
        }
    
    
    }
    
    
}
void IOserver::write_files()
{
    list<unstruct_bin_file>::iterator it_usb;
    list<unstruct_h5_file>::iterator it_ush;
    list<struct_h5_file>::iterator it_sh;
    list<h5_attr>::iterator it_attr;
    list<h5_dset>::iterator it_dset;
    
    list<unstruct_message>::iterator it_um;
    list<struct_message>::iterator it_sm;
    
    if(usb_files_.size()!=0)
    {
        //cout<< "wrinting usb file"<<endl;
        usb_files_.sort();
        for(it_usb=usb_files_.begin(); it_usb != usb_files_.end(); ++it_usb)
        {
            MPI_File ofile;
            MPI_Status status;
            
            long file_offset=0;
            if(io_node_rank_==0)file_offset += (*it_usb).headerSize;
            long temp_long;
            
            for(int k=0;k<(io_node_size_-1); k++)
            {
                if(io_node_rank_==k)
                {
                    temp_long = file_offset + (*it_usb).totalSize;
                    MPI_Send(&temp_long,1,MPI_LONG, k+1 , 0, io_node_comm_ );
                }
                if(io_node_rank_==k+1)
                {
                    MPI_Recv( &file_offset, 1, MPI_LONG, k, 0, io_node_comm_, &status);
                }
            }
            
            
            
            
            if((*it_usb).totalSize!=0)
            {
                MPI_File_open(io_node_comm_,(*it_usb).filename,MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&ofile);
                
                
                if((*it_usb).headerSize>0 && io_node_rank_==0)
                {
                    //MPI_File_set_view(ofile,file_offset+msg_offset,MPI_CHAR,MPI_CHAR,(char*)"native",MPI_INFO_NULL);
                    int headerOffset = 0;
                    MPI_File_write_at(ofile, headerOffset,(*it_usb).header,(*it_usb).headerSize, MPI_CHAR, &status);
                }
                (*it_usb).msg.sort();
                long msg_offset=0;
                for(it_um=(*it_usb).msg.begin();it_um != (*it_usb).msg.end();++it_um)
                {
                    MPI_File_set_view(ofile,file_offset+msg_offset,MPI_CHAR,MPI_CHAR,(char*)"native",MPI_INFO_NULL);
                    MPI_File_write_all(ofile,(*it_um).data,(*it_um).size,MPI_CHAR,&status);
                    msg_offset+=(*it_um).size;
                }
                MPI_File_close(&ofile);
            }
        }
    }
    
    
    if(ush_files_.size()!=0)
    {
         //cout<< "wrinting uh5 file"<<endl;
         ush_files_.sort();
         for(it_ush=ush_files_.begin(); it_ush != ush_files_.end(); ++it_ush)
         {
             long fileSize = 0;
             long file_offset = 0;
             long temp_long;
             MPI_Status status;
             
             for(int k=0;k<(io_node_size_-1); k++)
             {
                 if(io_node_rank_==k)
                 {
                     temp_long = file_offset + (*it_ush).totalSize;
                     MPI_Send(&temp_long,1,MPI_LONG, k+1 , 0, io_node_comm_ );
                 }
                 if(io_node_rank_==k+1)
                 {
                     MPI_Recv( &file_offset, 1, MPI_LONG, k, 0, io_node_comm_, &status);
                 }
             }
             
             if(io_node_rank_ == io_node_size_-1)fileSize = file_offset+(*it_ush).totalSize;
             
             MPI_Bcast(&fileSize,1,MPI::LONG,io_node_size_-1,io_node_comm_);
             
             if(fileSize!=0)
             {
                 //write the file!!
                 
                 hid_t plist_id,plist_id_file,file_id,filespace_id,memspace_id,dataset_id,attribute_id,dataspace_id;
                 hsize_t fsize = fileSize;
                 hsize_t msize,offsetf,offset;
                 
#ifdef H5_HAVE_PARALLEL
                 
                 MPI_Info info  = MPI_INFO_NULL;
                
                 plist_id_file = H5Pcreate(H5P_FILE_ACCESS);
                 H5Pset_fapl_mpio(plist_id_file,io_node_comm_,info);
                 file_id = H5Fcreate((*it_ush).filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_file);
                 filespace_id = H5Screate_simple(1,&fsize,NULL);
                 dataset_id = H5Dcreate(file_id, "data",(*it_ush).file_dtype , filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
                 H5Dclose(dataset_id);
                 H5Sclose(filespace_id);
                 
                 
                 (*it_ush).msg.sort();
                 offsetf=file_offset;
                 offset=0;
                 for(it_um=(*it_ush).msg.begin();it_um != (*it_ush).msg.end();++it_um)
                 {
                     msize = (*it_um).size;
                     dataset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
                     memspace_id = H5Screate_simple(1,&msize,NULL);
                     filespace_id = H5Dget_space(dataset_id);
                     H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&msize, NULL);
                     plist_id = H5Pcreate(H5P_DATASET_XFER);
                     H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
                     H5Dwrite(dataset_id,(*it_ush).mem_dtype, memspace_id, filespace_id, plist_id,(*it_um).data);
                     H5Pclose(plist_id);
                     H5Sclose(memspace_id);
                     H5Sclose(filespace_id);
                     H5Dclose(dataset_id);
                     
                     offsetf+=msize;
                 }
                 if(io_node_rank_==0)
                 {
                     if((*it_ush).attr.size() != 0)
                     {
                         for(it_attr = (*it_ush).attr.begin();it_attr != (*it_ush).attr.end();++it_attr)
                         {
                             if((*it_attr).size == 1)dataspace_id = H5Screate(H5S_SCALAR);
                             else dataspace_id = H5Screate_simple(1,&((*it_attr).size),NULL);
                             attribute_id = H5Acreate (file_id,(*it_attr).name.c_str() ,(*it_attr).dtype , dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
                             H5Awrite(attribute_id,(*it_attr).dtype  ,(*it_attr).attr );
                             H5Aclose(attribute_id);
                             H5Sclose(dataspace_id);
                        }
                     }
                     if((*it_ush).dset.size() != 0)
                     {
                         for(it_dset = (*it_ush).dset.begin();it_dset != (*it_ush).dset.end();++it_dset)
                         {
                             hsize_t sizeReverse[(*it_dset).dim];
                             for(int i=0;i<(*it_dset).dim;i++)sizeReverse[i]=(*it_dset).size[(*it_dset).dim-1-i];
                             dataspace_id = H5Screate_simple((*it_dset).dim,sizeReverse,NULL);
                             dataset_id = H5Dcreate(file_id, ((*it_dset).name).c_str(), (*it_dset).dtype, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                             H5Dwrite(dataset_id, (*it_dset).dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,(*it_dset).dset);
                             H5Dclose(dataset_id);
                             H5Sclose(dataspace_id);
                         }
                     }
                 }
                 
                 H5Fclose(file_id);
                 
#else
                 if(io_node_rank_==0)
                 {
                     plist_id = H5Pcreate(H5P_FILE_ACCESS);
                     file_id = H5Fcreate((*it_ush).filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                     H5Pclose(plist_id);
                     filespace_id = H5Screate_simple(1,&fsize,NULL);
                     dataset_id = H5Dcreate(file_id, "data",(*it_ush).file_dtype , filespace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
                     H5Dclose(dataset_id);
                     H5Sclose(filespace_id);
                     
                     
                     if((*it_ush).attr.size() != 0)
                     {
                         for(it_attr = (*it_ush).attr.begin();it_attr != (*it_ush).attr.end();++it_attr)
                         {
                             if((*it_attr).size == 1)dataspace_id = H5Screate(H5S_SCALAR);
                             else dataspace_id = H5Screate_simple(1,&((*it_attr).size),NULL);
                             attribute_id = H5Acreate (file_id,(*it_attr).name.c_str() ,(*it_attr).dtype , dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
                             H5Awrite(attribute_id,(*it_attr).dtype  ,(*it_attr).attr );
                             H5Aclose(attribute_id);
                             H5Sclose(dataspace_id);
                         }
                     }
                     
                     if((*it_ush).dset.size() != 0)
                     {
                         for(it_dset = (*it_ush).dset.begin();it_dset != (*it_ush).dset.end();++it_dset)
                         {
                             hsize_t sizeReverse[(*it_dset).dim];
                             for(int i=0;i<(*it_dset).dim;i++)sizeReverse[i]=(*it_dset).size[(*it_dset).dim-1-i];
                             dataspace_id = H5Screate_simple((*it_dset).dim,sizeReverse,NULL);
                             dataset_id = H5Dcreate(file_id, (*it_dset).name.c_str(), (*it_dset).dtype, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                             H5Dwrite(dataset_id, (*it_dset).dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,(*it_dset).dset);
                             H5Dclose(dataset_id);
                             H5Sclose(dataspace_id);
                         }
                     }
                     
                     /*
                      dataspace_id = H5Screate_simple(2,numPartsSize,NULL);
                      dataset_id = H5Dcreate(file_id, "localBoxOffset", REAL_TYPE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                      H5Dwrite(dataset_id, REAL_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
                      H5Dclose(dataset_id);
                      H5Sclose(dataspace_id);
                      */

                     
                     H5Fclose(file_id);
                 }
                 MPI_Barrier(io_node_comm_);
                 for(int p=0;p<io_node_size_;p++)
                 {
                     MPI_Barrier(io_node_comm_);
                     if(io_node_rank_==p)
                     {
                         plist_id = H5Pcreate(H5P_FILE_ACCESS);
                         file_id = H5Fopen((*it_ush).filename,H5F_ACC_RDWR,plist_id);
                         H5Pclose(plist_id);
                         (*it_ush).msg.sort();
                         offsetf=file_offset;
                         offset=0;
                         for(it_um=(*it_ush).msg.begin();it_um != (*it_ush).msg.end();++it_um)
                         {
                             msize = (*it_um).size;
                             dataset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
                             memspace_id = H5Screate_simple(1,&msize,NULL);
                             filespace_id = H5Dget_space(dataset_id);
                             H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offsetf, NULL,&msize, NULL);
                             plist_id = H5Pcreate(H5P_DATASET_XFER);
                             H5Dwrite(dataset_id,(*it_ush).mem_dtype, memspace_id, filespace_id, plist_id,(*it_um).data);
                             H5Pclose(plist_id);
                             H5Sclose(memspace_id);
                             H5Sclose(filespace_id);
                             H5Dclose(dataset_id);
                             
                             offsetf+=msize;
                         }
                         H5Fclose(file_id);
                     }
                 }
#endif
                 
             }
             
         }
        
    }
    
    if(sh_files_.size()!=0)
    {
        //cout<< "wrinting sh5 file"<<endl;
        sh_files_.sort();
        for(it_sh=sh_files_.begin(); it_sh != sh_files_.end(); ++it_sh)
        {
            hid_t plist_id,plist_id_file,file_id,filespace,dataset_id,dataspace_id,attribute_id,dset_id,memspace;
            hid_t mem_dtype_id,file_dtype_id,dtbase_id;
            hsize_t * components;
            herr_t status;
            
            
            //cout<<(*it_sh).components<<"  "<<(*it_sh).array_size<<endl;
            ///////////////////////////////
            // creat mem datatype
            ///////////////////////////////
            if((*it_sh).components == 1 && (*it_sh).array_size ==1)
            {
                mem_dtype_id = H5Tcopy((*it_sh).mem_dtype);
                status = H5Tset_order(mem_dtype_id, DATA_ORDER);
                components = new hsize_t[1]; //to be sure is allocated when freed
            }
            if((*it_sh).components == 1 && (*it_sh).array_size !=1)
            {
                components = new hsize_t[1];
                components[0] = (*it_sh).array_size;
                dtbase_id = H5Tcopy((*it_sh).mem_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                mem_dtype_id = H5Tarray_create(dtbase_id,1,components);
            }
            if((*it_sh).components != 1 && (*it_sh).array_size ==1)
            {
                components = new hsize_t[1];
                components[0] = (*it_sh).components;
                dtbase_id = H5Tcopy((*it_sh).mem_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                mem_dtype_id = H5Tarray_create(dtbase_id,1,components);
            }
            if((*it_sh).components != 1 && (*it_sh).array_size !=1)
            {
                components = new hsize_t[2];
                components[0] = (*it_sh).array_size;
                components[1] = (*it_sh).components;
                dtbase_id = H5Tcopy((*it_sh).mem_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                mem_dtype_id = H5Tarray_create(dtbase_id,2,components);
            }
            ///////////////////////////////
            ///////////////////////////////

            ///////////////////////////////
            // creat mem datatype
            ///////////////////////////////
            if((*it_sh).components == 1 && (*it_sh).array_size ==1)
            {
                file_dtype_id = H5Tcopy((*it_sh).file_dtype);
                status = H5Tset_order(mem_dtype_id, DATA_ORDER);
                components = new hsize_t[1]; //to be sure is allocated when freed
            }
            if((*it_sh).components == 1 && (*it_sh).array_size !=1)
            {
                components = new hsize_t[1];
                components[0] = (*it_sh).array_size;
                dtbase_id = H5Tcopy((*it_sh).file_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                file_dtype_id = H5Tarray_create(dtbase_id,1,components);
            }
            if((*it_sh).components != 1 && (*it_sh).array_size ==1)
            {
                components = new hsize_t[1];
                components[0] = (*it_sh).components;
                dtbase_id = H5Tcopy((*it_sh).file_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                file_dtype_id = H5Tarray_create(dtbase_id,1,components);
            }
            if((*it_sh).components != 1 && (*it_sh).array_size !=1)
            {
                components = new hsize_t[2];
                components[0] = (*it_sh).array_size;
                components[1] = (*it_sh).components;
                dtbase_id = H5Tcopy((*it_sh).file_dtype);
                status = H5Tset_order(dtbase_id, DATA_ORDER);
                file_dtype_id = H5Tarray_create(dtbase_id,2,components);
            }
            ///////////////////////////////
            ///////////////////////////////
            H5Tclose(dtbase_id);
            
            
#ifdef H5_HAVE_PARALLEL
            
           
            //cout << "arg wroking in parallel" <<endl;
            
            MPI_Info info  = MPI_INFO_NULL;
            
            hsize_t size[(*it_sh).dim];
            for(int i=0;i<(*it_sh).dim;i++)size[i] = (*it_sh).size[(*it_sh).dim-1-i];
            
            
            plist_id_file = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id_file,io_node_comm_,info);
            
            file_id = H5Fcreate((*it_sh).filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_file);
            filespace = H5Screate_simple((*it_sh).dim,size,NULL);
            dataset_id = H5Dcreate(file_id, "data",file_dtype_id , filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
            H5Sclose(filespace);
            H5Dclose(dataset_id);
            
            for(it_sm=(*it_sh).msg.begin();it_sm != (*it_sh).msg.end();++it_sm)
            {
                //offset,offsetf,local_size
                hsize_t offsetf[(*it_sh).dim];
                hsize_t local_size[(*it_sh).dim];
                for(int i =0;i<(*it_sh).dim;i++)
                {
                    offsetf[i]=(*it_sm).offset[(*it_sh).dim-1-i];
                    local_size[i]=(*it_sm).size[(*it_sh).dim-1-i];
                }
                offsetf[1] -= (*it_sh).offset;
                //cout<< offsetf[1] << endl;
                
                dset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
                filespace = H5Dget_space(dset_id);
                memspace = H5Screate_simple((*it_sh).dim,local_size,NULL);
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,local_size, NULL);
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
                H5Dwrite(dset_id, file_dtype_id, memspace, filespace, plist_id, (*it_sm).data);
                H5Pclose(plist_id);
                H5Sclose(filespace);
                H5Dclose(dset_id);
                
            }
            
            if(io_node_rank_==0)
            {
                if((*it_sh).attr.size() != 0)
                {
                    for(it_attr = (*it_sh).attr.begin();it_attr != (*it_sh).attr.end();++it_attr)
                    {
                        if((*it_attr).size == 1)dataspace_id = H5Screate(H5S_SCALAR);
                        else dataspace_id = H5Screate_simple(1,&((*it_attr).size),NULL);
                        attribute_id = H5Acreate (file_id,(*it_attr).name.c_str() ,(*it_attr).dtype , dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
                        H5Awrite(attribute_id,(*it_attr).dtype  ,(*it_attr).attr );
                        H5Aclose(attribute_id);
                        H5Sclose(dataspace_id);
                    }
                }
                
                if((*it_sh).dset.size() != 0)
                {
                    for(it_dset = (*it_sh).dset.begin();it_dset != (*it_sh).dset.end();++it_dset)
                    {
                        hsize_t sizeReverse[(*it_dset).dim];
                        for(int i=0;i<(*it_dset).dim;i++)sizeReverse[i]=(*it_dset).size[(*it_dset).dim-1-i];
                        dataspace_id = H5Screate_simple((*it_dset).dim,sizeReverse,NULL);
                        dataset_id = H5Dcreate(file_id, (*it_dset).name.c_str(), (*it_dset).dtype, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                        H5Dwrite(dataset_id, (*it_dset).dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,(*it_dset).dset);
                        H5Dclose(dataset_id);
                        H5Sclose(dataspace_id);
                    }
                }

            }
            
            H5Fclose(file_id);
            H5Pclose(plist_id_file);

            
            
            
#else
            if(io_node_rank_==0)
            {
                hsize_t size[(*it_sh).dim];
                for(int i=0;i<(*it_sh).dim;i++)size[i] = (*it_sh).size[(*it_sh).dim-1-i];
                
                plist_id = H5Pcreate(H5P_FILE_ACCESS);
                
                file_id = H5Fcreate((*it_sh).filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                filespace = H5Screate_simple((*it_sh).dim,size,NULL);
                dataset_id = H5Dcreate(file_id, "data",file_dtype_id , filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                
                H5Pclose(plist_id);
                H5Sclose(filespace);
                H5Dclose(dataset_id);
                
                
                if((*it_sh).attr.size() != 0)
                {
                    for(it_attr = (*it_sh).attr.begin();it_attr != (*it_sh).attr.end();++it_attr)
                    {
                        if((*it_attr).size == 1)dataspace_id = H5Screate(H5S_SCALAR);
                        else dataspace_id = H5Screate_simple(1,&((*it_attr).size),NULL);
                        attribute_id = H5Acreate (file_id,(*it_attr).name.c_str() ,(*it_attr).dtype , dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
                        H5Awrite(attribute_id,(*it_attr).dtype  ,(*it_attr).attr );
                        H5Aclose(attribute_id);
                        H5Sclose(dataspace_id);
                    }
                }
                
                if((*it_sh).dset.size() != 0)
                {
                    for(it_dset = (*it_sh).dset.begin();it_dset != (*it_sh).dset.end();++it_dset)
                    {
                        hsize_t sizeReverse[(*it_dset).dim];
                        for(int i=0;i<(*it_dset).dim;i++)sizeReverse[i]=(*it_dset).size[(*it_dset).dim-1-i];
                        dataspace_id = H5Screate_simple((*it_dset).dim,sizeReverse,NULL);
                        dataset_id = H5Dcreate(file_id, (*it_dset).name.c_str(), (*it_dset).dtype, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                        H5Dwrite(dataset_id, (*it_dset).dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,(*it_dset).dset);
                        H5Dclose(dataset_id);
                        H5Sclose(dataspace_id);
                    }
                }
                H5Fclose(file_id);

            }
            
            MPI_Barrier(io_node_comm_);
            for(int p=0;p<io_node_size_;p++)
            {
                MPI_Barrier(io_node_comm_);
                if(io_node_rank_==p)
                {
                    
                    plist_id = H5Pcreate(H5P_FILE_ACCESS);
                    file_id = H5Fopen((*it_sh).filename,H5F_ACC_RDWR,plist_id);
                    H5Pclose(plist_id);
                    
                    for(it_sm=(*it_sh).msg.begin();it_sm != (*it_sh).msg.end();++it_sm)
                    {
                        //offset,offsetf,local_size
                        hsize_t offsetf[(*it_sh).dim];
                        hsize_t local_size[(*it_sh).dim];
                        for(int i =0;i<(*it_sh).dim;i++)
                        {
                            offsetf[i]=(*it_sm).offset[(*it_sh).dim-1-i];
                            local_size[i]=(*it_sm).size[(*it_sh).dim-1-i];
                        }
                        offsetf[1] -= (*it_sh).offset;
                        //cout<< offsetf[1] << endl;
  
                        dset_id = H5Dopen(file_id, "/data", H5P_DEFAULT);
                        filespace = H5Dget_space(dset_id);
                        plist_id = H5Pcreate(H5P_DATASET_XFER);
                        memspace = H5Screate_simple((*it_sh).dim,local_size,NULL);
                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL,local_size, NULL);
                        H5Dwrite(dset_id, file_dtype_id, memspace, filespace, plist_id, (*it_sm).data);
                        H5Pclose(plist_id);
                        H5Sclose(filespace);
                        H5Dclose(dset_id);
                        
                    }
                    H5Fclose(file_id);
                }
            }
    
#endif
            H5Tclose(mem_dtype_id);
            H5Tclose(file_dtype_id);
        }
    }
    
    //cout<<"write done"<<endl;
}

#endif
