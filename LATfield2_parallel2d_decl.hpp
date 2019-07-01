#ifndef LATFIELD2_PARALLEL2D_DECL_HPP
#define LATFIELD2_PARALLEL2D_DECL_HPP

/*! \file LATfield2_parallel2d_decl.hpp
 \brief LATfield2_parallel2d_decl.hpp contains the class Parallel2d declaration.
 \author David Daverio, edited by Wessel Valkenburg
 */


#include <cstdlib>


//=============================================
//MPI PARALLELISM==============================
//=============================================

//#if PARALLEL_MPI   //Define for MPI parallelism

#include "mpi.h"
#define COUT if(parallel.isRoot())cout


/*! \class Parallel2d
 \brief LATfield2d underliying class for paralleization

 The parallel2d class is the handler of the parallelization of LATfield2d.
 LATfield2d distribute n-dimensional lattices into a 2-dimensional cartesian grid of MPI processes, a rod decomposition. The last dimension of the lattice is scattered into the first dimension of the process grid and the last-but-one dimension of the lattice is scattered into the second dimension of the process grid. This choice have been made to increase data locality of the ghost cells (halo), increases the efficiency of method to update them. Due to his scheme of parallelization, LATfield2d is only able to work with lattice of dimension bigger or equal to two.

 The geometry of the process grid (the size of the 2 dimensions), two layers of MPI communicator and simple communication methods are enbended in the "parallel" object, which is an instance of the class Parallel. This object is instantiated but not initialized within the library header, hence the users should never declare an instance of the Parallel class. but rather use directly its pre/defined instance "parallel".

 */

class Parallel2d{
  public :

  //CONSTRUCTOR AND DESTRUCTOR================
  Parallel2d();
  //Parallel2d(int proc_size0, int proc_size1);

  ~Parallel2d();

  //Communicators initialization (grid initialization)===================





  /*!
   Overall LATfield2 initialization when the output server is used. Should be the first call in any LATfield2 based application, as it initialize MPI (preprocessor define: -DEXTERNAL_IO)

   \param proc_size0    : size of the first dimension of the MPI process grid.
   \param proc_size1    : size of the second dimension of the MPI process grid.
   \param IO_total_size : number of MPI process reserved for the IO server.
   \param IO_node_size  : size of 1 goupe of process reserved for the IO server. Each group will write in a seperated file.
   */
  void initialize(int proc_size0, int proc_size1,int IO_total_size, int IO_node_size);

  /*!
   Overall LATfield2 initialization used when the output server is not used. Should be the first call in any LATfield2 based application, as it initialize MPI.

   \param proc_size0 : size of the first dimension of the MPI process grid.
   \param proc_size1 : size of the second dimension of the MPI process grid.
   */
  void initialize(int proc_size0, int proc_size1);

  //ABORT AND BARRIER===============================

  /*!
   Method to kill by force the executable. If one MPI process call this method the executable will be killed.
   */
  void abortForce();
  /*!
   Method to request to kill the executable. This method have to be call by every compute processes. It will wait every process had done the call before killing the executable.
   */
  void abortRequest();
  /*!
   Method to call MPI_Barrier. This barrier is only applied on the compute processes. Every compute process have to perform the call, otherwise the executable will not continue.
   */
  void barrier();

  //GLOBAL and DIRECTIONAL COMPUTE PROCESSES COMMUNICATIONS


  /*!
   Method to broadcast to every compute process a variable. Performed in lat_world_comm_ (compute processes world communicator).
   \param message: variable to send. the receivers will receive the value in that variable.
   \param from: rank (in lat_world_comm_) of the sender.
   */
  template<class Type> void broadcast(Type& message, int from);
  /*!
   Method to broadcast to every compute process a variable array. Performed in lat_world_comm_ (compute processes world communicator).
   \param message : pointer to the array to send. the receivers will receive the value in that variable.
   \param len     : length of the array.
   \param from    : rank (in lat_world_comm_) of the sender.
   */
  template<class Type> void broadcast(Type* array, int len, int from);
  /*!
   Method to perform a directional broadcast of a variable. The processes with grid_rank_[0]==from will broadcast the variable to every process which have same grid_rank_[1].
   \param message : variable to send. the receivers will receive the value in that variable.
   \param from    : grid_rank_[0] of the sender.
   */
  template<class Type> void broadcast_dim0(Type& message, int from);
  /*!
   Method to perform a directional broadcast of a variable. The processes with grid_rank_[0]==from will broadcast the variable to every process which have same grid_rank_[1].
   \param message : pointer to the array to send. the receivers will receive the value in that variable.
   \param len     : length of the array.
   \param from    : grid_rank_[0] of the sender.
   */
  template<class Type> void broadcast_dim0(Type* array, int len, int from);
  /*!
   Method to perform a directional broadcast of a variable. The processes with grid_rank_[1]==from will broadcast the variable to every process which have same grid_rank_[0].
   \param message : variable to send. the receivers will receive the value in that variable.
   \param from    : grid_rank_[1] of the sender.
   */
  template<class Type> void broadcast_dim1(Type& message, int from);
  /*!
   Method to perform a directional broadcast of a variable. The processes with grid_rank_[1]==from will broadcast the variable to every process which have same grid_rank_[0].
   \param message : pointer to the array to send. the receivers will receive the value in that variable.
   \param len     : length of the array.
   \param from    : grid_rank_[1] of the sender.
   */
  template<class Type> void broadcast_dim1(Type* array, int len, int from);




  /*!
   Method to sum a number over all the compute processes. Each process will have the result assigned in the input variable.
   \param number : variable to sum.
   */
  template<class Type> void sum(Type& number);
  /*!
   Method to sum an array of number over all the compute processes. Each process will have the result assigned in the input array.
   \param number : pointer to the array to sum.
   \param len    : size of the array.
   */
  template<class Type> void sum(Type* array, int len);

  /*!
   Method to perform a sum of a number over all the compute processes with same grid_rank_[1]. Each process will have the result assigned in the input variable.
   \param number : variable to sum.
   */
  template<class Type> void sum_dim0(Type& number);
  /*!
   Method to perform a sum of a number over all the compute processes with same grid_rank_[1]. Each process will have the result assigned in the input array.
   \param number : pointer to the array to sum.
   \param len    : size of the array.
   */
  template<class Type> void sum_dim0(Type* array, int len);
  /*!
   Method to perform a sum of a number over all the compute processes with same grid_rank_[0]. Each process will have the result assigned in the input variable.
   \param number: variable to sum.
   */
  template<class Type> void sum_dim1(Type& number);
  /*!
   Method to perform a sum of a number over all the compute processes with same grid_rank_[0]. Each process will have the result assigned in the input array.
   \param number : pointer to the array to sum.
   \param len    : size of the array.
   */
  template<class Type> void sum_dim1(Type* array, int len);
  /////////////////

  template<class Type> void sum_to(Type& number, int dest = 0);
  template<class Type> void sum_to(Type* array, int len, int dest = 0);
  template<class Type> void sum_dim0_to(Type& number, int dest = 0);
  template<class Type> void sum_dim0_to(Type* array, int len, int dest = 0);
  template<class Type> void sum_dim1_to(Type& number, int dest = 0);
  template<class Type> void sum_dim1_to(Type* array, int len, int dest = 0);


  /*!
   Method to find the maximum value of a variable across all the compute processes.
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void max(Type& number);
  /*!
   Method to find the maximum value of an array across all the compute processes.
   \param array : array of numbers to compare, the max value of each element will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void max(Type* array, int len);
  /*!
   Method to find the maximum value of a variable across all the compute processes with the same grid_rank_[1].
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void max_dim0(Type& number);
  /*!
   Method to find the maximum value of a variable across all the compute processes with the same grid_rank_[1].
   \param array : number to compare, the max value will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void max_dim0(Type* array, int len);
  /*!
   Method to find the maximum value of a variable across all the compute processes with the same grid_rank_[0].
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void max_dim1(Type& number);
  /*!
   Method to find the maximum value of a variable across all the compute processes with the same grid_rank_[0].
   \param array : number to compare, the max value will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void max_dim1(Type* array, int len);

  template<class Type> void max_to(Type& number, int dest = 0);
  template<class Type> void max_to(Type* array, int len, int dest = 0);
  template<class Type> void max_dim0_to(Type& number, int dest = 0);
  template<class Type> void max_dim0_to(Type* array, int len, int dest = 0);
  template<class Type> void max_dim1_to(Type& number, int dest = 0);
  template<class Type> void max_dim1_to(Type* array, int len, int dest = 0);

  /*!
   Method to find the minimal value of a variable across all the compute processes.
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void min(Type& number);
  /*!
   Method to find the minimal value of an array across all the compute processes.
   \param array : number to compare, the max value will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void min(Type* array, int len);
  /*!
   Method to find the minimal value of a variable across all the compute processes with the same grid_rank_[1].
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void min_dim0(Type& number);
  /*!
   Method to find the minimal value of a variable across all the compute processes with the same grid_rank_[1].
   \param array : number to compare, the max value will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void min_dim0(Type* array, int len);
  /*!
   Method to find the minimal value of a variable across all the compute processes with the same grid_rank_[0].
   \param number : number to compare, the max value will be assignent to this variable.
   */
  template<class Type> void min_dim1(Type& number);
  /*!
   Method to find the maximum value of a variable across all the compute processes with the same grid_rank_[0].
   \param array : number to compare, the max value will be assignent to this variable.
   \param len   : size of the array.
   */
  template<class Type> void min_dim1(Type* array, int len);

  template<class Type> void min_to(Type& number, int dest = 0);
  template<class Type> void min_to(Type* array, int len, int dest = 0);
  template<class Type> void min_dim0_to(Type& number, int dest = 0);
  template<class Type> void min_dim0_to(Type* array, int len, int dest = 0);
  template<class Type> void min_dim1_to(Type& number, int dest = 0);
  template<class Type> void min_dim1_to(Type* array, int len, int dest = 0);



  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the lat_world_comm communicator.
   \param message : variable to send.
   \param to      : rank of the receiver. (in lat_world_comm)
   */
  template<class Type> void send(Type& message, int to);
  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the lat_world_comm communicator.
   \param array : variable to send.
   \param len   : size of the array.
   \param to    : rank of the receiver. (in lat_world_comm)
   */
  template<class Type> void send(Type* array, int len, int to);
  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the directional communicator associated with the process caller. (direction=0)
   \param message : variable to send.
   \param to      : rank of the receiver. (grid_rank_[0])
   */
  template<class Type> void send_dim0(Type& message, int to);
  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the directional communicator associated with the process caller. (direction=0)
   \param array : variable to send.
   \param len   : size of the array.
   \param to    : rank of the receiver. (grid_rank_[0])
   */
  template<class Type> void send_dim0(Type* array, int len, int to);
  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the directional communicator associated with the process caller. (direction=1)
   \param message : variable to send.
   \param to      : rank of the receiver. (grid_rank_[1])
   */
  template<class Type> void send_dim1(Type& message, int to);
  /*!
   MPI send method on the compute processes. The method calls MPI_Send in the directional communicator associated with the process caller. (direction=1)
   \param array : variable to send.
   \param len   : size of the array.
   \param to    : rank of the receiver. (grid_rank_[1])
   */
  template<class Type> void send_dim1(Type* array, int len, int to);


  /*!
   MPI receive method on the compute processes. The method calls MPI_Recv in the lat_world_comm communicator.
   \param message : variable which will be assigned to the receive message.
   \param from    : rank of the sender. (in lat_world_comm_)
   */
  template<class Type> void receive(Type& message, int from);
  /*!
   MPI receive method on the compute processes. The method call MPI_Recv in the lat_world_comm communicator.
   \param message : variable which will be assigned to the receive message.
   \param len     : size of the array to be received.
   \param from    : rank of the sender. (in lat_world_comm_)
   */
  template<class Type> void receive(Type* array, int len, int from);
  /*!
   MPI receive method on the compute processes. The method call MPI_Recv in the directional communicator associated with the process caller. (direction=0)
   \param message : variable which will be assigned to the receive message.
   \param from    : rank of the sender. (grid_rank_[0])
   */
  template<class Type> void receive_dim0(Type& message, int from);
  /*!
   MPI receive method on the compute processes. The method call MPI_Recv in the directional communicator associated with the process caller. (direction=0)
   \param message : variable which will be assigned to the receive message.
   \param len     : size of the array to be received.
   \param from    : rank of the sender. (grid_rank_[0])
   */
  template<class Type> void receive_dim0(Type* array, int len, int from);
  /*!
   MPI receive method on the compute processes. The method call MPI_Recv in the directional communicator associated with the process caller. (direction=1)
   \param message : variable which will be assigned to the receive message.
   \param from    : rank of the sender. (grid_rank_[1])
   */
  template<class Type> void receive_dim1(Type& message, int from);
  /*!
   MPI receive method on the compute processes. The method call MPI_Recv in the directional communicator associated with the process caller. (direction=1)
   \param message : variable which will be assigned to the receive message.
   \param len     : size of the array to be received.
   \param from    : rank of the sender. (grid_rank_[1])
   */
  template<class Type> void receive_dim1(Type* array, int len, int from);

  /*!
   Method to send a message through dim0 of the process grid. Processes of grid_rank_[0]=N will send the message to the grid_rank_[0]=N+1, with a torus topology. Therefore each process will send and receive data.
   \param bufferSend : pointer to the data which will be sent.
   \param bufferRec  : pointer to the array where the receive data will be assigned.
   \param len        : size of the array bufferSend.
   */
  template<class Type> void sendUp_dim0(Type& bufferSend,Type& bufferRec, long len);
  /*!
   Method to send a message through dim0 of the processes grid. Processes of grid_rank_[0]=N will send the message to the grid_rank_[0]=N-1, with a torus topology. Therefore each process will send and receive data.
   \param bufferSend : pointer to the data which will be sent.
   \param bufferRec  : pointer to the array where the receive data will be assigned.
   \param len        : size of the array bufferSend.
   */
  template<class Type> void sendDown_dim0(Type& bufferSend,Type& bufferRec, long len);
  /*!
   Method to send 2 message through dim0 of the processes grid. Processes of grid_rank_[0]=N will send the bufferSendUp to the grid_rank_[0]=N+1, and the bufferSendDown to the grid_rank_[0]=N-1, with a torus topology. Therefore each process will send and receive 2 message.
   \param bufferSendUp   : pointer to the data which will be sent up.
   \param bufferRecUp    : pointer to the array where the receive down data will be assigned.
   \param lenUp          : size of the array bufferSendUp.
   \param bufferSendDown : pointer to the data which will be sent down.
   \param bufferRecDown  : pointer to the array where the receive up data will be assigned.
   \param lenDown        : size of the array bufferSendUp.
   */
  template<class Type> void sendUpDown_dim0(Type& bufferSendUp,Type& bufferRecUp, long lenUp, Type& bufferSendDown,Type& bufferRecDown, long lenDown );
  /*!
   Method to send a message through dim1 of the processes grid. Processes of grid_rank_[1]=N will send the message to the grid_rank_[1]=N+1, with a torus topology. Therefore each process will send and receive data.
   \param bufferSend : pointer to the data which will be sent.
   \param bufferRec  : pointer to the array where the receive data will be assigned.
   \param len        : size of the array bufferSend.
   */
  template<class Type> void sendUp_dim1(Type& bufferSend,Type& bufferRec, long len);
  /*!
   Method to send a message through dim1 of the processes grid. Processes of grid_rank_[1]=N will send the message to the grid_rank_[1]=N-1, with a torus topology. Therefore each process will send and receive data.
   \param bufferSend : pointer to the data which will be sent.
   \param bufferRec  : pointer to the array where the receive data will be assigned.
   \param len        : size of the array bufferSend.
   */
  template<class Type> void sendDown_dim1(Type& bufferSend,Type& bufferRec, long len);
  /*!
   Method to send 2 message through  dim1 of the processes grid. Processes of grid_rank_[1]=N will send the bufferSendUp to the grid_rank_[1]=N+1, and the bufferSendDown to the grid_rank_[1]=N-1, with a torus topology. Therefore each process will send and receive 2 message.
   \param bufferSendUp   : pointer to the data which will be sent up.
   \param bufferRecUp    : pointer to the array where the receive down data will be assigned.
   \param lenUp          : size of the array bufferSendUp.
   \param bufferSendDown : pointer to the data which will be sent down.
   \param bufferRecDown  : pointer to the array where the receive up data will be assigned.
   \param lenDown        : size of the array bufferSendUp.
   */
  template<class Type> void sendUpDown_dim1(Type& bufferSendUp,Type& bufferRecUp, long lenUp, Type& bufferSendDown,Type& bufferRecDown, long lenDown );



  //MISCELLANEOUS===================
  /*!
   \return lat_world_size_  the number of MPI process (compute processes)
   */
  int size() { return lat_world_size_; }
  /*!
   \return lat_world_rank_  rank of this process (in the compute world). This rank is set to -1 for IO processes.
   */
  int rank() { return lat_world_rank_; }
  /*!
   \return world_size_  the number of MPI process (compute + IOserver)
   */
  int world_size(){ return world_size_; }
  /*!
   \return world_rank_  rank of this process (in the world = compute + IOserver)
   */
  int world_rank(){ return world_rank_; }
  /*!
   \return grid_size_  array of size 2. Size of each dimension of the compute processes grid.
   */
  int *grid_size() { return grid_size_; }
  /*!
   \return grid_size_  array of size 2. Rank on each dimension of the compute proceses grid.
   */
  int *grid_rank() { return grid_rank_; }
  /*!
   \return root_  the rank of the process which is the root of the compute processes grid.
   */
  int root() { return root_; }
  /*!
   \return isRoot_  True for the compute root process, false if not the compute root process.
   */
  bool isRoot() { return isRoot_; }
  /*!
   \return last_proc_  array of size 2 containing the rank of the last process in each dimension of the compute processes grid.
   */
  bool * last_proc() {return last_proc_;}
  /*!
   \return lat_world_comm_  MPI_Comm, the communicator which contains all compute processes.
   */
  MPI_Comm lat_world_comm(){return lat_world_comm_;}
  /*!
   \return lat_world_group_  MPI_Group, the group which contains all compute processes.
   */
  MPI_Group lat_world_group(){return lat_world_group_;}
  /*!
   \return dim0_comm_  MPI_Comm array, array of directional communicator (dim 0, compute processes)
   */
  MPI_Comm *dim0_comm() {return dim0_comm_;}
  /*!
   \return dim1_comm_  MPI_Comm array, array of directional communicator (dim 1, compute processes)
   */
  MPI_Comm *dim1_comm() {return dim1_comm_;}
  /*!
   \return dim0_comm_  MPI_Group array, array of directional group (dim 0, compute processes)
   */
  MPI_Group *dim0_group() {return dim0_group_;}
  /*!
   \return dim1_comm_  MPI_Group array, array of directional group (dim 1, compute processes)
   */
  MPI_Group *dim1_group() {return dim1_group_;}

  /*!
   \return int, the rank in lat_world_comm_ for a given process in grid.
   \param n : rank in dim0_comm_
   \param m : rank in dim1_comm_
   */
  int grid2world(int n,int m) {return n + grid_size_[0]*m;}

#ifdef EXTERNAL_IO
  /*!
   \return isIO_  true if the process is reserved to the IO server, false if the process is a compute process.
   */
  bool isIO(){return isIO_;}
  MPI_Comm IOcomm(){return IO_comm_;}
  MPI_Group IOgroup(){return IO_group_;}
#endif

  /* for compatibility with OS X, when this library is linked in e.g. FalconIC */
  void PleaseNeverFinalizeMPI() { neverFinalizeMPI = true; };

private:
  //MEMBER VARIABLES================

  int lat_world_size_;//Number of processes
  int grid_size_[2]; //Number of processes for dim 0 and dim 1
  int lat_world_rank_; //Process ID
  int grid_rank_[2]; // Process ID in the 2d grid

  int root_;
  bool isRoot_;
  bool last_proc_[2];

  int world_rank_;
  int world_size_;

  MPI_Comm world_comm_,lat_world_comm_, *dim0_comm_, *dim1_comm_;
  MPI_Group world_group_,lat_world_group_, *dim0_group_,*dim1_group_ ;


#ifdef EXTERNAL_IO
  MPI_Comm IO_comm_;
  MPI_Group IO_group_;
#endif

#ifdef EXTERNAL_IO
  bool isIO_;
#endif

  /* for compatibility with OS X, when this library is linked in e.g. FalconIC */
  bool neverFinalizeMPI;


};
/* again, this here is not the actual compiled object, just the reference to its existence elsewhere. */
extern Parallel2d parallel;



#endif
