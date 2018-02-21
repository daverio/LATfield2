#ifndef LATFIELD2_PARALLEL2D_LAYER_DECL_HPP
#define LATFIELD2_PARALLEL2D_LAYER_DECL_HPP

/*! \file LATfield2_parallel2d_decl.hpp
 \brief LATfield2_parallel2d_decl.hpp contains the class Parallel2d declaration.
 \author David Daverio, edited by Wessel Valkenburg
 */


#include <cstdlib>




class Parallel2d_layer{
  public :



  //CONSTRUCTOR AND DESTRUCTOR================
  Parallel2d_layer();
  //Parallel2d_layer(int proc_size0, int proc_size1);

  ~Parallel2d_layer();


  void initialize_topLevel(int lat_world_size, int*  grid_size, int lat_world_rank,int * grid_rank,
                  int root, bool isRoot,bool * last_proc,
  												MPI_Comm lat_world_comm, MPI_Comm dim0_comm, MPI_Comm dim1_comm,
  												MPI_Group lat_world_group, MPI_Group dim0_group, MPI_Group dim1_group);


  void initialize(int lat_world_size, int lat_world_rank,
  								int * grid_size, int * grid_rank,
  								MPI_Comm lat_world_comm_up, MPI_Comm  dim0_comm_up, MPI_Comm dim1_comm_up,
  								MPI_Group lat_world_group_up, MPI_Group dim0_group_up, MPI_Group dim1_group_up,
  								Parallel2d_layer * layer_up, int level,
  								bool r_dim0_flag, bool r_dim1_flag);


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
  MPI_Comm dim0_comm() {return dim0_comm_;}
  /*!
   \return dim1_comm_  MPI_Comm array, array of directional communicator (dim 1, compute processes)
   */
  MPI_Comm dim1_comm() {return dim1_comm_;}
  /*!
   \return dim0_comm_  MPI_Group array, array of directional group (dim 0, compute processes)
   */
  MPI_Group dim0_group() {return dim0_group_;}
  /*!
   \return dim1_comm_  MPI_Group array, array of directional group (dim 1, compute processes)
   */
  MPI_Group dim1_group() {return dim1_group_;}

  bool isPartLayer() {return isPartLayer_;}

  /*!
   \return int, the rank in lat_world_comm_ for a given process in grid.
   \param n : rank in dim0_comm_
   \param m : rank in dim1_comm_
   */
  int grid2world(int n,int m) {return n + grid_size_[0]*m;}

  void abortForce(){MPI_Abort( lat_world_comm_, EXIT_FAILURE);}


private:
  //MEMBER VARIABLES================

  int lat_world_size_;//Number of processes
  int grid_size_[2]; //Number of processes for dim 0 and dim 1
  int lat_world_rank_; //Process ID
  int grid_rank_[2]; // Process ID in the 2d grid

  int root_;
  bool isRoot_;
  bool last_proc_[2];


  MPI_Comm lat_world_comm_, dim0_comm_, dim1_comm_;
  MPI_Group lat_world_group_, dim0_group_, dim1_group_ ;

  int level_;
  Parallel2d_layer * layer_up_;
  Parallel2d_layer * layer_down_;

  bool isPartLayer_;

};




#endif
