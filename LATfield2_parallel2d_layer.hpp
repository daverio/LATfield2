#ifndef LATFIELD2_PARALLEL2D_LAYER_HPP
#define LATFIELD2_PARALLEL2D_LAYER_HPP


#include "LATfield2_parallel2d_layer_decl.hpp"


Parallel2d_layer::Parallel2d_layer()
{
	isPartLayer_ = false;
}


Parallel2d_layer::~Parallel2d_layer()
{

}

void Parallel2d_layer::initialize_topLevel(int lat_world_size, int*  grid_size, int lat_world_rank, int * grid_rank,
																					int root,bool isRoot,bool * last_proc,
																					MPI_Comm lat_world_comm, MPI_Comm dim0_comm, MPI_Comm dim1_comm,
																					MPI_Group lat_world_group, MPI_Group dim0_group, MPI_Group dim1_group)
{

	  lat_world_size_=lat_world_size;//Number of processes
	  grid_size_[0]=grid_size[0];
		grid_size_[1]=grid_size[1];  //Number of processes for dim 0 and dim 1
	  lat_world_rank_=lat_world_rank; //Process ID
	  grid_rank_[0]=grid_rank[0];
		grid_rank_[1]=grid_rank[1]; // Process ID in the 2d grid

	  root_=root;
	  isRoot_=isRoot;
	  last_proc_[0]=last_proc[0];
		last_proc_[1]=last_proc[1];

		lat_world_comm_=lat_world_comm;
		dim0_comm_=dim0_comm;
		dim1_comm_=dim1_comm;
	  lat_world_group_=lat_world_group;
		dim0_group_=dim0_group;
		dim1_group_=dim1_group;

	  level_ = 0;
	 	layer_up_ = NULL;
	 	layer_down_ = NULL;

	  isPartLayer_=true;

}

void Parallel2d_layer::initialize(int lat_world_size, int lat_world_rank,
								int * grid_size, int * grid_rank,
								MPI_Comm lat_world_comm_up, MPI_Comm  dim0_comm_up, MPI_Comm dim1_comm_up,
								MPI_Group lat_world_group_up, MPI_Group dim0_group_up, MPI_Group dim1_group_up,
								Parallel2d_layer * layer_up, int level,
								bool r_dim0_flag, bool r_dim1_flag)
{
		level_ = level;
		layer_up_ = layer_up;
		layer_down_ = NULL;


		//first select the processes to generate the subgrids
		int comm_rank,comm_rank0,comm_rank1;
		int comm_size;
		int ratio_dim0,ratio_dim1;


		if(r_dim0_flag)
		{
			grid_size_[0] = grid_size[0]/2;
			ratio_dim0 = 2;
		}
		else
		{
			grid_size_[0] = grid_size[0];
			ratio_dim0 = 1;
		}

		if(r_dim1_flag)
		{
			grid_size_[1] = grid_size[1]/2;
			ratio_dim1 = 2;
		}
		else
		{
			grid_size_[1] = grid_size[1];
			ratio_dim1 = 1;
		}


		lat_world_size_ =  grid_size_[0]*grid_size_[1];


		int rang[grid_size_[1]][3];
		for(int m=0;m<grid_size_[1];m++)
		{
			rang[m][0]= ratio_dim1*m * grid_size[0];
			rang[m][1]= (ratio_dim1*m+1) * grid_size[0] - ratio_dim0;
			rang[m][2]= ratio_dim0;
		}

		MPI_Group_range_incl(lat_world_group_up,grid_size_[1],rang,&lat_world_group_);
		MPI_Comm_create(lat_world_comm_up, lat_world_group_,&lat_world_comm_);

		rang[0][0]= 0;
		rang[0][1]= grid_size[0] - ratio_dim0;
		rang[0][2]= ratio_dim0;

		MPI_Group_range_incl(dim0_group_up,1,rang,&dim0_group_);
		MPI_Comm_create(dim0_comm_up, dim0_group_,&dim0_comm_);

		rang[0][0]= 0;
		rang[0][1]= grid_size[1] - ratio_dim1;
		rang[0][2]= ratio_dim1;

		MPI_Group_range_incl(dim1_group_up,1,rang,&dim1_group_);
		MPI_Comm_create(dim1_comm_up, dim1_group_,&dim1_comm_);

		MPI_Group_rank(lat_world_group_, &comm_rank);
		if(comm_rank!=MPI_UNDEFINED)
		{
			lat_world_rank_ = comm_rank;
			isPartLayer_ = true;
			root_=0;
			if(lat_world_rank_=root_)isRoot_=true;
			else isRoot_=false;

			MPI_Group_rank(dim0_group_, &comm_rank0);
			MPI_Group_rank(dim1_group_, &comm_rank1);

			if(comm_rank0!=MPI_UNDEFINED)
			{
				grid_rank_[0]=comm_rank0;
				if(grid_rank_[0]==grid_size_[0]-1) last_proc_[0]=true;
				else last_proc_[0]=false;
			}
			else
			{
				cout<<"comm_rank0 not defined in intenal levels, exiting"<<endl;
				exit(0);
			}

			if(comm_rank1!=MPI_UNDEFINED)
			{
				grid_rank_[1]=comm_rank1;
				if(grid_rank_[1]==grid_size_[1]-1) last_proc_[1]=true;
				else last_proc_[1]=false;
			}
			else
			{
				cout<<"comm_rank1 not defined in intenal levels, exiting"<<endl;
				exit(0);
			}


		}
		else
		{
			lat_world_rank_ = -1;
			isPartLayer_ = false;
			root_ =-1;
			isRoot_=false;
			grid_rank_[0]=-1;
			grid_rank_[1]=-1;
			last_proc_[0]=false;
			last_proc_[1]=false;
		}

#ifdef DEBUG_MULTIGRID
	cout<< "Multigrid parallel level: "<< level << " . ranks: ( " << this->grid_rank()[0]<<" , "<<this->grid_rank()[1]<<" )"
			<< " uper layer ranks: ( " << layer_up_->grid_rank()[0]<<" , "<<layer_up_->grid_rank()[1] << " )"<< endl;
#endif


}


template<class Type> void Parallel2d_layer::sum(Type& number)
{
	sum( &number,1 );
}

template<class Type> void Parallel2d_layer::sum(Type* array, int len)
{
	//Gather numbers at root
	Type* gather;
	if( rank() == root() ) gather = new Type[len*size()];
	MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
			   gather, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);

	//Sum on root
	if( isRoot() )
	{
		int i, j;
		for(i=0; i<size(); i++)
		{
			if( i!=root() ) for(j=0; j<len; j++)
			{
				array[j] = array[j] + gather[len*i+j];
			}
		}
	}

	//Broadcast result
	MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);
	// Tidy up (bug found by MDP 12/4/06)
	if(rank() == root() ) delete[] gather;
}

template<class Type> void Parallel2d_layer::sum_dim0(Type& number)
{
	sum_dim0( &number,1 );
}

template<class Type> void Parallel2d_layer::sum_dim0(Type* array, int len)
{
	int i,j,comm_rank;

	//Gather numbers at root
	Type* gather;
	if( grid_rank()[0] == 0 ) gather = new Type[len*grid_size_[0]];

    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,gather, len*sizeof(Type), MPI_BYTE, 0,dim0_comm_);

	//Sum on root
	if( grid_rank()[0] == 0)
	{
		for(i=0; i<grid_size()[0]; i++)
		{
			if( i!=0 ) for(j=0; j<len; j++)
			{
				array[j] = array[j] + gather[len*i+j];
			}
		}
	}

	//Broadcast result

    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0, dim0_comm_);

	// Tidy up (bug found by MDP 12/4/06)
	if( grid_rank()[0] == 0 ) delete[] gather;
}
template<class Type> void Parallel2d_layer::sum_dim1(Type& number)
{
	sum_dim1( &number,1 );
}

template<class Type> void Parallel2d_layer::sum_dim1(Type* array, int len)
{
	int i,j,comm_rank;

	//Gather numbers at root
	Type* gather;
	if( grid_rank_[1] == 0 ) gather = new Type[len*grid_size_[1]];


    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,gather, len*sizeof(Type), MPI_BYTE, 0,dim1_comm_);

	//Sum on root
	if( grid_rank()[1] == 0)
	{
		for(i=0; i<grid_size()[1]; i++)
		{
			if( i!=0 ) for(j=0; j<len; j++)
			{
				array[j] = array[j] + gather[len*i+j];
			}
		}
	}


    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0, dim1_comm_);

	if( grid_rank()[1] == 0 ) delete[] gather;
}

template<class Type> void Parallel2d_layer::max(Type& number)
{
    max( &number,1 );
}

template<class Type> void Parallel2d_layer::max(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( rank() == root() ) gather = new Type[len*size()];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);

    //Find max on root
    if( isRoot() )
    {
        int i, j;
        for(i=0; i<size(); i++)
        {
            if( i!=root() ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] > array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( rank() == root() ) delete[] gather;
}

template<class Type> void Parallel2d_layer::max_dim0(Type& number)
{
    max_dim0( &number,1 );
}
template<class Type> void Parallel2d_layer::max_dim0(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( grid_rank_[0] == 0 ) gather = new Type[len*grid_size_[0]];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, 0,dim0_comm_);

    //Find max on root
    if( grid_rank_[0] == 0  )
    {
        int i, j;
        for(i=0; i<grid_size_[0]; i++)
        {
            if( i!=0 ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] > array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0,dim0_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( grid_rank_[0] == 0  ) delete[] gather;
}

template<class Type> void Parallel2d_layer::max_dim1(Type& number)
{
    max_dim1( &number,1 );
}
template<class Type> void Parallel2d_layer::max_dim1(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( grid_rank_[1] == 0 ) gather = new Type[len*grid_size_[1]];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, 0,dim1_comm_);

    //Find max on root
    if( grid_rank_[1] == 0  )
    {
        int i, j;
        for(i=0; i<grid_size_[1]; i++)
        {
            if( i!=0 ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] > array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0,dim1_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( grid_rank_[1] == 0  ) delete[] gather;
}

template<class Type> void Parallel2d_layer::min(Type& number)
{
    min( &number,1 );
}

template<class Type> void Parallel2d_layer::min(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( rank() == root() ) gather = new Type[len*size()];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);

    //Find min on root
    if( isRoot() )
    {
        int i, j;
        for(i=0; i<size(); i++)
        {
            if( i!=root() ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] < array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, this->root(), lat_world_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( rank() == root() ) delete[] gather;


}

template<class Type> void Parallel2d_layer::min_dim0(Type& number)
{
    min_dim0( &number,1 );
}
template<class Type> void Parallel2d_layer::min_dim0(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( grid_rank_[0] == 0 ) gather = new Type[len*grid_size_[0]];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, 0,dim0_comm_);

    //Find min on root
    if( grid_rank_[0] == 0  )
    {
        int i, j;
        for(i=0; i<grid_size_[0]; i++)
        {
            if( i!=0 ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] < array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0,dim0_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( grid_rank_[0] == 0  ) delete[] gather;
}

template<class Type> void Parallel2d_layer::min_dim1(Type& number)
{
    min_dim1( &number,1 );
}
template<class Type> void Parallel2d_layer::min_dim1(Type* array, int len)
{
    //Gather numbers at root
    Type* gather;
    if( grid_rank_[1] == 0 ) gather = new Type[len*grid_size_[1]];
    MPI_Gather( array, len*sizeof(Type), MPI_BYTE,
               gather, len*sizeof(Type), MPI_BYTE, 0,dim1_comm_);

    //Find min on root
    if( grid_rank_[1] == 0  )
    {
        int i, j;
        for(i=0; i<grid_size_[1]; i++)
        {
            if( i!=0 ) for(j=0; j<len; j++)
            {
                if( gather[len*i+j] < array[j] ) array[j] = gather[len*i+j];
            }
        }
    }

    //Broadcast result
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, 0,dim1_comm_);
    // Tidy up (bug found by MDP 12/4/06)
    if( grid_rank_[1] == 0  ) delete[] gather;
}


template<class Type> void Parallel2d_layer::broadcast(Type& message, int from)
{
	broadcast( &message, 1, from);
};

template<class Type> void Parallel2d_layer::broadcast(Type* array, int len, int from)
{
	MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, from, lat_world_comm_);
}

template<class Type> void Parallel2d_layer::broadcast_dim0(Type& message, int from)
{
	broadcast_dim0( &message, 1, from);
};

template<class Type> void Parallel2d_layer::broadcast_dim0(Type* array, int len, int from)
{
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, from, dim0_comm_);
}

template<class Type> void Parallel2d_layer::broadcast_dim1(Type& message, int from)
{
	broadcast_dim1( &message, 1, from);
};

template<class Type> void Parallel2d_layer::broadcast_dim1(Type* array, int len, int from)
{
    MPI_Bcast( array, len*sizeof(Type), MPI_BYTE, from, dim1_comm_);
}


template<class Type> void Parallel2d_layer::send(Type& message, int to)
{
	MPI_Send( &message, sizeof(Type), MPI_BYTE, to, 0, lat_world_comm_ );
}

template<class Type> void Parallel2d_layer::send(Type* array, int len, int to)
{
	MPI_Send( array, len*sizeof(Type), MPI_BYTE, to, 0, lat_world_comm_ );
}

template<class Type> void Parallel2d_layer::send_dim0(Type& message, int to)
{
    MPI_Send( &message, sizeof(Type), MPI_BYTE, to, 0, dim0_comm_ );
}

template<class Type> void Parallel2d_layer::send_dim0(Type* array, int len, int to)
{
    MPI_Send( array, len*sizeof(Type), MPI_BYTE, to, 0, dim0_comm_ );
}

template<class Type> void Parallel2d_layer::send_dim1(Type& message, int to)
{

    MPI_Send( &message, sizeof(Type), MPI_BYTE, to, 0, dim1_comm_ );
}

template<class Type> void Parallel2d_layer::send_dim1(Type* array, int len, int to)
{

    MPI_Send( array, len*sizeof(Type), MPI_BYTE, to, 0, dim1_comm_ );
}





template<class Type> void Parallel2d_layer::receive(Type& message, int from)
{
	MPI_Status  status;
	MPI_Recv( &message, sizeof(Type), MPI_BYTE, from, 0, lat_world_comm_, &status);
}

template<class Type> void Parallel2d_layer::receive(Type* array, int len, int from)
{
	MPI_Status  status;
	MPI_Recv( array, len*sizeof(Type), MPI_BYTE, from, 0, lat_world_comm_, &status);
}

template<class Type> void Parallel2d_layer::receive_dim0(Type& message, int from)
{

    MPI_Status  status;
    MPI_Recv( &message, sizeof(Type), MPI_BYTE, from, 0, dim0_comm_, &status);
}

template<class Type> void Parallel2d_layer::receive_dim0(Type* array, int len, int from)
{

    MPI_Status  status;
    MPI_Recv( array, len*sizeof(Type), MPI_BYTE, from, 0, dim0_comm_, &status);
}

template<class Type> void Parallel2d_layer::receive_dim1(Type& message, int from)
{

    MPI_Status  status;
    MPI_Recv( &message, sizeof(Type), MPI_BYTE, from, 0,  dim1_comm_, &status);
}

template<class Type> void Parallel2d_layer::receive_dim1(Type* array, int len, int from)
{

    MPI_Status  status;
    MPI_Recv( array, len*sizeof(Type), MPI_BYTE, from, 0, dim1_comm_, &status);
}


template<class Type> void Parallel2d_layer::sendUp_dim0(Type& bufferSend,Type& bufferRec, long len)
{
    if(this->grid_rank()[0]%2==0)
    {
        if(this->grid_rank()[0]!=this->grid_size()[0]-1)// si pas le dernier alors envoie au +1
        {
            this->send_dim0( bufferSend, len , this->grid_rank()[0]+1);
        }

        if(this->grid_rank()[0] != 0)     // si pas le premier alors recois du -1
        {
            this->receive_dim0( bufferRec, len , this->grid_rank()[0]-1);
        }
        else if(this->grid_size()[0]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->receive_dim0( bufferRec, len , this->grid_size()[0]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->receive_dim0( bufferRec, len , this->grid_rank()[0]-1);

        if(this->grid_rank()[0]!=this->grid_size()[0]-1)//si pas dernier alors envoie au +1
        {
            this->send_dim0( bufferSend, len , this->grid_rank()[0]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->send_dim0( bufferSend, len , 0);
        }

    }


    if(this->grid_size()[0]%2!=0)
    {

        if(this->grid_rank()[0]==this->grid_size()[0]-1)//dernier envoie au 0
        {
            this->send_dim0( bufferSend, len , 0);
        }
        if(this->grid_rank()[0]==0)//0 recoie du dernier
        {
            this->receive_dim0( bufferRec, len , this->grid_size()[0]-1);
        }
    }

}
template<class Type> void Parallel2d_layer::sendDown_dim0(Type& bufferSend,Type& bufferRec, long len)
{
    if(this->grid_rank()[0]%2==0)
    {
        if(this->grid_rank()[0]!=this->grid_size()[0]-1)// si pas le dernier alors envoie au +1
        {
            this->receive_dim0( bufferRec, len , this->grid_rank()[0]+1);
        }

        if(this->grid_rank()[0] != 0)     //si pas le premier alors recois du -1
        {
            this->send_dim0( bufferSend, len , this->grid_rank()[0]-1);
        }
        else if(this->grid_size()[0]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->send_dim0( bufferSend, len , this->grid_size()[0]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->send_dim0( bufferSend, len , this->grid_rank()[0]-1);

        if(this->grid_rank()[0]!=this->grid_size()[0]-1)//si pas dernier alors envoie au +1
        {
            this->receive_dim0( bufferRec, len , this->grid_rank()[0]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->receive_dim0( bufferRec, len , 0);
        }

    }


    if(this->grid_size()[0]%2!=0)
    {

        if(this->grid_rank()[0]==this->grid_size()[0]-1)//dernier envoie au 0
        {
            this->receive_dim0( bufferRec, len , 0);
        }
        if(this->grid_rank()[0]==0)//0 recoie du dernier
        {
            this->send_dim0( bufferSend, len , this->grid_size()[0]-1);
        }
    }

}

template<class Type> void Parallel2d_layer::sendUpDown_dim0(Type& bufferSendUp,Type& bufferRecUp, long lenUp, Type& bufferSendDown,Type& bufferRecDown, long lenDown )
{
    if(this->grid_rank()[0]%2==0)
    {
        if(this->grid_rank()[0]!=this->grid_size()[0]-1)// si pas le dernier alors envoie au +1
        {
            this->send_dim0( bufferSendUp, lenUp , this->grid_rank()[0]+1);
            this->receive_dim0( bufferRecDown, lenDown , this->grid_rank()[0]+1);
        }

        if(this->grid_rank()[0] != 0)     // si pas le premier alors recois du -1
        {
            this->receive_dim0( bufferRecUp, lenUp , this->grid_rank()[0]-1);
            this->send_dim0( bufferSendDown, lenDown , this->grid_rank()[0]-1);
        }
        else if(this->grid_size()[0]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->receive_dim0( bufferRecUp, lenUp , this->grid_size()[0]-1);
            this->send_dim0( bufferSendDown, lenDown , this->grid_size()[0]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->receive_dim0( bufferRecUp, lenUp , this->grid_rank()[0]-1);
        this->send_dim0( bufferSendDown, lenDown , this->grid_rank()[0]-1);

        if(this->grid_rank()[0]!=this->grid_size()[0]-1)//si pas dernier alors envoie au +1
        {
            this->send_dim0( bufferSendUp, lenUp , this->grid_rank()[0]+1);
            this->receive_dim0( bufferRecDown, lenDown , this->grid_rank()[0]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->send_dim0( bufferSendUp, lenUp , 0);
            this->receive_dim0( bufferRecDown, lenDown , 0);
        }

    }


    if(this->grid_size()[0]%2!=0)
    {

        if(this->grid_rank()[0]==this->grid_size()[0]-1)//dernier envoie au 0
        {
            this->send_dim0( bufferSendUp, lenUp , 0);
            this->receive_dim0( bufferRecDown, lenDown , 0);
        }
        if(this->grid_rank()[0]==0)//0 recoie du dernier
        {
            this->receive_dim0( bufferRecUp, lenUp , this->grid_size()[0]-1);
            this->send_dim0( bufferSendDown, lenDown , this->grid_size()[0]-1);
        }
    }

}

template<class Type> void Parallel2d_layer::sendUp_dim1(Type& bufferSend,Type& bufferRec, long len)
{
    if(this->grid_rank()[1]%2==0)
    {
        if(this->grid_rank()[1]!=this->grid_size()[1]-1)// si pas le dernier alors envoie au +1
        {
            this->send_dim1( bufferSend, len , this->grid_rank()[1]+1);
        }

        if(this->grid_rank()[1] != 0)     // si pas le premier alors recois du -1
        {
            this->receive_dim1( bufferRec, len , this->grid_rank()[1]-1);
        }
        else if(this->grid_size()[1]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->receive_dim1( bufferRec, len , this->grid_size()[1]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->receive_dim1( bufferRec, len , this->grid_rank()[1]-1);

        if(this->grid_rank()[1]!=this->grid_size()[1]-1)//si pas dernier alors envoie au +1
        {
            this->send_dim1( bufferSend, len , this->grid_rank()[1]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->send_dim1( bufferSend, len , 0);
        }

    }


    if(this->grid_size()[1]%2!=0)
    {

        if(this->grid_rank()[1]==this->grid_size()[1]-1)//dernier envoie au 0
        {
            this->send_dim1( bufferSend, len , 0);
        }
        if(this->grid_rank()[1]==0)//0 recoie du dernier
        {
            this->receive_dim1( bufferRec, len , this->grid_size()[1]-1);
        }
    }

}
template<class Type> void Parallel2d_layer::sendDown_dim1(Type& bufferSend,Type& bufferRec, long len)
{
    if(this->grid_rank()[1]%2==0)
    {
        if(this->grid_rank()[1]!=this->grid_size()[1]-1)// si pas le dernier alors envoie au +1
        {
            this->receive_dim1( bufferRec, len , this->grid_rank()[1]+1);
        }

        if(this->grid_rank()[1] != 0)     // si pas le premier alors recois du -1
        {
            this->send_dim1( bufferSend, len , this->grid_rank()[1]-1);
        }
        else if(this->grid_size()[1]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->send_dim1( bufferSend, len , this->grid_size()[1]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->send_dim1( bufferSend, len , this->grid_rank()[1]-1);

        if(this->grid_rank()[1]!=this->grid_size()[1]-1)//si pas dernier alors envoie au +1
        {
            this->receive_dim1( bufferRec, len , this->grid_rank()[1]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->receive_dim1( bufferRec, len , 0);
        }

    }


    if(this->grid_size()[1]%2!=0)
    {

        if(this->grid_rank()[1]==this->grid_size()[1]-1)//dernier envoie au 0
        {
            this->receive_dim1( bufferRec, len , 0);
        }
        if(this->grid_rank()[1]==0)//0 recoie du dernier
        {
            this->send_dim1( bufferSend, len , this->grid_size()[1]-1);
        }
    }

}

template<class Type> void Parallel2d_layer::sendUpDown_dim1(Type& bufferSendUp,Type& bufferRecUp, long lenUp, Type& bufferSendDown,Type& bufferRecDown, long lenDown )
{
    if(this->grid_rank()[1]%2==0)
    {
        if(this->grid_rank()[1]!=this->grid_size()[1]-1)// si pas le dernier alors envoie au +1
        {
            this->send_dim1( bufferSendUp, lenUp , this->grid_rank()[1]+1);
            this->receive_dim1 (bufferRecDown, lenDown  , this->grid_rank()[1]+1);

        }

        if(this->grid_rank()[1] != 0)     // si pas le premier alors recois du -1
        {
            this->receive_dim1( bufferRecUp, lenUp  , this->grid_rank()[1]-1);
            this->send_dim1( bufferSendDown, lenDown  , this->grid_rank()[1]-1);
        }
        else if(this->grid_size()[1]%2==0)  // si pair et = 0 alors recois du dernier
        {
            this->receive_dim1( bufferRecUp, lenUp  , this->grid_size()[1]-1);
            this->send_dim1( bufferSendDown, lenDown  , this->grid_size()[1]-1);
        }

    }
    else
    {
        //tous recoivent du -1
        this->receive_dim1( bufferRecUp, lenUp  , this->grid_rank()[1]-1);
        this->send_dim1( bufferSendDown, lenDown  , this->grid_rank()[1]-1);

        if(this->grid_rank()[1]!=this->grid_size()[1]-1)//si pas dernier alors envoie au +1
        {
            this->send_dim1( bufferSendUp, lenUp , this->grid_rank()[1]+1);
            this->receive_dim1 (bufferRecDown, lenDown  , this->grid_rank()[1]+1);
        }
        else //pair et dernier, donc enoie au 0
        {
            this->send_dim1( bufferSendUp, lenUp , 0);
            this->receive_dim1 (bufferRecDown, lenDown  , 0);
        }

    }


    if(this->grid_size()[1]%2!=0)
    {

        if(this->grid_rank()[1]==this->grid_size()[1]-1)//dernier envoie au 0
        {
            this->send_dim1( bufferSendUp, lenUp , 0 );
            this->receive_dim1 (bufferRecDown, lenDown  , 0);
        }
        if(this->grid_rank()[1]==0)//0 recoie du dernier
        {
            this->receive_dim1( bufferRecUp, lenUp  , this->grid_size()[1]-1);
            this->send_dim1( bufferSendDown, lenDown  , this->grid_size()[1]-1);
        }
    }
}
#endif
