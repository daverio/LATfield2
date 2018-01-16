#ifndef LATFIELD2_MULTIFIELD_HPP
#define LATFIELD2_MULTIFIELD_HPP


template<class FieldType>
class MultiField : public Field<FieldType>
{

public:

  MultiField(): parallel_layer_(0),Field<FieldType>()
  {
    //Field<FieldType>::status_=0;
    //parallel_layer_ = -1;
  }

  void initialize(MultiLAT& lattice, int components);
  void initialize(MultiLAT& lattice, int rows, int cols, int symmetry);
  void initialize(MultiLAT& lattice,int nMatrix ,int rows, int cols, int symmetry);

  void updateHalo();
  void updateHaloComms();

  void saveHDF5(string filename, string dataset_name);
  void saveHDF5(string filename){this->saveHDF5(filename, "data");}

  void loadHDF5(string filename,string dataset_name);
  void loadHDF5(string filename){this->loadHDF5(filename, "data");}


protected:

  int parallel_layer_;

};

template<class FieldType>
void MultiField<FieldType>::initialize(MultiLAT& lattice, int components)
{
  Field<FieldType>::initialize(lattice,components);
  parallel_layer_ =  lattice.parallel_layer();
}

template<class FieldType>
void MultiField<FieldType>::initialize(MultiLAT& lattice, int rows, int cols, int symmetry)
{
  Field<FieldType>::initialize(lattice,rows,cols,symmetry);
  parallel_layer_ =  lattice.parallel_layer();
}

template<class FieldType>
void MultiField<FieldType>::initialize(MultiLAT& lattice,int nMatrix ,int rows, int cols, int symmetry)
{
  Field<FieldType>::initialize(lattice,nMatrix,rows,cols,symmetry);
  parallel_layer_ =  lattice.parallel_layer();
}

template<class FieldType>
void MultiField<FieldType>::updateHalo()
{

  	Site site(*this->lattice_);
  	int copyfrom;

  	int  dim=this->lattice_->dim();
  	int* jump=new int[dim];
  	int* size=new int[dim];

  	for(int i=0; i<dim; i++)
  	{
  		jump[i]=this->lattice_->jump(i);
  		size[i]=this->lattice_->sizeLocal(i);
  	}

  	for(site.haloFirst(); site.haloTest(); site.haloNext())
  	{
  		//Work out where to copy from
  		copyfrom=site.index();
  		for(int i=0; i<dim; i++)
  		{
  			if( site.coordLocal(i)<0 )
  			{
  				copyfrom += jump[i] * size[i];
  			}
  			else if( site.coordLocal(i) >= size[i] )
  			{
  				copyfrom -= jump[i] * size[i];
  			}
  		}
  		//Copy data
  		for(int i=0; i<this->components_; i++)
  		{
  			this->data_[site.index()*this->components_+i] = this->data_[copyfrom*this->components_+i];
  			//memcpy(data_[site.index()*components_+i],data_[copyfrom*components_+i],sizeof_fieldType_);

  		}
  	}

  	delete[] jump;
  	delete[] size;

  	this->updateHaloComms();
}

template<class FieldType>
void MultiField<FieldType>::updateHaloComms()
{

  	int buffer_size0, buffer_size1,temp;
  	int i,j;

  	//Size of buffer : max size between 2 scatered dimension;
  	buffer_size0 = this->lattice_->halo() * this->components_*this->lattice_->jump(this->lattice_->dim()-1);
  	buffer_size1 = this->lattice_->halo() * this->components_*this->lattice_->jump(this->lattice_->dim()-2) *this->lattice_->sizeLocal(this->lattice_->dim()-1);

  	if(buffer_size0>buffer_size1)temp=buffer_size0;
  	else temp=buffer_size1;

  	FieldType* buffer_send = new FieldType[ temp ];
  	FieldType* buffer_rec = new FieldType[ temp ];

  	FieldType* pointer_send_up;
  	FieldType* pointer_send_down;
  	FieldType* pointer_rec_up;
  	FieldType* pointer_rec_down;


  	pointer_send_up = this->data_ + ((this->lattice_->halo()+1)*this->lattice_->jump(this->lattice_->dim()-1) - 2*this->lattice_->halo()*this->lattice_->jump(this->lattice_->dim()-2))*this->components_;
  	pointer_rec_up = this->data_ + ((this->lattice_->halo()+1)*this->lattice_->jump(this->lattice_->dim()-1) - this->lattice_->halo()*this->lattice_->jump(this->lattice_->dim()-2))*this->components_;

  	pointer_send_down = this->data_ + buffer_size0 + this->lattice_->jump(this->lattice_->dim()-2)*this->lattice_->halo()*this->components_   ;
  	pointer_rec_down = this->data_ + buffer_size0;

  	if(parallel.layer(parallel_layer_).grid_rank()[1]%2==0)
  	{


  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];

  			}
  		}

  		if(parallel.layer(parallel_layer_).grid_rank()[1]!=parallel.layer(parallel_layer_).grid_size()[1]-1)
  		{
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]+1);
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]+1);
  		}

  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				pointer_rec_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];
  			}
  		}


  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];
  			}
  		}


  		if(parallel.layer(parallel_layer_).grid_rank()[1] != 0)
  		{
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]-1);
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1,  parallel.layer(parallel_layer_).grid_rank()[1]-1);
  		}
  		else if(parallel.layer(parallel_layer_).grid_size()[1]%2==0)
  		{
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_size()[1]-1);
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1,  parallel.layer(parallel_layer_).grid_size()[1]-1);
  		}

  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				pointer_rec_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];
  			}
  		}


  	}
  	else
  	{
          for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];
  			}
  		}

  		parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]-1);
  		parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]-1);


  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				pointer_rec_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];
  			}
  		}


  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];
  			}
  		}

  		if(parallel.layer(parallel_layer_).grid_rank()[1]!=parallel.layer(parallel_layer_).grid_size()[1]-1)
  		{
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]+1);
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_rank()[1]+1);
  		}
  		else
  		{
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1,0);
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, 0);
  		}

  		for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  		{
  			for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  			{
  				pointer_rec_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];
  			}
  		}


  	}


  	if(parallel.layer(parallel_layer_).grid_size()[1]%2!=0)
  	{
  		if(parallel.layer(parallel_layer_).grid_rank()[1]==0)
  		{
  			for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  			{
  				for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  				{
  					buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];

  				}
  			}

  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, parallel.layer(parallel_layer_).grid_size()[1]-1);
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1,  parallel.layer(parallel_layer_).grid_size()[1]-1);

  			for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  			{
  				for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  				{
  					pointer_rec_down[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];
  				}
  			}

  		}
  		if(parallel.layer(parallel_layer_).grid_rank()[1]==parallel.layer(parallel_layer_).grid_size()[1]-1)
  		{
  			for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  			{
  				for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  				{
  					buffer_send[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)]= pointer_send_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_];
  				}
  			}
  			parallel.layer(parallel_layer_).receive_dim1( buffer_rec, buffer_size1,0);
  			parallel.layer(parallel_layer_).send_dim1( buffer_send, buffer_size1, 0);

  			for(j=0;j<this->lattice_->sizeLocal(this->lattice_->dim()-1);j++)
  			{
  				for(i=0;i<buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1);i++)
  				{
  					pointer_rec_up[i+j*this->lattice_->jump(this->lattice_->dim()-1)*this->components_] = buffer_rec[i+j*buffer_size1/this->lattice_->sizeLocal(this->lattice_->dim()-1)];

  				}
  			}
  		}
  	}



  	pointer_send_up = this->data_ + this->lattice_->sitesLocalGross() * this->components_ - 2*buffer_size0;
  	pointer_rec_up = this->data_ + this->lattice_->sitesLocalGross() * this->components_ - buffer_size0;

  	pointer_send_down = this->data_ + buffer_size0;
  	pointer_rec_down = this->data_;


  	if(parallel.layer(parallel_layer_).grid_rank()[0]%2==0)
  	{

  		if(parallel.layer(parallel_layer_).grid_rank()[0]!=parallel.layer(parallel_layer_).grid_size()[0]-1)
  		{
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_up, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]+1);
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_up, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]+1);
  		}




  		if(parallel.layer(parallel_layer_).grid_rank()[0] != 0)
  		{
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_down, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]-1);
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_down, buffer_size0,  parallel.layer(parallel_layer_).grid_rank()[0]-1);
  		}
  		else if(parallel.layer(parallel_layer_).grid_size()[0]%2==0)
  		{
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_down, buffer_size0, parallel.layer(parallel_layer_).grid_size()[0]-1);
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_down, buffer_size0,  parallel.layer(parallel_layer_).grid_size()[0]-1);
  		}


  	}
  	else
  	{

  		parallel.layer(parallel_layer_).receive_dim0( pointer_rec_down, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]-1);
  		parallel.layer(parallel_layer_).send_dim0( pointer_send_down, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]-1);


  		if(parallel.layer(parallel_layer_).grid_rank()[0]!=parallel.layer(parallel_layer_).grid_size()[0]-1)
  		{
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_up, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]+1);
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_up, buffer_size0, parallel.layer(parallel_layer_).grid_rank()[0]+1);
  		}
  		else
  		{
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_up, buffer_size0,0);
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_up, buffer_size0, 0);
  		}



  	}

  	if(parallel.layer(parallel_layer_).grid_size()[0]%2!=0)
  	{
  		if(parallel.layer(parallel_layer_).grid_rank()[0]==0)
  		{
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_down, buffer_size0, parallel.layer(parallel_layer_).grid_size()[0]-1);
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_down, buffer_size0,  parallel.layer(parallel_layer_).grid_size()[0]-1);
  		}
  		if(parallel.layer(parallel_layer_).grid_rank()[0]==parallel.layer(parallel_layer_).grid_size()[0]-1)
  		{
  			parallel.layer(parallel_layer_).receive_dim0( pointer_rec_up, buffer_size0,0);
  			parallel.layer(parallel_layer_).send_dim0( pointer_send_up, buffer_size0, 0);
  		}
  	}




  	delete[] buffer_send;
  	delete[] buffer_rec;
}

template <class FieldType>
void  MultiField<FieldType>::saveHDF5(string filename, string dataset_name)
{
#ifdef HDF5

  MPI_Barrier(parallel.lat_world_comm());

  if(parallel.layer(parallel_layer_).isPartLayer())
  {

    save_hdf5(this->data_,this->type_id,this->array_size,this->lattice_->coordSkip(),this->lattice_->size(),this->lattice_->sizeLocal(),this->lattice_->halo(),this->lattice_->dim(),this->components_,filename,dataset_name,parallel.layer(parallel_layer_).lat_world_comm());
  }
#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"saving file in binary: "<<filename<<"BIN"<<endl;
    this->write(filename+"BIN");
#endif

  MPI_Barrier(parallel.lat_world_comm());
}


template <class FieldType>
void  MultiField<FieldType>::loadHDF5(string filename, string dataset_name)
{
#ifdef HDF5
	if(parallel.layer(parallel_layer_).isPartLayer())load_hdf5(this->data_,this->lattice_->coordSkip(),this->lattice_->size(),this->lattice_->sizeLocal(),this->lattice_->halo(),this->lattice_->dim(),this->components_,this->filename,this->dataset_name,this->parallel.layer(parallel_layer_).lat_world_comm());
#else
    COUT<<"LATfield2d must be compiled with HDF5 (flag HDF5 turn on!)"<<endl;
    COUT<<"to be able to use hdf5 data format!!!)"<<endl;
    COUT<<"aborting...."<<endl;
#endif
}



#endif
