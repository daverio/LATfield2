#ifndef LATFIELD2_MULTIGRID_HPP
#define LATFIELD2_MULTIGRID_HPP



#include "LATfield2_Multigrid_decl.hpp"


MultiGrid::MultiGrid()
{

}

MultiGrid::~MultiGrid()
{

}

void MultiGrid::initialize(Lattice * lat_top, int levelNumber_max, int minGridperProc)
{
	//creat hierarchy of lattices;

	dim_ = lat_top->dim();
	lat_size_ = new int*[levelNumber_max];
	for(int i=0;i<levelNumber_max;i++)lat_size_[i] = new int[dim_];
	for(int i=0;i<dim_;i++)lat_size_[0][i]=lat_top->size(i);

	//determines nl_;
	bool flag_ok;
	for(int l = 1;l<levelNumber_max;l++)
	{
		flag_ok=true;
		for(int i=0;i<dim_;i++)
		{
			if(lat_size_[l-1][i] % 2 != 0 || lat_size_[l-1][i]/2 < minGridperProc)
			{
				flag_ok=false;
			}
		}

		if(flag_ok==false)
		{
			nl_ = l;
			break;
		}
		else
		{
			for(int i=0;i<dim_;i++)lat_size_[l][i] = lat_size_[l-1][i]/2;
		}
	}
	if(flag_ok==true)nl_=levelNumber_max;

	pLayer_ = new int[nl_];
	plr0_ = new bool[nl_];
	plr1_ = new bool[nl_];
	pLayer_[0] = 0;

	npl_=1;
	int pl_size[2];
	pl_size[0]=parallel.grid_size()[0];
	pl_size[1]=parallel.grid_size()[1];

	//bool plr0_[nl_],plr1_[nl_];

	for(int i=0;i<nl_;i++)
	{
		plr0_[i]=false;
		plr1_[i]=false;
	}

	for(int l = 1;l<nl_;l++)
	{
		if(lat_size_[l][dim_-1] / pl_size[0]  < minGridperProc)
		{
				if( pl_size[0]/2 >=2)
				{
					pl_size[0] /= 2;
					plr0_[npl_]=true;
				}
		}
		if(lat_size_[l][dim_-2] / pl_size[1]  < minGridperProc)
		{
			if( pl_size[1]/2 >=2)
			{
				pl_size[1] /= 2;
				plr1_[npl_]=true;
			}
		}
		if(plr0_[npl_] || plr1_[npl_])
		{
			npl_++;
		}
		pLayer_[l] = npl_-1;
	}



	parallel.MultiGrid_createLayers(nl_,npl_,pLayer_);
	parallel.MultiGrid_initTopLayer();
	for(int i=1;i<npl_;i++)if(  parallel.layers()[i-1].isPartLayer() )parallel.MultiGrid_initLayer(i,plr0_[i],plr1_[i]);

	lattice_ = new MultiLAT[nl_];
	for(int i=0;i<nl_;i++)
	{
		if(parallel.layer(pLayer_[i]).isPartLayer())lattice_[i].initialize(lat_top->dim(),lat_size_[i],lat_top->halo(),pLayer_[i]);
	}


	#ifdef DEBUG_MULTIGRID

		COUT<< "number of lattice layers: "<< nl_<<endl;
		COUT<< "number of parallel layers: "<< npl_<<endl;
		MPI_Barrier(parallel.lat_world_comm());
		COUT<< "list of layers:"<<endl;
		for(int i=0;i<nl_;i++)
		{
			COUT<< "level "<<i<<endl;
			COUT<< "sizes: "<< lat_size_[i][0] << " , " << lat_size_[i][1] << " , "<< lat_size_[i][2] << endl;
			COUT<< "parallel layer: "<< pLayer_[i] << endl;;
			COUT<<"-----------"<<endl;
		}
		MPI_Barrier(parallel.lat_world_comm());
	#endif



}

template<class FieldType>
void MultiGrid::initialize_Field(Field<FieldType> * fieldBase, MultiField<FieldType> *& field)
{

	field = new MultiField<FieldType>[nl_];

	for(int i=0;i<nl_;i++)
	{

			if(parallel.layer(pLayer_[i]).isPartLayer())field[i].initialize(lattice_[i],
																																			fieldBase->nMatrix(),
																																			fieldBase->rows(),
																																			fieldBase->cols(),
																																			fieldBase->symmetry());

			//COUT << "Lattice size: layer["<< i <<"] ("<< lattice_[i].size(0)<<","<< lattice_[i].size(1)<<","<< lattice_[i].size(2)<<");"<<endl;



	}
	for(int i=1;i<nl_;i++)
	{
		if(parallel.layer(pLayer_[i]).isPartLayer())field[i].alloc();
	}

	field[0].data_ = fieldBase->data_;
}


template<class FieldType>
void MultiGrid::restrict(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level, int method)
{
	if(method  == 0)
	{
		if(pLayer_[level]==pLayer_[level+1])restrict3d_spl(src,dst,level);
		else restrict3d_dpl(src,dst,level);
	}
	else
	{
		if(pLayer_[level]==pLayer_[level+1])restrict3d_spl_fw(src,dst,level);
		else restrict3d_dpl_fw(src,dst,level);
	}
}

template<class FieldType>
void MultiGrid::restrict3d_spl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		int lup = level+1;
		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);
		int components = src[level].components();
		for(xc.first();xc.test();xc.next())
		{
			if( xf.setCoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2) )
			{
				for(int c = 0; c< components;c++)dst[lup](xc,c) = src[level](xf,c);
			}
			else
			{
				cout<<"restriction issues...."<<endl;
			}
		}
	}
}

template<class FieldType>
void MultiGrid::restrict3d_spl_fw(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		src[level].updateHalo();
		int lup = level+1;
		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);
		int components = src[level].components();
		for(xc.first();xc.test();xc.next())
		{
			if( xf.setCoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2) )
			{
				for(int c = 0; c< components;c++)dst[lup](xc,c) = src[level](xf,c) * 0.125;

				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0,c) * 0.0625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0,c) * 0.0625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1,c) * 0.0625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1,c) * 0.0625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+2,c) * 0.0625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-2,c) * 0.0625;

				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1+2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1-2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1+2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1-2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+2,c) * 0.03125;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-2,c) * 0.03125;

				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1+2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1-2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1+2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1-2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1+2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1-2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1+2,c) * 0.015625;
				for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1-2,c) * 0.015625;
			}
			else
			{
				cout<<"restriction issues...."<<endl;
			}
		}
	}
}


template<class FieldType>
void MultiGrid::restrict3d_dpl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{

		int lup = level+1;
		int iref,jref,kref;
		int buffer_dimension[3];
		buffer_dimension[0] = lattice_[level].sizeLocal(0)/2;
		buffer_dimension[1] = lattice_[level].sizeLocal(1)/2;
		buffer_dimension[2] = lattice_[level].sizeLocal(2)/2;
		int components = src[level].components();
		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);
		int grid_rec_ranks[2];
		int send_rank;
		FieldType * buffer;
		long bufferSize = buffer_dimension[0]*buffer_dimension[1]*buffer_dimension[2]*components;
		buffer = (FieldType *)malloc(bufferSize * sizeof(FieldType));
		if(parallel.layer(pLayer_[lup]).isPartLayer())
		{
			for(int i = 0;i<buffer_dimension[0];i++)
			{
				for(int j = lattice_[lup].coordSkip()[1];j<lattice_[lup].coordSkip()[1] + buffer_dimension[1];j++)
				{
					for(int k = lattice_[lup].coordSkip()[0]; k<lattice_[lup].coordSkip()[0] + buffer_dimension[2];k++)
					{
						if(xc.setCoord(i,j,k))
						{
							if(xf.setCoord(2*i,2*j,2*k))
							{
								for(int c = 0; c< components;c++)dst[lup](xc,c) = src[level](xf,c);
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 2"<<endl;
							}
						}
						else
						{
							cout<< "MultiGrid::restrict3d_dpl : restriction error 1"<<endl;
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[lup]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
			if(plr1_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[1] = parallel.layer(pLayer_[lup]).grid_rank()[1]*2 + 1;
				grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
								{
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];

								}
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true && plr1_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[lup]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[lup]).grid_rank()[1]*2 + 1;
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
		}
		else
		{
			for(int i = 0;i<buffer_dimension[0];i++)
			{
				for(int j = 0;j<buffer_dimension[1];j++)
				{
					for(int k = 0;k<buffer_dimension[2];k++)
					{
						iref = i*2;
						jref = lattice_[level].coordSkip()[1]+j*2;
						kref = lattice_[level].coordSkip()[0]+k*2;

						if(xf.setCoord(iref,jref,kref))
						{
							for(int c = 0; c<components ; c++)
								buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))] = src[level](xf,c);
						}
						else
						{
							cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true) grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0]  - parallel.layer(pLayer_[level]).grid_rank()[0]%2;
			else grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0];
			if(plr1_[pLayer_[lup]] == true) grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1]  - parallel.layer(pLayer_[level]).grid_rank()[1]%2;
			else grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1];
			send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
			parallel.layer(pLayer_[level]).send(buffer,bufferSize,send_rank);
		}
		free(buffer);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template<class FieldType>
void MultiGrid::restrict3d_dpl_fw(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		src[level].updateHalo();
		int lup = level+1;
		int iref,jref,kref;
		int buffer_dimension[3];
		buffer_dimension[0] = lattice_[level].sizeLocal(0)/2;
		buffer_dimension[1] = lattice_[level].sizeLocal(1)/2;
		buffer_dimension[2] = lattice_[level].sizeLocal(2)/2;
		int components = src[level].components();
		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);
		int grid_rec_ranks[2];
		int send_rank;
		FieldType * buffer;
		long bufferSize = buffer_dimension[0]*buffer_dimension[1]*buffer_dimension[2]*components;
		buffer = (FieldType *)malloc(bufferSize * sizeof(FieldType));
		if(parallel.layer(pLayer_[lup]).isPartLayer())
		{
			for(int i = 0;i<buffer_dimension[0];i++)
			{
				for(int j = lattice_[lup].coordSkip()[1];j<lattice_[lup].coordSkip()[1] + buffer_dimension[1];j++)
				{
					for(int k = lattice_[lup].coordSkip()[0]; k<lattice_[lup].coordSkip()[0] + buffer_dimension[2];k++)
					{
						if(xc.setCoord(i,j,k))
						{
							if(xf.setCoord(2*i,2*j,2*k))
							{
								for(int c = 0; c< components;c++)dst[lup](xc,c) = src[level](xf,c) * 0.125;

								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0,c) * 0.0625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0,c) * 0.0625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1,c) * 0.0625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1,c) * 0.0625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+2,c) * 0.0625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-2,c) * 0.0625;

								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1+2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+1-2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1+2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-1-2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+2,c) * 0.03125;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-2,c) * 0.03125;

								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0+1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf+0-1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0+1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)dst[lup](xc,c) += src[level](xf-0-1-2,c) * 0.015625;
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 2"<<endl;
							}
						}
						else
						{
							cout<< "MultiGrid::restrict3d_dpl : restriction error 1"<<endl;
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[lup]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
			if(plr1_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[1] = parallel.layer(pLayer_[lup]).grid_rank()[1]*2 + 1;
				grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
								{
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];

								}
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true && plr1_[pLayer_[lup]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[lup]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[lup]).grid_rank()[1]*2 + 1;
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0];i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									dst[lup](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
							}
							else
							{
								cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
							}
						}
					}
				}
			}
		}
		else
		{
			for(int i = 0;i<buffer_dimension[0];i++)
			{
				for(int j = 0;j<buffer_dimension[1];j++)
				{
					for(int k = 0;k<buffer_dimension[2];k++)
					{
						iref = i*2;
						jref = lattice_[level].coordSkip()[1]+j*2;
						kref = lattice_[level].coordSkip()[0]+k*2;

						if(xf.setCoord(iref,jref,kref))
						{
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	= src[level](xf,c) * 0.125;

								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0,c) * 0.0625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0,c) * 0.0625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+1,c) * 0.0625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-1,c) * 0.0625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+2,c) * 0.0625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-2,c) * 0.0625;

								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0+1,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0-1,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0+1,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0-1,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+1+2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+1-2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-1+2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-1-2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0+2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0-2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0+2,c) * 0.03125;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0-2,c) * 0.03125;

								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0+1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0+1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0-1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf+0-1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0+1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0+1-2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0-1+2,c) * 0.015625;
								for(int c = 0; c< components;c++)
								 	buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]
																	+= src[level](xf-0-1-2,c) * 0.015625;
						}
						else
						{
							cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
						}
					}
				}
			}
			if(plr0_[pLayer_[lup]] == true) grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0]  - parallel.layer(pLayer_[level]).grid_rank()[0]%2;
			else grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0];
			if(plr1_[pLayer_[lup]] == true) grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1]  - parallel.layer(pLayer_[level]).grid_rank()[1]%2;
			else grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1];
			send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
			parallel.layer(pLayer_[level]).send(buffer,bufferSize,send_rank);
		}
		free(buffer);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template<class FieldType>
void MultiGrid::prolong(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	if(pLayer_[level]==pLayer_[level-1])prolong3d_spl(src,dst,level);
	else prolong3d_dpl(src,dst,level);
}

template<class FieldType>
void MultiGrid::prolong3d_spl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		int ldown = level-1;

		Site xc(lattice_[level]);
		Site xf(lattice_[ldown]);

		int components = src[level].components();

		src[level].updateHalo();

		for(xc.first();xc.test();xc.next())
		{
			if(xf.setCoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2))
			{

				for(int c = 0; c< components;c++)
					dst[ldown](xf,c) = src[level](xc,c);

				for(int i=0;i<3;i++)
					for(int c = 0; c< components;c++)
						dst[ldown](xf+i,c) = (src[level](xc,c) + src[level](xc+i,c))/2.0;


				for(int c = 0; c< components;c++)
						dst[ldown](xf+0+1,c) = (src[level](xc,c) + src[level](xc+0,c)
																			+ src[level](xc+1,c) + src[level](xc+0+1,c))/4.0;

				for(int c = 0; c< components;c++)
						dst[ldown](xf+0+2,c) = (src[level](xc,c) + src[level](xc+0,c)
																			+ src[level](xc+2,c) + src[level](xc+0+2,c))/4.0;

				for(int c = 0; c< components;c++)
							dst[ldown](xf+1+2,c) = (src[level](xc,c) + src[level](xc+1,c)
																				+ src[level](xc+2,c) + src[level](xc+1+2,c))/4.0;


				for(int c = 0; c< components;c++)
					dst[ldown](xf+0+1+2,c) = (src[level](xc,c) + src[level](xc+1,c)
																		+ src[level](xc+2,c) + src[level](xc+1+2,c)
																		+ src[level](xc+0,c) + src[level](xc+0+1,c)
																		+ src[level](xc+0+2,c) + src[level](xc+0+1+2,c))/8.0;

			}
			else
			{
				cout<< "MultiGrid::prolong3d_spl : restriction error 1"<<endl;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template<class FieldType>
void MultiGrid::prolong3d_dpl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(parallel.layer(pLayer_[level-1]).isPartLayer())
	{
		int ldown = level-1;

		int iref,jref,kref;
		int irefc,jrefc,krefc;

		int buffer_dimension[3];
		buffer_dimension[0] = lattice_[ldown].sizeLocal(0)/2+1;
		buffer_dimension[1] = lattice_[ldown].sizeLocal(1)/2+1;
		buffer_dimension[2] = lattice_[ldown].sizeLocal(2)/2+1;

		int components = src[ldown].components();

		Site xf(lattice_[ldown]);


		int grid_rec_ranks[2];
		int send_rank;

		FieldType * buffer;
		FieldType * buffer1;
		FieldType * buffer2;

		long bufferSize = buffer_dimension[0]*buffer_dimension[1]*buffer_dimension[2]*components;

		buffer  = (FieldType *)malloc(bufferSize * sizeof(FieldType));
		buffer1 = (FieldType *)malloc(bufferSize * sizeof(FieldType));
		buffer2 = (FieldType *)malloc(bufferSize * sizeof(FieldType));

		/*
		cout<< "player:" << pLayer_[level]
				<< " ,rank : " << parallel.layer(pLayer_[level]).rank()
				<< " ,grid_rank[0] : " << parallel.layer(pLayer_[level]).grid_rank()[0]
				<< " ,grid_rank[1] : " << parallel.layer(pLayer_[level]).grid_rank()[1]
				<< " ,buffer size: "<<bufferSize << " ,comp: "<<components
				<< " ,buffer_dimension[0] : "<< buffer_dimension[0]
				<< " ,buffer_dimension[1] : "<< buffer_dimension[1]
				<< " ,buffer_dimension[2] : "<< buffer_dimension[2]
				<<endl;
				*/

		if(parallel.layer(pLayer_[level]).isPartLayer())
		{
			src[level].updateHalo();

			Site xc(lattice_[level]);
			//Site xf(lattice_[ldown]);

			for(int i=0;i<buffer_dimension[0]-1;i++)
			{
				for(int j=lattice_[level].coordSkip()[1];j<lattice_[level].coordSkip()[1]+buffer_dimension[1]-1;j++)
				{
					for(int k=lattice_[level].coordSkip()[0];k<lattice_[level].coordSkip()[0]+buffer_dimension[2]-1;k++)
					{
						if(xc.setCoord(i,j,k))
						{
							if(xf.setCoord(i*2,j*2,k*2))
							{
								for(int c = 0; c< components;c++)
									dst[ldown](xf,c) = src[level](xc,c);

								for(int i=0;i<3;i++)
									for(int c = 0; c< components;c++)
										dst[ldown](xf+i,c) = (src[level](xc,c) + src[level](xc+i,c))/2.0;


								for(int c = 0; c< components;c++)
										dst[ldown](xf+0+1,c) = (src[level](xc,c) + src[level](xc+0,c)
																							+ src[level](xc+1,c) + src[level](xc+0+1,c))/4.0;

								for(int c = 0; c< components;c++)
										dst[ldown](xf+0+2,c) = (src[level](xc,c) + src[level](xc+0,c)
																							+ src[level](xc+2,c) + src[level](xc+0+2,c))/4.0;

								for(int c = 0; c< components;c++)
											dst[ldown](xf+1+2,c) = (src[level](xc,c) + src[level](xc+1,c)
																								+ src[level](xc+2,c) + src[level](xc+1+2,c))/4.0;

								for(int c = 0; c< components;c++)
									dst[ldown](xf+0+1+2,c) = (src[level](xc,c) + src[level](xc+1,c)
																						+ src[level](xc+2,c) + src[level](xc+1+2,c)
																						+ src[level](xc+0,c) + src[level](xc+0+1,c)
																						+ src[level](xc+0+2,c) + src[level](xc+0+1+2,c))/8.0;
							}
							else
							{
								cout<< "MultiGrid::prolong3d_dpl : prologation error 2"<<endl;
							}
						}
						else
						{
							cout<< "MultiGrid::prolong3d_dpl : prologation error 1"<<endl;
						}
					}
				}
			}
			if(plr0_[pLayer_[level]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[ldown]).grid_rank()[1];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[ldown]).grid_size()[0] * grid_rec_ranks[1];
				/*cout<< "level up: " << pLayer_[level]
						<< " ,rank up : " << parallel.layer(pLayer_[level]).rank()
						<< " ,grid_rank[0] up : " << parallel.layer(pLayer_[level]).grid_rank()[0]
						<< " ,grid_rank[1] up : " << parallel.layer(pLayer_[level]).grid_rank()[1]
						<< " ,level:" << pLayer_[ldown]
						<< " ,rank : " << parallel.layer(pLayer_[ldown]).rank()
						<< " ,grid_rank[0] : " << parallel.layer(pLayer_[ldown]).grid_rank()[0]
						<< " ,grid_rank[1] : " << parallel.layer(pLayer_[ldown]).grid_rank()[1]
						<< " ,rank from : " << send_rank
						<< " ,grid_rec_ranks[0] : " << grid_rec_ranks[0]
						<< " ,grid_rec_ranks[1] : " << grid_rec_ranks[1]
						<<endl;*/

				for(int i=0;i<buffer_dimension[0];i++)
				{
					for(int j=0;j<buffer_dimension[1];j++)
					{
						for(int k=0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j+lattice_[level].coordSkip()[1];
							kref = k+lattice_[level].coordSkip()[0] + buffer_dimension[2] - 1;
							if(xc.setCoordHalo(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]=src[level](xc,c);
							}
							else
							{
								cout<< "MultiGrid::prolong3d_dpl : prologation error 3"<<endl;
							}
						}
					}
				}
				parallel.layer(pLayer_[ldown]).send(buffer,bufferSize,send_rank);
			}

			if(plr1_[pLayer_[level]] == true)
			{
				grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1]*2 + 1;
				grid_rec_ranks[0] = parallel.layer(pLayer_[ldown]).grid_rank()[0];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[ldown]).grid_size()[0] * grid_rec_ranks[1];
				/*cout<< "level up: " << pLayer_[level]
						<< " ,rank up : " << parallel.layer(pLayer_[level]).rank()
						<< " ,grid_rank[0] up : " << parallel.layer(pLayer_[level]).grid_rank()[0]
						<< " ,grid_rank[1] up : " << parallel.layer(pLayer_[level]).grid_rank()[1]
						<< " ,level:" << pLayer_[ldown]
						<< " ,rank : " << parallel.layer(pLayer_[ldown]).rank()
						<< " ,grid_rank[0] : " << parallel.layer(pLayer_[ldown]).grid_rank()[0]
						<< " ,grid_rank[1] : " << parallel.layer(pLayer_[ldown]).grid_rank()[1]
						<< " ,rank from : " << send_rank
						<< " ,grid_rec_ranks[0] : " << grid_rec_ranks[0]
						<< " ,grid_rec_ranks[1] : " << grid_rec_ranks[1]
						<<endl;*/
				for(int i=0;i<buffer_dimension[0];i++)
				{
					for(int j=0;j<buffer_dimension[1];j++)
					{
						for(int k=0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j+lattice_[level].coordSkip()[1] + buffer_dimension[1] - 1;
							kref = k+lattice_[level].coordSkip()[0];
							if(xc.setCoordHalo(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									buffer1[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]=src[level](xc,c);
							}
							else
							{
								cout<< "MultiGrid::prolong3d_dpl : prologation error 4"<<endl;
							}
						}
					}
				}
				parallel.layer(pLayer_[ldown]).send(buffer1,bufferSize,send_rank);
			}

			if(plr0_[pLayer_[level]] == true && plr1_[pLayer_[level]] == true)
			{
				grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1]*2 + 1;
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[ldown]).grid_size()[0] * grid_rec_ranks[1];
				/*cout<< "level up: " << pLayer_[level]
						<< " ,rank up : " << parallel.layer(pLayer_[level]).rank()
						<< " ,grid_rank[0] up : " << parallel.layer(pLayer_[level]).grid_rank()[0]
						<< " ,grid_rank[1] up : " << parallel.layer(pLayer_[level]).grid_rank()[1]
						<< " ,level:" << pLayer_[ldown]
						<< " ,rank : " << parallel.layer(pLayer_[ldown]).rank()
						<< " ,grid_rank[0] : " << parallel.layer(pLayer_[ldown]).grid_rank()[0]
						<< " ,grid_rank[1] : " << parallel.layer(pLayer_[ldown]).grid_rank()[1]
						<< " ,rank from : " << send_rank
						<< " ,grid_rec_ranks[0] : " << grid_rec_ranks[0]
						<< " ,grid_rec_ranks[1] : " << grid_rec_ranks[1]
						<<endl;*/
				for(int i=0;i<buffer_dimension[0];i++)
				{
					for(int j=0;j<buffer_dimension[1];j++)
					{
						for(int k=0;k<buffer_dimension[2];k++)
						{
							iref = i;
							jref = j+lattice_[level].coordSkip()[1] + buffer_dimension[1] - 1;
							kref = k+lattice_[level].coordSkip()[0] + buffer_dimension[2] - 1;
							if(xc.setCoordHalo(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									buffer2[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))]=src[level](xc,c);
							}
							else
							{
								cout<< "MultiGrid::prolong3d_dpl : prologation error 5"<<endl;
							}
						}
					}
				}
				parallel.layer(pLayer_[ldown]).send(buffer2,bufferSize,send_rank);
			}
		}
		else
		{
			if(plr0_[pLayer_[level]] == true) grid_rec_ranks[0] = parallel.layer(pLayer_[ldown]).grid_rank()[0]  - parallel.layer(pLayer_[ldown]).grid_rank()[0]%2;
			else grid_rec_ranks[0] = parallel.layer(pLayer_[ldown]).grid_rank()[0];

			if(plr1_[pLayer_[level]] == true) grid_rec_ranks[1] = parallel.layer(pLayer_[ldown]).grid_rank()[1]  - parallel.layer(pLayer_[ldown]).grid_rank()[1]%2;
			else grid_rec_ranks[1] = parallel.layer(pLayer_[ldown]).grid_rank()[1];

			send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[ldown]).grid_size()[0] * grid_rec_ranks[1];

			//cout<<"hihi"<<endl;
/*
			cout<< "== Receiver: "
					<< " ,level:" << pLayer_[ldown]
					<< " ,rank : " << parallel.layer(pLayer_[ldown]).rank()
					<< " ,grid_rank[0] : " << parallel.layer(pLayer_[ldown]).grid_rank()[0]
					<< " ,grid_rank[1] : " << parallel.layer(pLayer_[ldown]).grid_rank()[1]
					<< " ,rank to : " << send_rank
					<< " ,grid_rec_ranks[0] : " << grid_rec_ranks[0]
					<< " ,grid_rec_ranks[1] : " << grid_rec_ranks[1]
					<<endl;
*/
			parallel.layer(pLayer_[ldown]).receive(buffer,bufferSize,send_rank);

			int map[7];
			int index;
			//c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))
			map[0]=components; //+0
			map[1]=components * buffer_dimension[0]; //+1
			map[2]=components * buffer_dimension[0] * buffer_dimension[1]; //+2
			map[3]=components*(1+buffer_dimension[0]); //+0+1
			map[4]=components*(1+buffer_dimension[0]*buffer_dimension[1]); //+0+2
			map[5]=components*(buffer_dimension[0]*(1+buffer_dimension[1])); //+1+2
			map[6]=components*(1+buffer_dimension[0]*(1+buffer_dimension[1])); //+0+1+2
			/*
			for(int i=0;i<buffer_dimension[0];i++)
			{
				for(int j=0;j<buffer_dimension[1];j++)
				{
					for(int k=0;k<buffer_dimension[2];k++)
					{
						iref = i*2;
						jref = lattice_[ldown].coordSkip()[1] + 2*j;
						kref = lattice_[ldown].coordSkip()[0] + 2*k;
						index = components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k));

						if(iref*pow(2,ldown) != buffer[index] || jref*pow(2,ldown) != buffer[1+index] || kref*pow(2,ldown) != buffer[2+index] )
						cout<<"buffer error 0: "
								<< i<<" "<<j<<" "<<k<<" ; "
								<< iref<<" "<<jref<<" "<<kref<< " ; "
								<<buffer[index]<<" "<<buffer[1+index]<<" "<<buffer[2+index]
								<<endl;
					}
				}
			}*/
			for(int i=0;i<buffer_dimension[0]-1;i++)
			{
				for(int j=0;j<buffer_dimension[1]-1;j++)
				{
					for(int k=0;k<buffer_dimension[2]-1;k++)
					{
						iref = i*2;
						jref = lattice_[ldown].coordSkip()[1] + 2*j;
						kref = lattice_[ldown].coordSkip()[0] + 2*k;
						if(xf.setCoord(iref,jref,kref))
						{
							index = components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k));
							/*
							cout<<"buffer: "
									<< i<<" "<<j<<" "<<k<<" ; "
									<< iref<<" "<<jref<<" "<<kref<< " ; "
									<<buffer[index]<<" "<<buffer[1+index]<<" "<<buffer[2+index]<<" , "
									<<buffer[index+map[0]]<<" "<<buffer[1+index+map[1]]<<" "<<buffer[2+index+map[2]]
									<<endl;
							*/

							for(int c = 0; c< components;c++)
								dst[ldown](xf,c) = buffer[c+index];
								//field[ldown](xf,c) = field[level](xc,c);

							for(int i=0;i<3;i++)
								for(int c = 0; c< components;c++)
								dst[ldown](xf+i,c) = (buffer[c+index] + buffer[c+index+map[i]])/2.0;
							//		field[ldown](xf+i,c) = (field[level](xc,c) + field[level](xc+i,c))/2.0;


							for(int c = 0; c< components;c++)
								dst[ldown](xf+0+1,c) = (buffer[c+index] + buffer[c+index+map[0]]
																				 	+ buffer[c+index+map[1]] + buffer[c+index+map[3]])/4.0;
							//		field[ldown](xf+0+1,c) = (field[level](xc,c) + field[level](xc+0,c)
							//															+ field[level](xc+1,c) + field[level](xc+0+1,c))/4.0;

							for(int c = 0; c< components;c++)
								dst[ldown](xf+0+2,c) = (buffer[c+index] + buffer[c+index+map[0]]
																				 	+ buffer[c+index+map[2]] + buffer[c+index+map[4]])/4.0;
							//		field[ldown](xf+0+2,c) = (field[level](xc,c) + field[level](xc+0,c)
							//															+ field[level](xc+2,c) + field[level](xc+0+2,c))/4.0;

							for(int c = 0; c< components;c++)
								dst[ldown](xf+1+2,c) = (buffer[c+index] + buffer[c+index+map[1]]
																				 	+ buffer[c+index+map[2]] + buffer[c+index+map[5]])/4.0;
							//			field[ldown](xf+1+2,c) = (field[level](xc,c) + field[level](xc+1,c)
							//																+ field[level](xc+2,c) + field[level](xc+1+2,c))/4.0;

							for(int c = 0; c< components;c++)
								dst[ldown](xf+0+1+2,c) = (buffer[c+index] + buffer[c+index+map[0]]
																						+ buffer[c+index+map[1]] + buffer[c+index+map[2]]
																						+ buffer[c+index+map[3]] + buffer[c+index+map[4]]
																				 		+ buffer[c+index+map[5]] + buffer[c+index+map[6]])/8.0;
							//	field[ldown](xf+0+1+2,c) = (field[level](xc,c) + field[level](xc+1,c)
							//														+ field[level](xc+2,c) + field[level](xc+1+2,c)
							//														+ field[level](xc+0,c) + field[level](xc+0+1,c)
							//														+ field[level](xc+0+2,c) + field[level](xc+0+1+2,c))/8.0;
						}
						else
						{
							cout<< "MultiGrid::prolong3d_dpl : prologation error 6 : "
							<< i<<" "<<j<<" "<<k
							<< iref<<" "<<jref<<" "<<kref<<endl;
						}
					}
				}
			}

		}


		free(buffer);
		free(buffer1);
		free(buffer2);


	}
	MPI_Barrier(MPI_COMM_WORLD);
}



#endif
