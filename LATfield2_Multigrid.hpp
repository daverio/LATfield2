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
	plr0_ = new int[nl_];
	plr1_ = new int[nl_];
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

	#ifdef DEBUG_MULTIGRID
	/*
		cout<< "number of lattice layers: "<< nl_<<endl;
		cout<< "number of parallel layers: "<< npl_<<endl;
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
		*/
	#endif


	parallel.MultiGrid_createLevels(npl_);
	parallel.MultiGrid_initTopLevel();
	for(int i=1;i<npl_;i++)if(  parallel.layers()[i-1].isPartLayer() )parallel.MultiGrid_initLevel(i,plr0_[i],plr1_[i]);

	lattice_ = new MultiLAT[nl_];
	for(int i=0;i<nl_;i++)
	{
		if(parallel.layer(pLayer_[i]).isPartLayer())lattice_[i].initialize(lat_top->dim(),lat_size_[i],lat_top->halo(),pLayer_[i]);
	}

}

template<class FieldType>
void MultiGrid::intitialize_Field(Field<FieldType> * fieldBase, MultiField<FieldType> *& field)
{

	field = new MultiField<FieldType>[nl_];

	for(int i=0;i<nl_;i++)
	{

			if(parallel.layer(pLayer_[i]).isPartLayer())field[i].initialize(lattice_[i],
																																			fieldBase->nMatrix(),
																																			fieldBase->rows(),
																																			fieldBase->cols(),
																																			fieldBase->symmetry());

			COUT << "Lattice size: layer["<< i <<"] ("<< lattice_[i].size(0)<<","<< lattice_[i].size(1)<<","<< lattice_[i].size(2)<<");"<<endl;



	}
	for(int i=1;i<nl_;i++)
	{
		if(parallel.layer(pLayer_[i]).isPartLayer())field[i].alloc();
	}

	field[0].data_ = fieldBase->data_;
}


template<class FieldType>
void MultiGrid::restrict(MultiField<FieldType> *& field, int level)
{
	if(pLayer_[level]==pLayer_[level+1])restrict3d_spl(field,level);
	else restrict3d_dpl(field,level);
}

template<class FieldType>
void MultiGrid::restrict3d_spl(MultiField<FieldType> *& field, int level)
{
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		int lup = level+1;

		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);

		int components = field[level].components();

		for(xc.first();xc.test();xc.next())
		{
			if(xf.setcoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2))
			{
				for(int c = 0; c< components;c++)field[lup](xc,c) = field[level](xf,c)
			}
			else
			{
				cout<<"restriction issues...."<<endl;
			}

		}
	}
}

template<class FieldType>
void MultiGrid::restrict3d_dpl(MultiField<FieldType> *& field, int level)
{
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		int lup = level+1;

		int iref,jref,kref;

		int buffer_dimension[3]
		buffer_dimension[0] = lattice_[level].sizeLocal(0)/2;
		buffer_dimension[1] = lattice_[level].sizeLocal(1)/2;
		buffer_dimension[2] = lattice_[level].sizeLocal(2)/2;

		if(plr1_[pLayer_[lup]]==true)buffer_dimension[1]/2;
		if(plr0_[pLayer_[lup]]==true)buffer_dimension[2]/2;

		int components = field[level].components();

		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);

		int grid_rec_ranks[2];
		int send_rank;

		FieldType * buffer;
		long bufferSize = buffer_dimension[0]*buffer_dimension[1]*buffer_dimension[3];

		buffer = (FieldType *)malloc(sizeof(FieldType));


		if(parallel.layer(pLayer_[lup]).isPartLayer())
		{
			for(int i = 0;i<buffer_dimension[0],i++)
			{
				for(int j = lattice_[lup].coordSkip()[1];j<lattice_[lup].coordSkip()[1] + buffer_dimension[1];j++)
				{
					for(int k = lattice_[lup].coordSkip()[0]; k<lattice_[lup].coordSkip()[0] + buffer_dimension[2],k++)
					{
						if(xc.setCoord(i,j,k))
						{
							if(xf.setCoord(2*i,2*j,2*k))
							{
								for(int c = 0; c< components;c++)field[lup](xc,c) = field[level](xf,c)
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
				grid_rec_ranks[0] = parallel.level(pLayer_[lup]).grid_rank()[0]*2 + 1;
				if(plr1_[pLayer_[lup]] == true) grid_rec_ranks[1] = parallel.level(pLayer_[lup]).grid_rank()[1]*2 + 1;
				else grid_rec_ranks[1] = parallel.level(pLayer_[lup]).grid_rank()[1];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0],i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2],k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									field[level](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
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
				grid_rec_ranks[1] = parallel.level(pLayer_[lup]).grid_rank()[1]*2 + 1;
				if(plr0_[pLayer_[lup]] == true) grid_rec_ranks[0] = parallel.level(pLayer_[lup]).grid_rank()[0]*2 + 1;
				else grid_rec_ranks[0] = parallel.level(pLayer_[lup]).grid_rank()[0];
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0],i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2],k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									field[level](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
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
				grid_rec_ranks[0] = parallel.level(pLayer_[lup]).grid_rank()[0]*2 + 1;
				grid_rec_ranks[1] = parallel.level(pLayer_[lup]).grid_rank()[1]*2 + 1;
				send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];
				parallel.layer(pLayer_[level]).receive(buffer,bufferSize,send_rank);
				for(int i = 0;i<buffer_dimension[0],i++)
				{
					for(int j = 0;j<buffer_dimension[1];j++)
					{
						for(int k = 0;k<buffer_dimension[2],k++)
						{
							iref = i;
							jref = j + lattice_[lup].coordSkip()[1] + buffer_dimension[1];
							kref = k + lattice_[lup].coordSkip()[0] + buffer_dimension[2];

							if(xc.setCoord(iref,jref,kref))
							{
								for(int c = 0; c< components;c++)
									field[level](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
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
			for(int i = 0;i<buffer_dimension[0],i++)
			{
				for(int j = 0;j<buffer_dimension[1];j++)
				{
					for(int k = 0;k<buffer_dimension[2],k++)
					{
						iref = i*2;
						jref = lattice_[level].coordSkip()[1]+j*2;
						kref = lattice_[level].coordSkip()[0]+k*2;

						if(xc.setCoord(iref,jref,kref))
						{
							for(int c = 0; c< components;c++)
								field[level](xc,c) = buffer[c+components*(i+buffer_dimension[0]*(j+buffer_dimension[1]*k))];
						}
						else
						{
							cout<< "MultiGrid::restrict3d_dpl : restriction error 4"<<endl;
						}
					}
				}
			}

			if(parallel.layer(pLayer_[level]).grid_rank() )

			if(plr0_[pLayer_[lup]] == true) grid_rec_ranks[0] = ( parallel.layer(pLayer_[level]).grid_rank()[0]  - (parallel.layer(pLayer_[level]).grid_rank()[0]%2) ) / 2.0);
			else grid_rec_ranks[0] = parallel.layer(pLayer_[level]).grid_rank()[0];

			if(plr1_[pLayer_[lup]] == true) grid_rec_ranks[1] = ( parallel.layer(pLayer_[level]).grid_rank()[1]  - (parallel.layer(pLayer_[level]).grid_rank()[1]%2) ) / 2.0);
			else grid_rec_ranks[1] = parallel.layer(pLayer_[level]).grid_rank()[1];

			send_rank = grid_rec_ranks[0] +  parallel.layer(pLayer_[level]).grid_size()[0] * grid_rec_ranks[1];

			parallel.layer(pLayer_[level]).send(buffer,bufferSize,send_rank);

		}
		free(buffer);
	}
}

template<class FieldType>
void MultiGrid::prolonge3d(MultiField<FieldType> *& field, int level)
{
	if(pLayer_[level]==pLayer_[level-1])restrict3d_spl(field,level);
	else restrict3d_dpl(field,level);
}

template<class FieldType>
void MultiGrid::prolonge3d_spl(MultiField<FieldType> *& field, int level)
{
	if(parallel.layer(pLayer_[level]).isPartLayer())
	{
		int ldown = level-1;

		Site xc(lattice_[level]);
		Site xf(lattice_[ldown]);

		field[level].updateHalo();

		for(xc.first();xc.test();xc.next())
		{
			if(xf.setCoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2))
			{
				for(int i=0;i<3;i++)
					for(int c = 0; c< components;c++)
						field[ldown](xf+i,c) = (field[level](xc,c) + field[level](xc+i,c))/2.0;

				for(int i=0;i<2;i++)
					for(int j=i+1;j<3;j++)
						for(int c = 0; c< components;c++)
							field[ldown](xf+i+j,c) = (field[level](xc,c) + field[level](xc+i,c)
																			+ field[level](xc+j,c) + field[level](xc+i+j,c))/4.0;

				for(int c = 0; c< components;c++)
					field[ldown](xf+0+1+2,c) = (field[level](xc-0,c) + field[level](xc-0+1,c)
																		+ field[level](xc-0+2,c) + field[level](xc-0+1+2,c)
																		+ field[level](xc+0,c) + field[level](xc+0+1,c)
																		+ field[level](xc+0+2,c) + field[level](xc+0+1+2,c))/8.0;
			}
			else
			{
				cout<< "MultiGrid::prolonge3d_spl : restriction error 1"<<endl;
			}
		}



	}
}

template<class FieldType>
void MultiGrid::prolonge3d_dpl(MultiField<FieldType> *& field, int level)
{

}



#endif
