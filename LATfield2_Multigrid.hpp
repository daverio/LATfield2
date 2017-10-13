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

	lLayer_ = new int[nl_];
	lLayer_[0] = 0;

	npl_=1;
	int pl_size[2];
	pl_size[0]=parallel.grid_size()[0];
	pl_size[1]=parallel.grid_size()[1];

	bool lr0[nl_],lr1[nl_];

	for(int i=0;i<nl_;i++)
	{
		lr0[i]=false;
		lr1[i]=false;
	}

	for(int l = 1;l<nl_;l++)
	{
		if(lat_size_[l][dim_-1] / pl_size[0]  < minGridperProc)
		{
				if( pl_size[0]/2 >=2)
				{
					pl_size[0] /= 2;
					lr0[npl_]=true;
				}
		}
		if(lat_size_[l][dim_-2] / pl_size[1]  < minGridperProc)
		{
			if( pl_size[1]/2 >=2)
			{
				pl_size[1] /= 2;
				lr1[npl_]=true;
			}
		}
		if(lr0[npl_] || lr1[npl_])
		{
			npl_++;
		}
		lLayer_[l] = npl_-1;
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
			COUT<< "parallel layer: "<< lLayer_[i] << endl;;
			COUT<<"-----------"<<endl;

		}

		MPI_Barrier(parallel.lat_world_comm());
		*/
	#endif


	parallel.MultiGrid_createLevels(npl_);
	parallel.MultiGrid_initTopLevel();
	for(int i=1;i<npl_;i++)if(  parallel.layers()[i-1].isPartLayer() )parallel.MultiGrid_initLevel(i,lr0[i],lr1[i]);

	lattice_ = new MultiLAT[nl_];
	for(int i=0;i<nl_;i++)
	{
		if(parallel.layer(lLayer_[i]).isPartLayer())lattice_[i].initialize(lat_top->dim(),lat_size_[i],lat_top->halo(),lLayer_[i]);
	}

}

template<class FieldType>
void MultiGrid::intitialize_Field(Field<FieldType> * fieldBase, MultiField<FieldType> *& field)
{

	field = new MultiField<FieldType>[nl_];

	for(int i=0;i<nl_;i++)
	{

			if(parallel.layer(lLayer_[i]).isPartLayer())field[i].initialize(lattice_[i],
																																			fieldBase->nMatrix(),
																																			fieldBase->rows(),
																																			fieldBase->cols(),
																																			fieldBase->symmetry());

			COUT << "Lattice size: layer["<< i <<"] ("<< lattice_[i].size(0)<<","<< lattice_[i].size(1)<<","<< lattice_[i].size(2)<<");"<<endl;



	}
	for(int i=1;i<nl_;i++)
	{
		if(parallel.layer(lLayer_[i]).isPartLayer())field[i].alloc();
	}

	field[0].data_ = fieldBase->data_;
}


template<class FieldType>
void MultiGrid::restrict(MultiField<FieldType> *& field, int level)
{
	if(lLayer_[level]==lLayer_[level+1])restrict3d_spl(field,level);
	else restrict3d_dpl(field,level);
}

template<class FieldType>
void MultiGrid::restrict3d_spl(MultiField<FieldType> *& field, int level)
{
	if(parallel.layer(lLayer_[level]).isPartLayer())
	{
		int lup = level+1;

		Site xc(lattice_[lup]);
		Site xf(lattice_[level]);

		Site xadd(lattice_[level]);

		for(xc.first();xc.test();xc.next())
		{
			if(xf.setcoord(xc.coord(0)*2,xc.coord(1)*2,xc.coord(2)*2))
			{
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xf,c) / 8.0;
				xadd = xf+0;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;
				xadd = xf-0;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;
				xadd = xf+1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;
				xadd = xf-1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;
				xadd = xf+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;
				xadd = xf-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 16.0;

				xadd = xf+0+1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf+0-1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf+0+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf+0-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-0+1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-0-1;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-0+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-0-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf+1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf+1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;
				xadd = xf-1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 32.0;

				xadd = xf+0+1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf+0+1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf+0-1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf+0-1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf-0+1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf-0+1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf-0-1+2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
				xadd = xf-0-1-2;
				for(int c = 0; c< field[level].components();c++)field[lup](xc,c) = field[level](xadd,c) / 64.0;
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
	if(parallel.layer(lLayer_[level]).isPartLayer())
	{
		if(parallel.layer(lLayer_[level]).isPartLayer())
		{
			//1 do the local part

			//recieves 1)

			//unpack 1)

			//recieves 2)

			//unpack 2)

			//recieves 3)

			//unpack 3)

		}
		else
		{
			//compute buffer
			

			//send it


		}









	}
}



#endif
