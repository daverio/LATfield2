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

	pl_number_=1;
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
					lr0[pl_number_]=true;
				}
		}
		if(lat_size_[l][dim_-2] / pl_size[1]  < minGridperProc)
		{
			if( pl_size[1]/2 >=2)
			{
				pl_size[1] /= 2;
				lr1[pl_number_]=true;
			}
		}

		if(lr0[pl_number_] || lr1[pl_number_])
		{
			pl_number_++;
		}
	}

	#ifdef DEBUG_MULTIGRID
		cout<< "number of lattice layers: "<< nl_<<endl;
		cout<< "number of parallel layers: "<< pl_number_<<endl;
		MPI_Barrier(parallel.lat_world_comm());
	#endif


	parallel.MultiGrid_createLevels(pl_number_);
	parallel.MultiGrid_initTopLevel();
	for(int i=1;i<pl_number_;i++)if(  parallel.layers()[i-1].isPartLayer() )parallel.MultiGrid_initLevel(i,lr0[i],lr1[i]);



}


#endif
