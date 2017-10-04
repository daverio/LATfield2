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
	#endif


	parallel.MultiGrid_createLevels(npl_);
	parallel.MultiGrid_initTopLevel();
	for(int i=1;i<npl_;i++)if(  parallel.layers()[i-1].isPartLayer() )parallel.MultiGrid_initLevel(i,lr0[i],lr1[i]);

	lattice_ = new MultiLAT[nl_];
	for(int i=0;i<nl_;i++)
	{
		if(parallel.layer(lLayer_[i]).isPartLayer())lattice_[i].initialize(3,16,1,lLayer_[i]);
	}

}







//lattice for MultiGrid



void MultiLAT::initialize(int dim, const int size, int halo, int parallel_layer)
{
	status_=0;
	arch_saved_=false;
	int* sizeArray=new int[dim];
	for(int i=0; i<dim; i++) { sizeArray[i]=size; }
	this->initialize(dim, sizeArray, halo,parallel_layer);
	delete[] sizeArray;
}

void MultiLAT::initialize(int dim, const int* size, int halo, int parallel_layer)
{

		int i;

		if((status_ & initialized) > 0)
	    {
			delete[] size_;
			delete[] sizeLocal_;
			delete[] jump_;
	        delete[] sizeLocalAllProcDim0_;
	        delete[] sizeLocalAllProcDim1_;
	    }
		//Store input lattice properties
		dim_ =dim;
		size_=new int[dim_];
		for(i=0;i<dim_;i++) size_[i]=size[i];
		halo_=halo;

		parallel_layer_ = parallel_layer;

		//Calculate local size
		sizeLocal_=new int[dim_];
		sizeLocal_[dim_-1]=int(ceil( (parallel.layer(parallel_layer).grid_size()[0]-parallel.layer(parallel_layer).grid_rank()[0])*size_[dim_-1]/float(parallel.layer(parallel_layer).grid_size()[0]) ));
		sizeLocal_[dim_-1]-=int(ceil((parallel.layer(parallel_layer).grid_size()[0]-parallel.layer(parallel_layer).grid_rank()[0]-1)*size_[dim_-1]/float(parallel.layer(parallel_layer).grid_size()[0]) ));
		sizeLocal_[dim_-2]=int(ceil( (parallel.layer(parallel_layer).grid_size()[1]-parallel.layer(parallel_layer).grid_rank()[1])*size_[dim_-2]/float(parallel.layer(parallel_layer).grid_size()[1]) ));
		sizeLocal_[dim_-2]-=int(ceil((parallel.layer(parallel_layer).grid_size()[1]-parallel.layer(parallel_layer).grid_rank()[1]-1)*size_[dim_-2]/float(parallel.layer(parallel_layer).grid_size()[1]) ));
		for(i=0;i<dim_-2;i++) sizeLocal_[i]=size_[i];

		//cout<< "layer: "<< parallel_layer <<  "  psize0: "<<parallel.layer(parallel_layer).grid_size()[0]<<endl;
		//cout<< "layer: "<< parallel_layer <<  "  psize1: "<<parallel.layer(parallel_layer).grid_size()[1]<<endl;
		
	  sizeLocalAllProcDim0_ = new int[parallel.layer(parallel_layer).grid_size()[0]];
		sizeLocalAllProcDim1_ = new int[parallel.layer(parallel_layer).grid_size()[1]];

		for(i=0;i<parallel.layer(parallel_layer).grid_size()[0];i++)
	    {
		    sizeLocalAllProcDim0_[i] = int(ceil( (parallel.layer(parallel_layer).grid_size()[0]-i)*size_[dim_-1]/float(parallel.layer(parallel_layer).grid_size()[0]) ));
		    sizeLocalAllProcDim0_[i] -= int(ceil((parallel.layer(parallel_layer).grid_size()[0]-i-1)*size_[dim_-1]/float(parallel.layer(parallel_layer).grid_size()[0]) ));
	    }
		for(i=0;i<parallel.layer(parallel_layer).grid_size()[1];i++)
	    {
		    sizeLocalAllProcDim1_[i] = int(ceil( (parallel.layer(parallel_layer).grid_size()[1]-i)*size_[dim_-2]/float(parallel.layer(parallel_layer).grid_size()[1]) ));
		    sizeLocalAllProcDim1_[i] -= int(ceil((parallel.layer(parallel_layer).grid_size()[1]-i-1)*size_[dim_-2]/float(parallel.layer(parallel_layer).grid_size()[1]) ));
	    }

		//Calculate index jumps
		jump_=new long[dim_];
		jump_[0]=1;
		for(i=1;i<dim_;i++) jump_[i]=jump_[i-1]*(sizeLocal_[i-1]+2*halo_);

		//Calculate number of sites in lattice
		sitesLocal_=1;
		sitesLocalGross_=1;
		for(i=0;i<dim_;i++)
	    {
			sitesLocal_*=sizeLocal_[i];
			sitesLocalGross_*=sizeLocal_[i]+(2*halo_);
	    }
		sites_=sitesLocal_;
		sitesGross_=sitesLocalGross_;
		parallel.layer(parallel_layer).sum(sites_);
		parallel.layer(parallel_layer).sum(sitesGross_);

		//Calculate index of first and last local sites on lattice
		siteFirst_=0;
		siteLast_=sitesLocalGross_-1;
		for(i=0;i<dim_;i++)
	    {
			siteFirst_+=jump_[i]*halo_;
			siteLast_-=jump_[i]*halo_;
	    }

		////calculate coordSkip



		//Get each processor to tell the others in his dim0_group its local sizeLocal_[dim-1])
		int* sizes_dim0=new int[parallel.layer(parallel_layer).grid_size()[0]];
		for(i=0;i<parallel.layer(parallel_layer).grid_size()[0];i++)
		{
			if(i==parallel.layer(parallel_layer).grid_rank()[0]) { sizes_dim0[i]=sizeLocal_[dim_-1]; }
			parallel.layer(parallel_layer).broadcast_dim0(sizes_dim0[i],i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[0]=0;
		for(i=0; i<parallel.layer(parallel_layer).grid_rank()[0]; i++)coordSkip_[0]+=sizes_dim0[i];



		//Get each processor to tell the others in his dim1_group its local sizeLocal_[dim-2])
		int* sizes_dim1=new int[parallel.layer(parallel_layer).grid_size()[1]];
		for(i=0;i<parallel.layer(parallel_layer).grid_size()[1];i++)
		{
			if(i==parallel.layer(parallel_layer).grid_rank()[1]) { sizes_dim1[i]=sizeLocal_[dim_-2]; }
			parallel.layer(parallel_layer).broadcast_dim1(sizes_dim1[i],i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[1]=0;
		for(i=0; i<parallel.layer(parallel_layer).grid_rank()[1]; i++)coordSkip_[1]+=sizes_dim1[i];




		////calculate sitesSkip : sitesskip used for fastread , fastload (function witch need to be coded :-) )

		sitesSkip_=coordSkip_[0];
		for(i=0;i<dim_-1;i++)sitesSkip_*=size_[i];

		long siteSkiptemp= coordSkip_[1]*sizeLocal_[dim_-1];
		for(i=0;i<dim_-2;i++)siteSkiptemp*=size_[i];

		sitesSkip_+=siteSkiptemp;

		//calculate sitesSkip2d :

		int* sizes1 = new int[parallel.layer(parallel_layer).grid_size()[1]];
		int* sizes0 = new int[parallel.layer(parallel_layer).grid_size()[0]];
		long* offset1 = new long[parallel.layer(parallel_layer).grid_size()[1]];
		long* offset0 = new long[parallel.layer(parallel_layer).grid_size()[0]];

		int b=1;
		int n;
		for(i=0;i<dim_-2;i++)b*=sizeLocal_[i];

		//calulate offset in dim-2
		for(i=0;i<parallel.layer(parallel_layer).grid_size()[1];i++)sizes1[i]=sizes_dim1[i] * b;
		for(n=0;n<parallel.layer(parallel_layer).grid_size()[1];n++)offset1[n]=0;

		for(n=1;n<parallel.layer(parallel_layer).grid_size()[1];n++)
		{
			for(i=0;i<n;i++)offset1[n]+=sizes1[i];
		}

		//calulate offset in dim-1
		for(i=0;i<parallel.layer(parallel_layer).grid_size()[0];i++)sizes0[i]=size_[dim_-2] * sizes_dim0[i] * b;
		for(n=0;n<parallel.layer(parallel_layer).grid_size()[0];n++)offset0[n]=0;

		for(n=1;n<parallel.layer(parallel_layer).grid_size()[0];n++)
		{
			for(i=0;i<n;i++)offset0[n]+=sizes0[i];
		}

		sitesSkip2d_ = offset0[parallel.layer(parallel_layer).grid_rank()[0]] + offset1[parallel.layer(parallel_layer).grid_rank()[1]];

		//Set status
		status_ = status_ | initialized;

		//Free memory
		delete[] sizes_dim0;
		delete[] sizes_dim1;

	  // WV: Free more memory!
	  delete[] sizes1;
	  delete[] sizes0;
	  delete[] offset1;
	  delete[] offset0;

}




#endif
