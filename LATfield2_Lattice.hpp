#ifndef LATFIELD2_LATTICE_HPP
#define LATFIELD2_LATTICE_HPP


/*! \file LATfield2_Lattice.hpp
 \brief Lattice class definition
 
 LATfield2_Lattice.hpp contain the class Lattice definition.
 
 */ 



/*! \class Lattice  
 
    \brief The Lattice class describe a cartesian mesh (with 2 or more dimensions). The updateHalo method of the Field class generate the periodicity.
 
 
    It store the global and local geometry of the mesh. The last 2 dimension of the lattice are scattered into the MPI processes grid. 
 
 
 */
class Lattice
	{
	public:
        //! Constructor.
		Lattice();
        
        /*!
         Constructor with initialization
         \sa initialize(int dim, const int* size, int halo);
         \param dim : number of dimension
         \param size : array containing the size of each dimension.
         \param halo : size of the halo (ghost cells, same for each dimension)
         */ 
		Lattice(int dim, const int* size, int halo); 
        
        /*!
         Constructor with initialization
         \sa initialize(int dim, const int size, int halo);
         \param dim : number of dimension
         \param size : size of each dimension (same for each dimension)
         \param halo : size of the halo (same for each dimension)
         */ 
		Lattice(int dim, const int size, int halo);  
		
        //! Destructor.
		~Lattice();
		
        /*!
         Initialization of a dim-dimensional lattice, the size of each dimension is set by the second parameter: int *size. The ghost cell number (halo) is the same for each dimension.  
         
         \param dim : number of dimension
         \param size : array containing the size of each dimension.
         \param halo : size of the halo (same for each dimension)
         */ 
        void initialize(int dim, const int* size, int halo);
        
        /*!
         Initialization of a dim-dimensional lattice, each dimension have the same size. The ghost cell number (halo) is the same for each dimension. 
         
         
         \param dim : number of dimension
         \param size : size of each dimension (same for each dimension)
         \param halo : size of the halo (same for each dimension)
         */ 
		void initialize(int dim, const int size, int halo);
        
#ifdef FFT3D
        
        
        /*!
         Initialization of a lattice for Fourier space in case of real to complex transform. The Fourier space lattice size is defined according to the real space one. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
         \param lat_real : pointer to a real space lattice.
         \param halo : size of the halo (same for each dimension)
         */ 
		void initializeRealFFT(Lattice & lat_real, int halo);
        
        /*!
         Initialization of a lattice for Fourier space in case of complex to complex transform. The Fourier space lattice size is defined according to the real space one.. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
         \param lat_real : pointer to a real space lattice.
         \param halo : size of the halo (same for each dimension)
         */  
		void initializeComplexFFT(Lattice & lat_real, int halo);
#endif
        
    
        
        /*!
         \return int. Number of dimensions of the lattice.
         */
        int  dim();
        /*!
         \return int. Size of the halo (ghost cells).
         */
		int  halo();
        
        /*!
         \return int*. Pointer to the array of the size of each dimension of the lattice.
         */
        int* size();               
        
        /*!
         Function which returns the size of a given dimension of the lattice.
         \param direction : asked dimension.
         \return int. Global size of the lattice in the given dimension.
         */
        int  size(int direction);  //Size in a particular dimension
        
        /*!
         \return int*. Pointer to the array of the size of each dimension of the sublattice stored in this MPI process.
         */
        int* sizeLocal();               //Local version
        
        /*!
         Function which returns the size of a given dimension of the sublattice stored in this MPI process.
         \param direction : asked dimension.
         \return int. Global size of the sublattice (of this MPI process) in the given dimension.
         */
		int  sizeLocal(int direction);  //Local version
        
        
        /*!
         \return long. Number of sites on the lattice (excluding halo sites).
         */
        long  sites();              //Number of (global) sites
        /*!
         \return long. Number of sites on the lattice (including halo sites).
         */
		long  sitesGross();         //Number of (global) sites including halo
        
        /*!
         \return long. Number of sites (excluding halo sites) of the sublattice stored in this MPI process.
         */
        long  sitesLocal();              //Local version
        
        /*!
         \return long. Number of sites (including halo sites) of the sublattice stored in this MPI process.
         */
		long  sitesLocalGross();         //Local version
        
		/*!
         \return long. Array index of the first site which is not within the halo.
         */
        long siteFirst();
        
        /*!
         \return long. Array index of the last site which is not within the halo.
         */
		long siteLast();
		
		
		
		/*!
         Function which return the number of data_ array elements to jump to move to the next site in the given direction. (does not take into account the number of component of the fields, therefor should be multiplied by Field.components().) Should not be used by user.
         \param direction : asked direction.
         \return long. Number of array elements to jump.
         */
		long  jump(int direction);       //Number of sites jumped to move in direction
		
        /*!
         \return long. Number of sites before first local site in lattice. Should not be used by users.
         */
        long  sitesSkip();               //Number of sites before first local site in lattice (say in a file)
		
        /*!
         \return long. Number of sites before first local site in lattice. Should not be used by users.
         */
        long  sitesSkip2d();
        
        /*!
         \return *long. Pointer to an array which store the last 2 dimensions coordinate of the first local(in this MPI process) sites. Index 0 is for dim-1, index 1 is for dim-2/
         */
		long*  coordSkip();              //Number to add to coord[dim_-1] to get global value
		
		/*!
         Function which save in serial and in ASCII the global and local description of the Lattice. Usefull to read a file writen by fast_save or fast_write methods of the Field class.
         \param filename : filename of the architectur file.
         */
        void save_arch(const string filename);
        
        /*!
         \return return true if the description of the lattice has been written on disk.
         
         \sa save_arch(const string filename)
         */
		bool is_arch_saved();
        
        int getRank(int* coord) ; //return the world rank of the process who get the lattices site "coord"
        int getRankDim0(int coord) ;
        int getRankDim1(int coord) ;
        
        
	private:
		int        status_;
		static int initialized;
		
		//Global variables==============
		int  dim_;              //ok//Number of dimensions
		int* size_;             //ok//Number of lattice sites in each direction
		long  sites_;            //ok//Number of sites in lattice
		long  sitesGross_;       //ok//Number of sites in lattice plus halos
		int  halo_;             //ok//Number of sites extra in each direction
		
		//Local variables===============
		int* sizeLocal_;       //ok//Number of local lattice sites in each direction
		int* sizeLocalAllProcDim0_;
        int* sizeLocalAllProcDim1_;
        long  sitesLocal_;      //ok//Number of local sites in lattice
		long  sitesLocalGross_; //ok//Number of local sites in lattice plus halo
		long* jump_;            //ok//Jumps needed to move up one in each direction
		
    
        long  siteFirst_;       //ok//Index of first local site in lattice
		
        long  siteLast_;        //ok//Index of last local site in lattice
		long  sitesSkip_;      //Number of global lattice sites before first local site (say in a file)
		long  sitesSkip2d_;
		long  coordSkip_[2];       //Number to add to coord[dim_-1] and coord[dim_-2] to get global value
		
		
		//save variable for fast save
		int arch_saved_;
		
	};

//CONSTANTS=====================

int Lattice::initialized = 1;        //Status flag for initialized

//CONSTRUCTORS===================

Lattice::Lattice()
{
	status_=0;
	arch_saved_=false;
}

Lattice::Lattice(int dim, const int* size, int halo)
{
	status_=0;
	arch_saved_=false;
	this->initialize(dim, size, halo);
}  


Lattice::Lattice(int dim, const int size, int halo)
{
	status_=0;
	arch_saved_=false;
	int* sizeArray=new int[dim];
	for(int i=0; i<dim; i++) { sizeArray[i]=size; }
	this->initialize(dim, sizeArray, halo);
	delete[] sizeArray;
}



//DESTRUCTOR=========================

Lattice::~Lattice() 
{ 
	if((status_ & initialized) > 0)
    {
		delete[] size_;
		delete[] sizeLocal_;
		delete[] jump_;
        delete[] sizeLocalAllProcDim0_;
        delete[] sizeLocalAllProcDim1_;
    }
}
//INITIALIZE=========================

void Lattice::initialize(int dim, const int size, int halo)
{
	int* sizeArray=new int[dim];
	for(int i=0; i<dim; i++) { sizeArray[i]=size; }
	this->initialize(dim, sizeArray, halo);
}

void Lattice::initialize(int dim, const int* size, int halo)
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
	
	
	//Calculate local size
	sizeLocal_=new int[dim_];
	sizeLocal_[dim_-1]=int(ceil( (parallel.grid_size()[0]-parallel.grid_rank()[0])*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	sizeLocal_[dim_-1]-=int(ceil((parallel.grid_size()[0]-parallel.grid_rank()[0]-1)*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	sizeLocal_[dim_-2]=int(ceil( (parallel.grid_size()[1]-parallel.grid_rank()[1])*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	sizeLocal_[dim_-2]-=int(ceil((parallel.grid_size()[1]-parallel.grid_rank()[1]-1)*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	for(i=0;i<dim_-2;i++) sizeLocal_[i]=size_[i];
    
    sizeLocalAllProcDim0_ = new int[parallel.grid_size()[0]];
	sizeLocalAllProcDim1_ = new int[parallel.grid_size()[1]];
    
	for(i=0;i<parallel.grid_size()[0];i++)
    {
	    sizeLocalAllProcDim0_[i] = int(ceil( (parallel.grid_size()[0]-i)*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	    sizeLocalAllProcDim0_[i] -= int(ceil((parallel.grid_size()[0]-i-1)*size_[dim_-1]/float(parallel.grid_size()[0]) ));	
    }
	for(i=0;i<parallel.grid_size()[1];i++)
    {
	    sizeLocalAllProcDim1_[i] = int(ceil( (parallel.grid_size()[1]-i)*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	    sizeLocalAllProcDim1_[i] -= int(ceil((parallel.grid_size()[1]-i-1)*size_[dim_-2]/float(parallel.grid_size()[1]) ));	
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
	parallel.sum(sites_);
	parallel.sum(sitesGross_);
	
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
	int* sizes_dim0=new int[parallel.grid_size()[0]];
	for(i=0;i<parallel.grid_size()[0];i++)
	{
		if(i==parallel.grid_rank()[0]) { sizes_dim0[i]=sizeLocal_[dim_-1]; }
		parallel.broadcast_dim0(sizes_dim0[i],i);
	}
	//Sum up sizes for the processors of less than or equal rank
	coordSkip_[0]=0;
	for(i=0; i<parallel.grid_rank()[0]; i++)coordSkip_[0]+=sizes_dim0[i];
	
	
	
	//Get each processor to tell the others in his dim1_group its local sizeLocal_[dim-2])
	int* sizes_dim1=new int[parallel.grid_size()[1]];
	for(i=0;i<parallel.grid_size()[1];i++)
	{
		if(i==parallel.grid_rank()[1]) { sizes_dim1[i]=sizeLocal_[dim_-2]; }
		parallel.broadcast_dim1(sizes_dim1[i],i);
	}
	//Sum up sizes for the processors of less than or equal rank
	coordSkip_[1]=0;
	for(i=0; i<parallel.grid_rank()[1]; i++)coordSkip_[1]+=sizes_dim1[i];
	
	
	
	
	////calculate sitesSkip : sitesskip used for fastread , fastload (function witch need to be coded :-) )
	
	sitesSkip_=coordSkip_[0];
	for(i=0;i<dim_-1;i++)sitesSkip_*=size_[i];
	
	long siteSkiptemp= coordSkip_[1]*sizeLocal_[dim_-1];
	for(i=0;i<dim_-2;i++)siteSkiptemp*=size_[i];
	
	sitesSkip_+=siteSkiptemp;
	
	//calculate sitesSkip2d : 
	
	int* sizes1 = new int[parallel.grid_size()[1]];
	int* sizes0 = new int[parallel.grid_size()[0]];
	long* offset1 = new long[parallel.grid_size()[1]];
	long* offset0 = new long[parallel.grid_size()[0]];
	
	int b=1;
	int n;
	for(i=0;i<dim_-2;i++)b*=sizeLocal_[i];
	
	//calulate offset in dim-2
	for(i=0;i<parallel.grid_size()[1];i++)sizes1[i]=sizes_dim1[i] * b;
	for(n=0;n<parallel.grid_size()[1];n++)offset1[n]=0;
	
	for(n=1;n<parallel.grid_size()[1];n++)
	{
		for(i=0;i<n;i++)offset1[n]+=sizes1[i];
	}
	
	//calulate offset in dim-1
	for(i=0;i<parallel.grid_size()[0];i++)sizes0[i]=size_[dim_-2] * sizes_dim0[i] * b;
	for(n=0;n<parallel.grid_size()[0];n++)offset0[n]=0;
	
	for(n=1;n<parallel.grid_size()[0];n++)
	{
		for(i=0;i<n;i++)offset0[n]+=sizes0[i];
	}
	
	sitesSkip2d_ = offset0[parallel.grid_rank()[0]] + offset1[parallel.grid_rank()[1]];
	
	//Set status
	status_ = status_ | initialized;
	
	//Free memory
	delete[] sizes_dim0;
	delete[] sizes_dim1;
	
	
}
#ifdef FFT3D
void Lattice::initializeRealFFT(Lattice & lat_real, int halo)
{
	
	if(lat_real.dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	if(lat_real.size(0)!=lat_real.size(1) | lat_real.size(2)!=lat_real.size(1))
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	int lat_size[3];
	
	lat_size[0]=lat_real.size(0)/2+1;
	lat_size[1]=lat_real.size(0);
	lat_size[2]=lat_real.size(0);
	
	this->initialize(3, lat_size, halo);
}
void Lattice::initializeComplexFFT(Lattice & lat_real, int halo)
{
	
	if(lat_real.dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	if(lat_real.size(0)!=lat_real.size(1) | lat_real.size(2)!=lat_real.size(1))
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	int lat_size[3];
	
	lat_size[0]=lat_real.size(0);
	lat_size[1]=lat_real.size(0);
	lat_size[2]=lat_real.size(0);
	
	this->initialize(3, lat_size, halo);
}
#endif


void Lattice::save_arch(const string filename)
{
	fstream file;
	int p,i;
	
	for( p=0; p<parallel.size(); p++ ) 
	{
		if( parallel.rank()==p )
		{
			if( parallel.rank()==0)
			{
				//Truncate file if first process
				file.open(filename.c_str(), fstream::out | fstream::trunc);
			}
			else
			{
				//Append to file if another process
				file.open(filename.c_str(), fstream::out | fstream::app);
			}
			if(!file.is_open())
			{
				cerr<<"Latfield::Lattice::save_arch - Could not open file for writing"<<endl;
				cerr<<"Latfield::Lattice::save_arch - File: "<<filename<<endl;
				parallel.abortRequest();
			}
			
			if( parallel.rank()==0)
			{
				file<<"# Architerctur of the lattice"<<endl;
				file<<"# Number of dimension :"<<endl;
				file<<this->dim()<<endl;
				file<<"############################"<<endl;
			}
			file<<"############################"<<endl;
			file<<"#  Ranks of processor, world, dim-1, dim-2 :"<<endl;
			file<<parallel.rank()<<" "<<parallel.grid_rank()[0]<<" "<<parallel.grid_rank()[1] <<endl;
			file<<"#  Local size :"<<endl;
			for(i=0;i<this->dim();i++)
			{
				file<<this->sizeLocal()[i]<<endl;
			}
			file<<"############################"<<endl;
			
			
			file.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	arch_saved_=true;
	if(parallel.rank()==0)cout<<"Architecture saved in : "<<filename<<endl;
	
}

int Lattice::getRank(int* coord)
{
    int n,m;
    int temp;
    bool flag;
    
    temp=0;
    flag = true;
    for(n=0;flag;n++)
    {
        temp+= sizeLocalAllProcDim0_[n];
        if(temp>coord[dim_-1])flag=false;
    }
    n-=1;
    temp=0;
    flag = true;
    for(m=0;flag;m++)
    {
        temp+= sizeLocalAllProcDim1_[m];
        if(temp>coord[dim_-2])flag=false;
    }
    m-=1;
    
    
    return parallel.grid2world(n,m);
    
}
int Lattice::getRankDim0(int coord)
{
    int n;
    int temp;
    bool flag;
    
    temp=0;
    flag = true;
    for(n=0;flag;n++)
    {
        temp+= sizeLocalAllProcDim0_[n];
        if(temp>coord)flag=false;
    }
    return n-1;
}
int Lattice::getRankDim1(int coord)
{
    int m;
    int temp;
    bool flag;
    temp=0;
    flag = true;
    for(m=0;flag;m++)
    {
        temp+= sizeLocalAllProcDim1_[m];
        if(temp>coord)flag=false;
    }
    return  m-1;
}


//MISCELLANEOUS======================

bool Lattice::is_arch_saved() {return arch_saved_;};
int  Lattice::dim() { return dim_; };
int* Lattice::size() { return size_; };
int  Lattice::size(int i) { return size_[i]; }
long  Lattice::sites() { return sites_; }
long  Lattice::sitesGross() { return sitesGross_; }
int  Lattice::halo() { return halo_; }

int* Lattice::sizeLocal() { return sizeLocal_; };
int  Lattice::sizeLocal(int i) { return sizeLocal_[i]; }
long  Lattice::sitesLocal() { return sitesLocal_; }
long  Lattice::sitesLocalGross() { return sitesLocalGross_; }

long  Lattice::jump(int i) { return jump_[i]; }
long  Lattice::sitesSkip() { return sitesSkip_; }
long  Lattice::sitesSkip2d() { return sitesSkip2d_; }
long*  Lattice::coordSkip() { return coordSkip_; }
long Lattice::siteFirst() { return siteFirst_; }
long Lattice::siteLast() { return siteLast_; }

#endif

