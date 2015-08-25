
#ifndef LATFIELD2_PLANFFT_HPP
#define LATFIELD2_PLANFFT_HPP

/*! \file LATfield2_PlanFFT.hpp
 \brief LATfield2_PlanFFT.hpp contains the class PlanFFT definition.
 \author David Daverio
 */ 



#include "fftw3.h"

#ifdef SINGLE
#define MPI_DATA_PREC MPI_FLOAT
#endif

#ifndef SINGLE
#define MPI_DATA_PREC MPI_DOUBLE
#endif


const int FFT_FORWARD = FFTW_FORWARD;
const int FFT_BACKWARD = FFTW_BACKWARD;
const int FFT_IN_PLACE = 16;
const int FFT_OUT_OF_PLACE = -16;


/*! \class temporaryMemFFT
 \brief A class wich handle the additional memory needed by the class PlanFFT; Should never be used, internal temporary memroy handler, if real need, can be hacked, but could comflict with the PlanFFT!
 */
class temporaryMemFFT
	{
	public:
		temporaryMemFFT();
		~temporaryMemFFT();
		temporaryMemFFT(long size);
		
		int setTemp(long size);
		
#ifdef SINGLE
		fftwf_complex* temp1(){return temp1_;}
		fftwf_complex* temp2(){return temp2_;}
#endif
		
#ifndef SINGLE
        fftw_complex * temp1(){return temp1_;}
        fftw_complex * temp2(){return temp2_;}
#endif
		
	private:
#ifdef SINGLE
		fftwf_complex * temp1_;
		fftwf_complex * temp2_;
#endif
		
#ifndef SINGLE
		fftw_complex * temp1_;
		fftw_complex * temp2_;
#endif
		long allocated_; //number of variable stored (bit = allocated*sizeof(fftw(f)_complex))
	};
temporaryMemFFT::temporaryMemFFT()
{
	allocated_=0;
}
temporaryMemFFT::~temporaryMemFFT()
{
    if(allocated_!=0){
#ifdef SINGLE
	fftwf_free(temp1_);
	fftwf_free(temp2_);	
#endif
#ifndef SINGLE
	fftw_free(temp1_);
	fftw_free(temp2_);
#endif
	}
}

temporaryMemFFT::temporaryMemFFT(long size)
{
#ifdef SINGLE
	temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
	temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
#endif
	
#ifndef SINGLE
	temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif
	allocated_=size;
}
int temporaryMemFFT::setTemp(long size)
{
	if(size>allocated_)
	{
#ifdef SINGLE
        if(allocated_!=0){
		fftwf_free(temp1_);
		fftwf_free(temp2_);
        }
		temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		
		//debug cout<<"("<< parallel.grid_rank()[0]<< ","<<parallel.grid_rank()[1] <<"): called temporary resize. Old size: " <<allocated_<<" , new size: "<< size<<endl;
		
#endif
		
#ifndef SINGLE
        if(allocated_!=0){
		fftw_free(temp1_);
		fftw_free(temp2_);
        }
		temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif
		allocated_ = size ;
	

	}
	return 1;
}


//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;


/*! \class PlanFFT  
 
 \brief Class which handle Fourier transforms of fields (real/complex, single/double precision) on cubic lattices. (currently implemented only for 3d)
 
  
 This class allow to perform Fourier transform of real and complex fields. See the FFTs example to have have a short intro of the usage.  See the PoissonSolver example for more advanced usage (as linking several field to the same Fourier image)
 
 One should understand that first a plan is created then executed (in the FFTW fashion). The plan links two fields, one on Fourier space, one on real space. Boths field will be allocated by the planer when the plan is initilized. But need to be initialized before passing them to the planer. 
 
 One need to be carefull to corretly define the lattice and field.
 \sa void Lattice::initializeRealFFT(Lattice & lat_real, int halo);
 \sa void Lattice::initializeComplexFFT(Lattice & lat_real, int halo);
 \sa cKSite class
 \sa rKSite class
 
 */
template<class compType>
class PlanFFT
	{
	public:
        //! Constructor.
		PlanFFT();
		
#ifndef SINGLE
		
        /*!
         Constructor with initialization for complex to complex transform.
         
         \param rfield : real space field
         \param kfield : Fourier space field
         \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both Fourier and real space field point to the same data array.
         
         \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
         
         */ 
		PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
		/*!
         Initialization for complex to complex transform.
         
         \param rfield : real space field
         \param kfield : Fourier space field
         \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both Fourier and real space field point to the same data array.
         */
        void initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
		
		/*!
         Constructor with initialization for real to complex transform.
        
         \param rfield : real space field
         \param kfield : Fourier space field
         \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both Fourier and real space fields point to the same data array.
         
         \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
         */
		PlanFFT(Field<double>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
        /*!
         Initialization for real to complex transform.
         
         \param rfield : real space field
         \param kfield : Fourier space field
         \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both Fourier and real space fields point to the same data array.
         */
		void initialize(Field<double>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
		
		
		
#endif
		
#ifdef SINGLE		
		
		PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
		void initialize(Field<compType> * rfield,Field<compType> * kfield,const int mem_type = FFT_OUT_OF_PLACE);
		
		
		PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
		void initialize(Field<float>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
		
		
#endif
		
		
		/*!
         Execute the Fourier transform.
         
         \param fft_type : dirrection of the transform. Can be FFT_BACKWARD or FFT_FORWARD.
         */
		void execute(int fft_type);
        
	private:
		bool status_;
		bool type_;
		int mem_type_;
		
		static bool R2C;
		static bool C2C;
		
		static bool initialized;
		
		//data description variable, fftw plan, and temp
		
		int components_;
		int rSize_[3];
		int kSize_[3];
		int rJump_[3];
		int kJump_[3];
		int rSizeLocal_[3];
		int kSizeLocal_[3];
		int r2cSize_;  //of 0 dimension
		int r2cSizeLocal_; //of r2csize on proc_dim[1]
		int r2cSizeLocal_as_;
		int rHalo_;
		int kHalo_;
		
#ifdef SINGLE
		float * rData_; //pointer to start of data (halo skip)
		fftwf_complex * cData_; //pointer to start of data (halo skip)
		fftwf_complex * kData_; //pointer to start of data (halo skip)
		fftwf_complex * temp_;
		fftwf_complex * temp1_;//needed if field got more than 1 component
        
		
		fftwf_plan fPlan_i_;
		fftwf_plan fPlan_j_;
		fftwf_plan fPlan_k_;
		fftwf_plan fPlan_k_real_;
		
		fftwf_plan bPlan_i_;
		fftwf_plan bPlan_j_;
		fftwf_plan bPlan_j_real_;
		fftwf_plan bPlan_k_;
		
		//transpostion fonction
		
		// forward real to complex
		// first transopsition
		void transpose_0_2( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
		void transpose_0_2_last_proc( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
		void implement_local_0_last_proc( fftwf_complex * in, fftwf_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size);
		// second transposition
		void transpose_1_2(fftwf_complex * in , fftwf_complex * out  ,int dim_i,int dim_j ,int dim_k);
		//third transposition 
		void transpose_back_0_3(fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components,int comp);
		void implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k,int halo,int components,int comp);
		//backward real to complex
		void b_arrange_data_0(fftwf_complex *in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp);
		void b_transpose_back_0_1(fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size);
		void b_implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k);
        
		
#endif
#ifndef SINGLE
        
        
        double * rData_; //pointer to start of data (halo skip)
		fftw_complex * cData_; //pointer to start of data (halo skip)
		fftw_complex * kData_; //pointer to start of data (halo skip)
		fftw_complex * temp_;
		fftw_complex * temp1_;//needed if field got more than 1 component
        
		
		fftw_plan fPlan_i_;
		fftw_plan fPlan_j_;
		fftw_plan fPlan_k_;
		fftw_plan fPlan_k_real_;
		
		fftw_plan bPlan_i_;
		fftw_plan bPlan_j_;
		fftw_plan bPlan_j_real_;
		fftw_plan bPlan_k_;
		
		///transpostion fonction
		
		/// forward real to complex
		// first transopsition
		void transpose_0_2( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k);
		void transpose_0_2_last_proc( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k);
		void implement_local_0_last_proc( fftw_complex * in, fftw_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size);
		// second transposition
		void transpose_1_2(fftw_complex * in , fftw_complex * out  ,int dim_i,int dim_j ,int dim_k);
		//third transposition 
		void transpose_back_0_3(fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components,int comp);
		void implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k,int halo,int components,int comp);
		////backward real to complex
		void b_arrange_data_0(fftw_complex *in, fftw_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp);
		void b_transpose_back_0_1(fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size);
		void b_implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k);
        
        
		
		
#endif
		
	};

//constants
template<class compType>
bool PlanFFT<compType>::initialized = true;
template<class compType>
bool PlanFFT<compType>::R2C=false;
template<class compType>
bool PlanFFT<compType>::C2C=true;




template<class compType>
PlanFFT<compType>::PlanFFT()
{
	status_ = false;
}


#ifdef SINGLE



template<class compType>
PlanFFT<compType>::PlanFFT(Field<compType>*  rfield,Field<compType>* kfield,const int mem_type)
{
	status_ = false;
	initialize(rfield,kfield,mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<compType>*  rfield,Field<compType>*  kfield,const int mem_type )
{
	type_ = C2C;
	mem_type_=mem_type;
    
	//general variable
	if(rfield->components() != kfield->components())
	{
		cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
		cerr<<"Latfield2d::PlanFFT::initialize : coordinate and Fourier space fields have not the same number of components"<<endl;
		cerr<<"Latfield2d : Abort Process Requested"<<endl;
		
	}
	else components_ = rfield->components();
	
	for(int i = 0; i<3; i++)
	{
		rSize_[i]=rfield->lattice().size(i);
		kSize_[i]=kfield->lattice().size(i);
		rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
		kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
		rJump_[i]=rfield->lattice().jump(i);
		kJump_[i]=kfield->lattice().jump(i);
	}
	
	rHalo_ = rfield->lattice().halo();
	kHalo_ = kfield->lattice().halo();
    
	tempMemory.setTemp((long)(rSize_[0])  * (long)(rSizeLocal_[1]) * (long)(rSizeLocal_[2]));
	
	temp_  = tempMemory.temp1();
	temp1_ = tempMemory.temp2();
    
	if(rfield->lattice().dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	
	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
 		}
		parallel.abortForce();
	}
    
	//initialization of fftw plan
    
    
	//Forward plan
	fPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,cData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_FORWARD,FFTW_ESTIMATE);
	//Backward plan
    
	bPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,kData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	bPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE); 
	bPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE);
    
	//allocation of field
	
	long rfield_size = rfield->lattice().sitesLocalGross();
	long kfield_size = kfield->lattice().sitesLocalGross(); 
	
	if(mem_type_ == FFT_IN_PLACE)
	{
		if(rfield_size>=kfield_size)
		{		    
			rfield->alloc();
			kfield->data() = (Imag*)rfield->data();
		}
		else
		{
			kfield->alloc();
			rfield->data() = (Imag*)kfield->data();
		}
	}
	if(mem_type_ == FFT_OUT_OF_PLACE)
	{
		rfield->alloc();
		kfield->alloc();
	}
    
	//Pointer to data
	
	rData_ = (float*)rfield->data(); //to be sure that rData is instantiate !
	cData_ = (fftwf_complex*)rfield->data();
	cData_ += rfield->lattice().siteFirst()*components_;
	kData_ = (fftwf_complex*)kfield->data();
	kData_ += kfield->lattice().siteFirst()*components_;
}

template<class compType>
PlanFFT<compType>::PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type )
{
	status_ = false;
	initialize(rfield,kfield,mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<float>*  rfield,Field<compType>*   kfield,const int mem_type )
{
	type_ = R2C;
	mem_type_=mem_type;
	
	//general variable
	if(rfield->components() != kfield->components())
	{
		cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
		cerr<<"Latfield2d::PlanFFT::initialize : coordinate and Fourier space fields have not the same number of components"<<endl;
		cerr<<"Latfield2d : Abort Process Requested"<<endl;
		
	}
	else components_ = rfield->components();
	
	for(int i = 0; i<3; i++)
	{
		rSize_[i]=rfield->lattice().size(i);
		kSize_[i]=kfield->lattice().size(i);
		rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
		kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
		rJump_[i]=rfield->lattice().jump(i);
		kJump_[i]=kfield->lattice().jump(i);
	}
	r2cSize_ = rfield->lattice().size(0)/2 + 1;
	r2cSizeLocal_as_ = rfield->lattice().sizeLocal(0)/(2*parallel.grid_size()[1]);
	if(parallel.last_proc()[1]) r2cSizeLocal_ = r2cSizeLocal_as_ + 1;
	else r2cSizeLocal_ = r2cSizeLocal_as_;
	rHalo_ = rfield->lattice().halo();
	kHalo_ = kfield->lattice().halo();
	
	
	
    
	
	tempMemory.setTemp((r2cSize_+2) * (rSizeLocal_[1]+2) * (rSizeLocal_[2]+2));
    
    
	
	temp_  = tempMemory.temp1();
	temp1_ = tempMemory.temp2();
    
	
	
	
	if(rfield->lattice().dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	//COUT << rSize_[0] << "   "<< rSize_[1]<< "   "<< rSize_[2] <<endl;
	
	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	//create the fftw_plan
	
	fPlan_i_ = fftwf_plan_many_dft_r2c(1,&rSize_[0],rSizeLocal_[1] ,rData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_real_ =  fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);
	
	bPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],kSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
	bPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_BACKWARD,FFTW_ESTIMATE); 
	bPlan_j_real_ =  fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
	bPlan_i_ = fftwf_plan_many_dft_c2r(1,&rSize_[0],rSizeLocal_[1] ,temp1_,NULL,1, r2cSize_,rData_,NULL,components_,rJump_[1]*components_,FFTW_ESTIMATE);
	
	
	//allocation of field
	
	long rfield_size = rfield->lattice().sitesLocalGross();
	long kfield_size = kfield->lattice().sitesLocalGross()*2; //*2 for complex type
	
	if(mem_type_==FFT_IN_PLACE)
	{
		if(rfield_size>kfield_size)
		{
			rfield->alloc();
			kfield->data() = (Imag *)rfield->data();
		}
		else
		{
			kfield->alloc();
			rfield->data() = (float *)kfield->data();
		}
	}
	if(mem_type_ == FFT_OUT_OF_PLACE)
	{
		rfield->alloc();
		kfield->alloc();
	}
	
	
	//Pointer to data
	
	rData_ = rfield->data();
	rData_ += rfield->lattice().siteFirst()*components_;//usefull point !!
	cData_ = (fftwf_complex*)kfield->data(); //to be sure that cData is instantiate !
	kData_ = (fftwf_complex*)kfield->data();
	kData_ += kfield->lattice().siteFirst()*components_;
	
	
	
}

#endif

#ifndef SINGLE

template<class compType>
PlanFFT<compType>::PlanFFT(Field<compType>*  rfield,Field<compType>* kfield,const int mem_type)
{
	status_ = false;
	initialize(rfield,kfield,mem_type);
}

template<class compType>
void PlanFFT<compType>::initialize(Field<compType>*  rfield,Field<compType>*  kfield,const int mem_type )
{
	type_ = C2C;
	mem_type_=mem_type;
    
	//general variable
	if(rfield->components() != kfield->components())
	{
		cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
		cerr<<"Latfield2d::PlanFFT::initialize : coordinate and Fourier space fields have not the same number of components"<<endl;
		cerr<<"Latfield2d : Abort Process Requested"<<endl;
		
	}
	else components_ = rfield->components();
	
	for(int i = 0; i<3; i++)
	{
		rSize_[i]=rfield->lattice().size(i);
		kSize_[i]=kfield->lattice().size(i);
		rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
		kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
		rJump_[i]=rfield->lattice().jump(i);
		kJump_[i]=kfield->lattice().jump(i);
	}
	
	rHalo_ = rfield->lattice().halo();
	kHalo_ = kfield->lattice().halo();
    
	tempMemory.setTemp((long)(rSize_[0])  * (long)(rSizeLocal_[1]) * (long)(rSizeLocal_[2]));
	
	temp_  = tempMemory.temp1();
	temp1_ = tempMemory.temp2();
    
	if(rfield->lattice().dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	
	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
 		}
		parallel.abortForce();
	}
    
	//initialization of fftw plan
    
    
	//Forward plan
	fPlan_i_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,cData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_FORWARD,FFTW_ESTIMATE);
	//Backward plan
    
	bPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,kData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	bPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE); 
	bPlan_i_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE);
    
	//allocation of field
	
	long rfield_size = rfield->lattice().sitesLocalGross();
	long kfield_size = kfield->lattice().sitesLocalGross(); 
	
	if(mem_type_ == FFT_IN_PLACE)
	{
		if(rfield_size>=kfield_size)
		{		    
			rfield->alloc();
			kfield->data() = (Imag*)rfield->data();
		}
		else
		{
			kfield->alloc();
			rfield->data() = (Imag*)kfield->data();
		}
	}
	if(mem_type_ == FFT_OUT_OF_PLACE)
	{
		rfield->alloc();
		kfield->alloc();
	}
    
	//Pointer to data
	
	rData_ = (double*)rfield->data(); //to be sure that rData is instantiate !
	cData_ = (fftw_complex*)rfield->data();
	cData_ += rfield->lattice().siteFirst()*components_;
	kData_ = (fftw_complex*)kfield->data();
	kData_ += kfield->lattice().siteFirst()*components_;
	
}



template<class compType>
PlanFFT<compType>::PlanFFT(Field<double>* rfield, Field<compType>*  kfield,const int mem_type )
{
	status_ = false;
	initialize(rfield,kfield,mem_type);
}



template<class compType>
void PlanFFT<compType>::initialize(Field<double>*  rfield,Field<compType>*   kfield,const int mem_type )
{
	type_ = R2C;
	mem_type_=mem_type;
	
	//general variable
	if(rfield->components() != kfield->components())
	{
		cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
		cerr<<"Latfield2d::PlanFFT::initialize : coordinate and Fourier space fields have not the same number of components"<<endl;
		cerr<<"Latfield2d : Abort Process Requested"<<endl;
		
	}
	else components_ = rfield->components();
	
	for(int i = 0; i<3; i++)
	{
		rSize_[i]=rfield->lattice().size(i);
		kSize_[i]=kfield->lattice().size(i);
		rSizeLocal_[i]=rfield->lattice().sizeLocal(i);
		kSizeLocal_[i]=kfield->lattice().sizeLocal(i);
		rJump_[i]=rfield->lattice().jump(i);
		kJump_[i]=kfield->lattice().jump(i);
	}
	r2cSize_ = rfield->lattice().size(0)/2 + 1;
	r2cSizeLocal_as_ = rfield->lattice().sizeLocal(0)/(2*parallel.grid_size()[1]);
	if(parallel.last_proc()[1]) r2cSizeLocal_ = r2cSizeLocal_as_ + 1;
	else r2cSizeLocal_ = r2cSizeLocal_as_;
	rHalo_ = rfield->lattice().halo();
	kHalo_ = kfield->lattice().halo();
	
	
	
    
	
	tempMemory.setTemp((r2cSize_+2) * (rSizeLocal_[1]+2) * (rSizeLocal_[2]+2));
    
    
	
	temp_  = tempMemory.temp1();
	temp1_ = tempMemory.temp2();
    
	
	
	
	if(rfield->lattice().dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	//COUT << rSize_[0] << "   "<< rSize_[1]<< "   "<< rSize_[2] <<endl;
	
	if(rSize_[0]!=rSize_[1] | rSize_[1]!=rSize_[2])
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::PlanFFT::initialize : real lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	//create the fftw_plan
	
	
	
    
	fPlan_i_ = fftw_plan_many_dft_r2c(1,&rSize_[0],rSizeLocal_[1] ,rData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_FORWARD,FFTW_ESTIMATE); 
	fPlan_k_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);
	
	bPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],kSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
	bPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_BACKWARD,FFTW_ESTIMATE); 
	bPlan_j_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
	bPlan_i_ = fftw_plan_many_dft_c2r(1,&rSize_[0],rSizeLocal_[1] ,temp1_,NULL,1, r2cSize_,rData_,NULL,components_,rJump_[1]*components_,FFTW_ESTIMATE);
	
	//allocation of field
	
	long rfield_size = rfield->lattice().sitesLocalGross();
	long kfield_size = kfield->lattice().sitesLocalGross()*2; //*2 for complex type
	
	if(mem_type_==FFT_IN_PLACE)
	{
		if(rfield_size>kfield_size)
		{
			rfield->alloc();
			kfield->data() = (Imag *)rfield->data();
		}
		else
		{
			kfield->alloc();
			rfield->data() = (double *)kfield->data();
		}
	}
	if(mem_type_ == FFT_OUT_OF_PLACE)
	{
		rfield->alloc();
		kfield->alloc();
	}
	
	
	//Pointer to data
	
	rData_ = rfield->data();
	rData_ += rfield->lattice().siteFirst()*components_;//usefull point !!
	cData_ = (fftw_complex*)kfield->data(); //to be sure that cData is instantiate !
	kData_ = (fftw_complex*)kfield->data();
	kData_ += kfield->lattice().siteFirst()*components_;
	
	
	
}

#endif

template<class compType>
void PlanFFT<compType>::execute(int fft_type)
{
    temp_  = tempMemory.temp1();
	temp1_ = tempMemory.temp2();
    
	
    //#ifdef SINGLE
	if(type_ == R2C)
	{
		
		if(fft_type == FFT_FORWARD)
		{
			int i,j,k;
			int comp;
			int comm_rank;
			
#ifdef SINGLE
			float *p_in;
			fftwf_complex *p_out;		
#else
			double *p_in;
			fftw_complex *p_out;
#endif	
			
			for(comp=0;comp<components_;comp++)
			{
                
                
                //execute first dimension fft, prepar rData in temp_ to be send via AlltoAll + gather
                
#ifdef SINGLE
				for(int l = 0;l< rSizeLocal_[2] ;l++)
				{
					p_in = &rData_[rJump_[2]*l*components_ + comp];
					p_out = &temp_[l*rSizeLocal_[1]];
					fftwf_execute_dft_r2c(fPlan_i_,p_in,p_out);
				}
                
#else
				for(int l = 0;l< rSizeLocal_[2] ;l++)
				{
					p_in = &rData_[rJump_[2]*l*components_ + comp];
					p_out = &temp_[l*rSizeLocal_[1]];
					fftw_execute_dft_r2c(fPlan_i_,p_in,p_out);
				}
#endif
                /*
                 //verif step 1
                 for(int proc=0; proc<parallel.size(); proc++)
                 {
                 if(proc==parallel.rank())
                 {
                 for(k=0 ; k < rSizeLocal_[2] ;k++)
                 {
                 for(j=0 ; j < rSizeLocal_[1] ;j++)
                 {
                 for(i=0 ; i <  (r2cSize_) ;i++)
                 {
                 //cout << "("<<parallel.grid_rank()[0]<<","<< parallel.grid_rank()[1]<< "),("<< i <<","<< j + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0]<<" , "<<temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1]<<endl;
                 //temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0] = k + rSizeLocal_[2]*parallel.grid_rank()[0];
                 //temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1] = 0;//j + rSizeLocal_[1]*parallel.grid_rank()[1];
                 //cout <<rSizeLocal_[1]<<" "<< parallel.grid_rank()[1] <<"("<< i <<","<< j + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0]<<" , "<<temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1]<<endl;
                 
                 
                 }
                 }
                 }
                 }
                 MPI_Barrier(parallel.lat_world_comm());
                 
                 }	*/		
				
				
				MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Gather(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
				
				
				
				
				if(parallel.last_proc()[1])
				{
					for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2_last_proc(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
					implement_local_0_last_proc(&temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],temp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1]);
				}
				else for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
				
#ifdef SINGLE
				fftwf_execute(fPlan_j_);
#else
				fftw_execute(fPlan_j_);
#endif
				//verif step 2
				/*
                 for(k=0 ; k < rSizeLocal_[2] ;k++)
                 {
                 for(j=0 ; j < rSize_[1] ;j++)
                 {
                 for(i=0 ; i <  (r2cSizeLocal_) ;i++)
                 {
                 
                 //cout << "("<< i + r2cSizeLocal_as_*parallel.grid_rank()[1] <<","<< j  <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[i+r2cSizeLocal_*(k+rSizeLocal_[2]*j)][0]<<" , "<<temp_[i+r2cSizeLocal_*(k+rSizeLocal_[2]*j)][1]<<endl;							
                 temp_[i+r2cSizeLocal_*(k+rSizeLocal_[2]*j)][1] = k + rSizeLocal_[2]*parallel.grid_rank()[0] ;
                 temp_[i+r2cSizeLocal_*(k+rSizeLocal_[2]*j)][0] = i + r2cSizeLocal_as_*parallel.grid_rank()[1];
                 
                 }
                 }
                 }*/
				
				
                
				
				MPI_Barrier(parallel.lat_world_comm());
				MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
				MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
				
				
				for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2]);
                
				/*
                 //verif transpos
                 for(int proc=0; proc<parallel.size(); proc++)
                 {
                 if(proc==parallel.rank())
                 {
                 for(k=0 ; k < rSize_[2] ;k++)
                 {
                 for(j=0 ; j < rSizeLocal_[2] ;j++)
                 {
                 for(i=0 ; i < r2cSizeLocal_ ;i++)
                 {
                 cout  << "("<<parallel.grid_rank()[0]<<","<< parallel.grid_rank()[1]<< "),("<< i + r2cSizeLocal_as_*parallel.grid_rank()[1] <<","<< j+rSizeLocal_[2]*parallel.grid_rank()[0]  <<","<< k<<")"<< temp_[i+r2cSizeLocal_*(j+rSizeLocal_[2]*k)][0]<<" , "<<temp_[i+r2cSizeLocal_*(j+rSizeLocal_[2]*k)][1]<<endl;
                 temp_[i+r2cSizeLocal_*(j+rSizeLocal_[2]*k)][0] = 1;
                 temp_[i+r2cSizeLocal_*(j+rSizeLocal_[2]*k)][1] = 0;
                 }
                 }
                 }
                 }
                 MPI_Barrier(parallel.lat_world_comm());
                 }*/
                
                
#ifdef SINGLE
				for(int l=0;l<rSizeLocal_[2];l++)
                {
				    fftwf_execute_dft(fPlan_k_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
                }
                
				if(parallel.last_proc()[1])fftwf_execute(fPlan_k_real_);
#else 
				for(int l=0;l<rSizeLocal_[2];l++)
                {
				    fftw_execute_dft(fPlan_k_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
                }
                
				if(parallel.last_proc()[1])fftw_execute(fPlan_k_real_);
                
#endif		
                
				
				
				MPI_Alltoall(temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Scatter(&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
                
                
				transpose_back_0_3(temp_, kData_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1],kHalo_,components_,comp);
				implement_0(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], kData_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1],kHalo_,components_,comp);
				
				
			}
			
			
		}
		if(fft_type == FFT_BACKWARD)
		{	
			int i,j,k,comp; 
			int comm_rank;
            
#ifdef SINGLE	
			float *p_out;
			fftwf_complex *p_in;
#else
			double *p_out;
			fftw_complex *p_in;
#endif
			
			for(comp=0;comp<components_;comp++)
			{
				b_arrange_data_0(kData_, temp_,kSizeLocal_[0],kSizeLocal_[1] ,kSizeLocal_[2], kHalo_, components_, comp);
				MPI_Barrier(parallel.lat_world_comm());
				
				
				
				
                
				MPI_Alltoall(temp_,2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Gather(&temp_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
                MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
                
                
				
				if(parallel.last_proc()[1])
				{
					for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2_last_proc(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
					implement_local_0_last_proc(&temp1_[rSize_[0]/2*rSizeLocal_[1]*rSizeLocal_[2]],temp_,rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_,parallel.grid_size()[1]);
				}
				else for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_],rSizeLocal_[1],rSizeLocal_[2],r2cSizeLocal_as_);
				
				
				
                
				
				
#ifdef SINGLE				
				fftwf_execute(bPlan_k_);
#else
				fftw_execute(bPlan_k_);
#endif				
                
				
				MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
				MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);
                
                
				for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2]);
				
                
                
#ifdef SINGLE				
                
				for(int l=0;l<rSizeLocal_[2];l++)
                {
				    fftwf_execute_dft(bPlan_j_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
                }
                
				if(parallel.last_proc()[1])fftwf_execute(bPlan_j_real_);
				
#else
				for(int l=0;l<rSizeLocal_[2];l++)
                {
				    fftw_execute_dft(bPlan_j_,&temp_[l*r2cSizeLocal_],&temp1_[l*r2cSizeLocal_as_]);
                }
                
				if(parallel.last_proc()[1])fftw_execute(bPlan_j_real_);
#endif				
                
                
				
                
                
				MPI_Alltoall(temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*r2cSizeLocal_as_, MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Scatter(&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]][0], 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC, &temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]][0] , 2*rSizeLocal_[1]*rSizeLocal_[2], MPI_DATA_PREC ,parallel.grid_size()[1]-1, parallel.dim1_comm()[parallel.grid_rank()[0]]);
				MPI_Barrier(parallel.dim1_comm()[parallel.grid_rank()[0]]);
                
                
				
				b_transpose_back_0_1(temp_, temp1_,r2cSize_,r2cSizeLocal_as_,rSizeLocal_[2],rSizeLocal_[1],parallel.grid_size()[1]);
				b_implement_0(&temp_[(r2cSize_-1)*rSizeLocal_[1]*rSizeLocal_[2]], temp1_,r2cSize_,rSizeLocal_[2],rSizeLocal_[1]);
				
                
#ifdef SINGLE
				for(int l = 0;l< rSizeLocal_[2] ;l++)
				{
					p_in = &temp1_[ l*r2cSize_*rSizeLocal_[1] ];
					p_out = &rData_[l*rJump_[2]*components_ + comp];
					fftwf_execute_dft_c2r(bPlan_i_,p_in,p_out);
				}
#else
				for(int l = 0;l< rSizeLocal_[2] ;l++)
				{
					p_in = &temp1_[ l*r2cSize_*rSizeLocal_[1] ];
					p_out = &rData_[l*rJump_[2]*components_ + comp];
					fftw_execute_dft_c2r(bPlan_i_,p_in,p_out);
				}
#endif
				
				
			}
			
			
            
		}
		
		
	}
	if(type_ == C2C)
	{

        if(fft_type == FFT_FORWARD)
        {
            int i,j,k;
            int comp;
            int comm_rank;
            
#ifdef SINGLE
            fftwf_complex *p_in;
            fftwf_complex *p_out;		
#else
            fftw_complex *p_in;
            fftw_complex *p_out;
#endif	
            
            for(comp=0;comp<components_;comp++)
            {
                for(int l = 0;l< rSizeLocal_[2] ;l++)
                {
                    
                    p_in =  &cData_[rJump_[2]*l*components_ + comp];
                    p_out = &temp_[l*rSizeLocal_[1]];
#ifdef SINGLE
                    fftwf_execute_dft(fPlan_i_,p_in,p_out);
                    
#else
                    fftw_execute_dft(fPlan_i_,p_in,p_out);
#endif
                }
                
                
                /*
				//verif step 1
				for(int proc=0; proc<parallel.size(); proc++)
				{
					if(proc==parallel.rank())
					{
                        //cout << "proc number : "<<  proc <<endl;
                        for(k=0 ; k < rSizeLocal_[2] ;k++)
					    {
                            for(j=0 ; j < rSizeLocal_[1] ;j++)
                            {
                                for(i=0 ; i <  rSize_[0] ;i++)
                                {
                                    //cout << "("<<parallel.grid_rank()[0]<<","<< parallel.grid_rank()[1]<< "),("<< i <<","<< j + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0]<<" , "<<temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1]<<endl;
                                    //temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0] = k + rSizeLocal_[2]*parallel.grid_rank()[0];
                                    //temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1] = 0;//j + rSizeLocal_[1]*parallel.grid_rank()[1];
                                    //cout <<rSizeLocal_[1]<<" "<< parallel.grid_rank()[1] <<"("<< i <<","<< j + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][0]<<" , "<<temp_[j+rSizeLocal_[1]*(k+rSizeLocal_[2] *i)][1]<<endl;
                                    
                                    
                                }
                            }
					    }
					}
					MPI_Barrier(parallel.lat_world_comm());
                    
				}			
				*/
                
                
                
                MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
                
                for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[1]);
                
#ifdef SINGLE
                fftwf_execute(fPlan_j_);
#else
                fftw_execute(fPlan_j_);
#endif
                
                
                //verif step 2
                /*
                 for(k=0 ; k < rSizeLocal_[2] ;k++)
                 { 
                 for(j=0 ; j < rSize_[1] ;j++)
                 {
                 for(i=0 ; i <  (rSizeLocal_[1]) ;i++)
                 {
                 
                 //cout << "("<< i + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< j  <<","<< k+rSizeLocal_[2]*parallel.grid_rank()[0]<<")"<< temp_[i+rSizeLocal_[1]*(k+rSizeLocal_[2]*j)][0]<<" , "<<temp_[i+rSizeLocal_[1]*(k+rSizeLocal_[2]*j)][1]<<endl;							
                 //temp_[i+rSizeLocal_[1]*(k+rSizeLocal_[2]*j)][1] = k + rSizeLocal_[2]*parallel.grid_rank()[0] ;
                 //temp_[i+rSizeLocal_[1]*(k+rSizeLocal_[2]*j)][0] = i + rSizeLocal_[1]*parallel.grid_rank()[1];
                 
                 }
                 }
                 } 
                 */
                
                
                
                MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
                
                for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);
				
                
                
                
                /*
                 //verif transpos
                 for(int proc=0; proc<parallel.size(); proc++)
                 {
                 if(proc==parallel.rank())
                 {
                 for(k=0 ; k < rSize_[2] ;k++)
                 {
                 for(j=0 ; j < rSizeLocal_[2] ;j++)
                 {
                 for(i=0 ; i < rSizeLocal_[1] ;i++)
                 {
                 //cout  << "("<<parallel.grid_rank()[0]<<","<< parallel.grid_rank()[1]<< "),("<< i + rSizeLocal_[1]*parallel.grid_rank()[1] <<","<< j+rSizeLocal_[2]*parallel.grid_rank()[0]  <<","<< k<<")"<< temp_[i+rSizeLocal_[1]*(j+rSizeLocal_[2]*k)][0]<<" , "<<temp_[i+rSizeLocal_[1]*(j+rSizeLocal_[2]*k)][1]<<endl;
                 //temp_[i+rSizeLocal_[1]*(j+rSizeLocal_[2]*k)][0] = j+rSizeLocal_[2]*parallel.grid_rank()[0] ;
                 //temp_[i+rSizeLocal_[1]*(j+rSizeLocal_[2]*k)][1] = 0;
                 }
                 }
                 }
                 }
                 MPI_Barrier(parallel.lat_world_comm());
                 }
                 */
                
                
                for(int l = 0;l< rSizeLocal_[2] ;l++)
                {
                    
                    p_in = &temp_[l*rSizeLocal_[1]];
                    p_out =  &kData_[rJump_[2]*l*components_ + comp];
                    
#ifdef SINGLE
                    fftwf_execute_dft(fPlan_k_,p_in,p_out);    
#else
                    fftw_execute_dft(fPlan_k_,p_in,p_out);
#endif
                }
                
                
            }
            
        }
        if(fft_type == FFT_BACKWARD)
        {
            
            int i,j,k;
            int comp;
            int comm_rank;
            
#ifdef SINGLE
            fftwf_complex *p_in;
            fftwf_complex *p_out;		
#else
            fftw_complex *p_in;
            fftw_complex *p_out;
#endif	
            
            //STEP 1 : SAME AS STEP ONE OF FORWARD
            
            for(comp=0;comp<components_;comp++)
            {
                for(int l = 0;l< rSizeLocal_[2] ;l++)
                {
                    
                    p_in =  &kData_[rJump_[2]*l*components_ + comp];
                    p_out = &temp_[l*rSizeLocal_[1]];
#ifdef SINGLE
                    fftwf_execute_dft(bPlan_k_,p_in,p_out);
                    
#else
                    fftw_execute_dft(bPlan_k_,p_in,p_out);
#endif
                }
                
                
                //step 2 : same as step 4 of forward
                
                
                MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
                
                for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);
                
                
#ifdef SINGLE
                fftwf_execute(bPlan_j_);
#else
                fftw_execute(bPlan_j_);
#endif
                
                
                MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);
                
                for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[1]);
                
                
                for(int l = 0;l< rSizeLocal_[2] ;l++)
                {
                    
                    p_in = &temp_[l*rSizeLocal_[1]];
                    p_out =  &cData_[rJump_[2]*l*components_ + comp];
                    
#ifdef SINGLE
                    fftwf_execute_dft(bPlan_i_,p_in,p_out);    
#else
                    fftw_execute_dft(bPlan_i_,p_in,p_out);
#endif
                }
            }
            
        }
	
	}	
    
	
    
}

//transposition function
#ifdef SINGLE

/////
template<class compType>
void PlanFFT<compType>::transpose_0_2( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k)
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				
				out[k+dim_k*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[k+dim_k*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
				
			}
		}
	}
	
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k)
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[k+(dim_k+1)*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[k+(dim_k+1)*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( fftwf_complex * in, fftwf_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size)
{
	int i_in,i_out,j,rank;
	for(i_in=0;i_in<proc_dim_i;i_in++)
	{
		for(j=0;j<proc_dim_j;j++)
		{
			for(rank=0;rank<proc_size;rank++)
			{
				i_out=i_in+rank*proc_dim_i;
				out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][0]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][0];
				out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][1]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][1];
			}
		}
	}
}


template<class compType>
void PlanFFT<compType>::transpose_1_2(fftwf_complex * in , fftwf_complex * out ,int dim_i,int dim_j ,int dim_k )
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[i+dim_i*(k+j*dim_k)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[i+dim_i*(k+j*dim_k)][1]=in[i+dim_i*(j+k*dim_j)][1];
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::transpose_back_0_3( fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components, int comp)
{
	int i,j,k,l, i_t, j_t, k_t;
	int r2c_halo = r2c + 2*halo;
	int local_size_k_halo = local_size_k + 2*halo;
	for (i=0;i<local_r2c;i++)
	{
		for(k=0;k<local_size_k;k++)
		{
			for(j=0;j<local_size_j;j++)
			{
				for(l=0;l<proc_size;l++)
				{
					i_t = i + l*local_r2c;
					j_t = j ;
					k_t = k ;
					out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
					out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
				}
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k, int halo,int components, int comp)
{
	int i,j,k;
	i=r2c_size-1;
	int r2c_halo = r2c_size + 2*halo;
	int local_size_k_halo = local_size_k + 2*halo;
	
	for(j=0;j<local_size_j;j++)
	{
		for(k=0;k<local_size_k;k++)
		{
			out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][0]=in[j + local_size_j *k][0];
			out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][1]=in[j + local_size_j *k][1];
		}
	}
	
}


template<class compType>
void PlanFFT<compType>::b_arrange_data_0(fftwf_complex *in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp)
{
	int i,j,k;
	int jump_i=(dim_i+ 2 *khalo);
	int jump_j=dim_j+ 2 *khalo;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[j + dim_j * (k + dim_k * i)][0]=in[comp+components*(i + jump_i * (j + jump_j*k))][0];
				out[j + dim_j * (k + dim_k * i)][1]=in[comp+components*(i + jump_i * (j + jump_j*k))][1];
			}
		}
	}
	
}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size)
{
	int i,j,k,l, i_t, j_t, k_t;
	
	for (i=0;i<local_r2c;i++)
	{
		for(k=0;k<local_size_k;k++)
		{
			for(j=0;j<local_size_j;j++)
			{
				for(l=0;l<proc_size;l++)
				{
					i_t = i + l*local_r2c;
					j_t = j ;
					k_t = k ;
					out[i_t + r2c * (k_t + local_size_k * j_t)][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
					out[i_t + r2c * (k_t + local_size_k * j_t)][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
				}
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::b_implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k)
{
	int i,j,k;
	i=r2c_size-1;
	
	
	for(j=0;j<local_size_j;j++)
	{
		for(k=0;k<local_size_k;k++)
		{
			out[i + r2c_size * (k + local_size_k *j)][0]=in[j + local_size_j *k][0];
			out[i + r2c_size * (k + local_size_k *j)][1]=in[j + local_size_j *k][1];
		}
	}
	
}


#endif

#ifndef SINGLE

/////
template<class compType>
void PlanFFT<compType>::transpose_0_2( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k)
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				
				out[k+dim_k*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[k+dim_k*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
				
			}
		}
	}
	
}

template<class compType>
void PlanFFT<compType>::transpose_0_2_last_proc( fftw_complex * in, fftw_complex * out,int dim_i,int dim_j ,int dim_k)
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[k+(dim_k+1)*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[k+(dim_k+1)*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::implement_local_0_last_proc( fftw_complex * in, fftw_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size)
{
	int i_in,i_out,j,rank;
	for(i_in=0;i_in<proc_dim_i;i_in++)
	{
		for(j=0;j<proc_dim_j;j++)
		{
			for(rank=0;rank<proc_size;rank++)
			{
				i_out=i_in+rank*proc_dim_i;
				out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][0]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][0];
				out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][1]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][1];
			}
		}
	}
}


template<class compType>
void PlanFFT<compType>::transpose_1_2(fftw_complex * in , fftw_complex * out ,int dim_i,int dim_j ,int dim_k )
{
	int i,j,k;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[i+dim_i*(k+j*dim_k)][0]=in[i+dim_i*(j+k*dim_j)][0];
				out[i+dim_i*(k+j*dim_k)][1]=in[i+dim_i*(j+k*dim_j)][1];
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::transpose_back_0_3( fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components, int comp)
{
	int i,j,k,l, i_t, j_t, k_t;
	int r2c_halo = r2c + 2*halo;
	int local_size_k_halo = local_size_k + 2*halo;
	for (i=0;i<local_r2c;i++)
	{
		for(k=0;k<local_size_k;k++)
		{
			for(j=0;j<local_size_j;j++)
			{
				for(l=0;l<proc_size;l++)
				{
					i_t = i + l*local_r2c;
					j_t = j ;
					k_t = k ;
					out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
					out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
				}
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k, int halo,int components, int comp)
{
	int i,j,k;
	i=r2c_size-1;
	int r2c_halo = r2c_size + 2*halo;
	int local_size_k_halo = local_size_k + 2*halo;
	
	for(j=0;j<local_size_j;j++)
	{
		for(k=0;k<local_size_k;k++)
		{
			out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][0]=in[j + local_size_j *k][0];
			out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][1]=in[j + local_size_j *k][1];
		}
	}
	
}


template<class compType>
void PlanFFT<compType>::b_arrange_data_0(fftw_complex *in, fftw_complex * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp)
{
	int i,j,k;
	int jump_i=(dim_i+ 2 *khalo);
	int jump_j=dim_j+ 2 *khalo;
	for(i=0;i<dim_i;i++)
	{
		for(j=0;j<dim_j;j++)
		{
			for(k=0;k<dim_k;k++)
			{
				out[j + dim_j * (k + dim_k * i)][0]=in[comp+components*(i + jump_i * (j + jump_j*k))][0];
				out[j + dim_j * (k + dim_k * i)][1]=in[comp+components*(i + jump_i * (j + jump_j*k))][1];
			}
		}
	}
	
}

template<class compType>
void PlanFFT<compType>::b_transpose_back_0_1( fftw_complex * in, fftw_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size)
{
	int i,j,k,l, i_t, j_t, k_t;
	
	for (i=0;i<local_r2c;i++)
	{
		for(k=0;k<local_size_k;k++)
		{
			for(j=0;j<local_size_j;j++)
			{
				for(l=0;l<proc_size;l++)
				{
					i_t = i + l*local_r2c;
					j_t = j ;
					k_t = k ;
					out[i_t + r2c * (k_t + local_size_k * j_t)][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
					out[i_t + r2c * (k_t + local_size_k * j_t)][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
				}
			}
		}
	}
}

template<class compType>
void PlanFFT<compType>::b_implement_0(fftw_complex * in, fftw_complex * out,int r2c_size,int local_size_j,int local_size_k)
{
	int i,j,k;
	i=r2c_size-1;
	
	
	for(j=0;j<local_size_j;j++)
	{
		for(k=0;k<local_size_k;k++)
		{
			out[i + r2c_size * (k + local_size_k *j)][0]=in[j + local_size_j *k][0];
			out[i + r2c_size * (k + local_size_k *j)][1]=in[j + local_size_j *k][1];
		}
	}
	
}
#endif


#endif
