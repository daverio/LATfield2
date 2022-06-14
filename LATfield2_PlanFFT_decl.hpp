#ifndef LATFIELD2_PLANFFT_DECL_HPP
#define LATFIELD2_PLANFFT_DECL_HPP

#ifdef SINGLE
#define MPI_DATA_PREC MPI_FLOAT
#endif

#ifndef SINGLE
#define MPI_DATA_PREC MPI_DOUBLE
#endif

#ifndef NULLFFTWPLAN
#ifndef SINGLE
#define NULLFFTWPLAN static_cast<fftw_plan>(NULL)
#else
#define NULLFFTWPLAN static_cast<fftwf_plan>(NULL)
#endif

#endif


/* these are not the actual compiler variables, these are
   only the declaration that they exist somewhere at all. */

extern  const int FFT_FORWARD;
extern  const int FFT_BACKWARD;
extern  const int FFT_IN_PLACE;
extern  const int FFT_OUT_OF_PLACE;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*! \class temporaryMemFFT
 \brief A class wich handle the additional memory needed by the class PlanFFT_CPU; No documentation!
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


extern  temporaryMemFFT tempMemory;

/*! \class PlanFFT

 \brief Class which handle fourier transforms of fields on 3d cubic lattices.
 This class allow to perform fourier transform of real and complex fields. See poissonSolver example to have have a short intro of usage.

 One should understand that first a plan is created then execute (in the FFTW fashion). The plan link to fields, one on fourier space, one on real space. Both field will be allocated by the planer. But need to be initialized.

 One need to be carefull to corretly define the lattice and field.
 \sa void Lattice::initializeRealFFT(Lattice & lat_real, int halo);
 \sa void Lattice::initializeComplexFFT(Lattice & lat_real, int halo);
 For more detail see the QuickStart guide.
 */
template<class compType>
class PlanFFT
{
public:
  //! Constructor.
  PlanFFT();

  //! Destructor.
  ~PlanFFT();

#ifndef SINGLE
  /*!
   Constructor with initialization for complex to complex tranform.
   \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
  /*!
   initialization for complex to complex tranform.
   For more detail see the QuickStart guide.
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  void initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);

  /*!
   Constructor with initialization for real to complex tranform.
   \sa initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  PlanFFT(Field<double>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
  /*!
   initialization for real to complex tranform.
   For more detail see the QuickStart guide.
   \param rfield : real space field
   \param kfield : fourier space field
   \param mem_type : memory type (FFT_OUT_OF_PLACE or FFT_IN_PLACE). In place mean that both fourier and real space field point to the same data array.
   */
  void initialize(Field<double>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);



#endif

#ifdef SINGLE

  PlanFFT(Field<compType>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
  void initialize(Field<compType>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);


  PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type = FFT_OUT_OF_PLACE);
  void initialize(Field<float>*  rfield,Field<compType>*   kfield,const int mem_type = FFT_OUT_OF_PLACE);


#endif



  void execute(int fft_type);

private:
  void PrintPlans() {
#ifndef SINGLE
    std::cout << fPlan_i_ << " "; fftw_print_plan(fPlan_i_); std::cout << std::endl;
    std::cout << fPlan_j_ << " "; fftw_print_plan(fPlan_j_); std::cout << std::endl;
    std::cout << fPlan_k_ << " "; fftw_print_plan(fPlan_k_); std::cout << std::endl;
    std::cout << fPlan_k_real_ << " "; fftw_print_plan(fPlan_k_real_); std::cout << std::endl;
    std::cout << bPlan_k_ << " "; fftw_print_plan(bPlan_k_); std::cout << std::endl;
    std::cout << bPlan_j_ << " "; fftw_print_plan(bPlan_j_); std::cout << std::endl;
    std::cout << bPlan_j_real_ << " "; fftw_print_plan(bPlan_j_real_); std::cout << std::endl;
    std::cout << bPlan_i_ << " "; fftw_print_plan(bPlan_i_); std::cout << std::endl;
#else
    std::cout << fPlan_i_ << " "; fftwf_print_plan(fPlan_i_); std::cout << std::endl;
    std::cout << fPlan_j_ << " "; fftwf_print_plan(fPlan_j_); std::cout << std::endl;
    std::cout << fPlan_k_ << " "; fftwf_print_plan(fPlan_k_); std::cout << std::endl;
    std::cout << fPlan_k_real_ << " "; fftwf_print_plan(fPlan_k_real_); std::cout << std::endl;
    std::cout << bPlan_k_ << " "; fftwf_print_plan(bPlan_k_); std::cout << std::endl;
    std::cout << bPlan_j_ << " "; fftwf_print_plan(bPlan_j_); std::cout << std::endl;
    std::cout << bPlan_j_real_ << " "; fftwf_print_plan(bPlan_j_real_); std::cout << std::endl;
    std::cout << bPlan_i_ << " "; fftwf_print_plan(bPlan_i_); std::cout << std::endl;
#endif


}

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

  ///transpostion fonction

  /// forward real to complex
  // first transopsition
  void transpose_0_2( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
  void transpose_0_2_last_proc( fftwf_complex * in, fftwf_complex * out,int dim_i,int dim_j ,int dim_k);
  void implement_local_0_last_proc( fftwf_complex * in, fftwf_complex * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size);
  // second transposition
  void transpose_1_2(fftwf_complex * in , fftwf_complex * out  ,int dim_i,int dim_j ,int dim_k);
  //third transposition
  void transpose_back_0_3(fftwf_complex * in, fftwf_complex * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components,int comp);
  void implement_0(fftwf_complex * in, fftwf_complex * out,int r2c_size,int local_size_j,int local_size_k,int halo,int components,int comp);
  ////backward real to complex
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
PlanFFT<compType>::~PlanFFT() {
#ifndef SINGLE
  //if (fPlan_i_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_i_); }
  //if (fPlan_j_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_j_); }
  //if (fPlan_k_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_k_); }
  //if (fPlan_k_real_ != NULLFFTWPLAN) { fftw_destroy_plan(fPlan_k_real_); }
  //if (bPlan_k_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_k_); }
  //if (bPlan_j_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_j_); }
  //if (bPlan_j_real_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_j_real_); }
  //if (bPlan_i_ != NULLFFTWPLAN) { fftw_destroy_plan(bPlan_i_); }
#else
  //if (fPlan_i_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_i_); }
  //if (fPlan_j_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_j_); }
  //if (fPlan_k_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_k_); }
  //if (fPlan_k_real_ != NULLFFTWPLAN) { fftwf_destroy_plan(fPlan_k_real_); }
  //if (bPlan_k_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_k_); }
  //if (bPlan_j_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_j_); }
  //if (bPlan_j_real_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_j_real_); }
  //if (bPlan_i_ != NULLFFTWPLAN) { fftwf_destroy_plan(bPlan_i_); }
#endif



}


template<class compType>
PlanFFT<compType>::PlanFFT() :
fPlan_i_(NULLFFTWPLAN),
fPlan_j_(NULLFFTWPLAN),
fPlan_k_(NULLFFTWPLAN),
fPlan_k_real_(NULLFFTWPLAN),
bPlan_k_(NULLFFTWPLAN),
bPlan_j_(NULLFFTWPLAN),
bPlan_j_real_(NULLFFTWPLAN),
bPlan_i_(NULLFFTWPLAN)
{
  status_ = false;
}


#ifdef SINGLE



template<class compType>
PlanFFT<compType>::PlanFFT(Field<compType>*  rfield,Field<compType>* kfield,const int mem_type) : PlanFFT()
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

  COUT<<"INITIALIZING COMPLEX FFT"<<endl;

  if(rfield->components() != kfield->components())
  {
    cerr<<"Latfield2d::PlanFFT::initialize : fft curently work only for fields with same number of components"<<endl;
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
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
  r2cSize_ = 0;
  r2cSizeLocal_as_ = 0;
  r2cSizeLocal_ = 0;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();

  /////from latfield2d_IO
  tempMemory.setTemp((long)(rSize_[0]+2)  * (long)(rSizeLocal_[1]+2) * (long)(rSizeLocal_[2]+2));

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

    //Pointer to data

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

  	rData_ = (float*)rfield->data(); //to be sure that rData is instantiate !
  	cData_ = (fftwf_complex*)rfield->data();
  	cData_ += rfield->lattice().siteFirst()*components_;
  	kData_ = (fftwf_complex*)kfield->data();
  	kData_ += kfield->lattice().siteFirst()*components_;

  	//Forward plan
  	fPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,cData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	fPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);
  	fPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,kData_,NULL,components_, rJump_[1]*components_,FFTW_FORWARD,FFTW_ESTIMATE);
  	//Backward plan

  	bPlan_k_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,kData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  	bPlan_j_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1]*rSizeLocal_[2],temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  	bPlan_i_ = fftwf_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,cData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE| FFTW_PRESERVE_INPUT);

  	//allocation of field






  ///end of from latfield2d_IO



}

template<class compType>
PlanFFT<compType>::PlanFFT(Field<float>* rfield, Field<compType>*  kfield,const int mem_type )  : PlanFFT()
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
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
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

//  PrintPlans();

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
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
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
  r2cSize_ = 0;
  r2cSizeLocal_as_ = 0;
  r2cSizeLocal_ = 0;
  rHalo_ = rfield->lattice().halo();
  kHalo_ = kfield->lattice().halo();

  /////from latfield2d_IO
  tempMemory.setTemp((long)(rSize_[0]+10)  * (long)(rSizeLocal_[1]+10) * (long)(rSizeLocal_[2]+10));

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
  	bPlan_i_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[1] ,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,cData_,NULL,components_, rJump_[1]*components_,FFTW_BACKWARD,FFTW_ESTIMATE);

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
    cerr<<"Latfield2d::PlanFFT::initialize : coordinate and fourier space fields have not the same number of components"<<endl;
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

  fPlan_i_ = fftw_plan_many_dft_r2c(1,&rSize_[0],rSizeLocal_[1] ,rData_,NULL,components_, rJump_[1]*components_,temp_,NULL,rSizeLocal_[1]*rSizeLocal_[2],1,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  fPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  fPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_FORWARD,FFTW_ESTIMATE);
  fPlan_k_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_FORWARD,FFTW_ESTIMATE);

  bPlan_k_ = fftw_plan_many_dft(1,&rSize_[0],kSizeLocal_[2]*r2cSizeLocal_,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,temp_,NULL,kSizeLocal_[2]*r2cSizeLocal_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  bPlan_j_ = fftw_plan_many_dft(1,&rSize_[0],r2cSizeLocal_as_,temp_,NULL,rSizeLocal_[2]*r2cSizeLocal_,1,temp1_,NULL,rSizeLocal_[2]*r2cSizeLocal_as_,1,FFTW_BACKWARD,FFTW_ESTIMATE);
  bPlan_j_real_ =  fftw_plan_many_dft(1,&rSize_[0],rSizeLocal_[2],&temp_[r2cSizeLocal_as_],NULL,rSizeLocal_[2]*r2cSizeLocal_,r2cSizeLocal_,&temp1_[r2cSizeLocal_as_*rSizeLocal_[2]*rSize_[0]],NULL,rSizeLocal_[2],1,FFTW_BACKWARD,FFTW_ESTIMATE);
  bPlan_i_ = fftw_plan_many_dft_c2r(1,&rSize_[0],rSizeLocal_[1] ,temp1_,NULL,1, r2cSize_,rData_,NULL,components_,rJump_[1]*components_,FFTW_ESTIMATE);

//  PrintPlans();
}

#endif

template<class compType>
void PlanFFT<compType>::execute(int fft_type)
{

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

        MPI_Barrier(parallel.lat_world_comm());
        MPI_Alltoall(temp_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, temp1_, (2*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_), MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);
        MPI_Barrier(parallel.dim0_comm()[parallel.grid_rank()[1]]);


        for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*r2cSizeLocal_], r2cSizeLocal_,rSizeLocal_[2],rSizeLocal_[2]);

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


  				   MPI_Alltoall(temp_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, temp1_, 2* rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim1_comm()[parallel.grid_rank()[0]]);

  				   for(i=0;i<parallel.grid_size()[1];i++)transpose_0_2(&temp1_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[1]*rSizeLocal_[2]*rSizeLocal_[1]],rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[1]);

  #ifdef SINGLE
  				   fftwf_execute(fPlan_j_);
  #else
  				   fftw_execute(fPlan_j_);
  #endif


  				   MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);

  				   for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);


  				  for(int l = 0;l< rSizeLocal_[2] ;l++)
            //for(int l = 0;l< 13 ;l++)
  					{

  					    p_in  = &temp_[l*rSizeLocal_[1]];
  					    p_out = &kData_[rJump_[2]*l*components_ + comp];
                //p_out = &kData_[0];
  #ifdef SINGLE
  					    fftwf_execute_dft(fPlan_k_,p_in,p_out);
  #else
  					    fftw_execute_dft(fPlan_k_,p_in,p_out);
  #endif
  					}

            //p_in = NULL;
            //p_out = NULL;
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


  				  //MPI_Alltoall(temp_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1],MPI_DATA_PREC,temp1_,2*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1], MPI_DATA_PREC, parallel.dim0_comm()[parallel.grid_rank()[1]]);

  				   //for(i=0;i<parallel.grid_size()[0];i++)transpose_1_2(&temp1_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]],&temp_[i*rSizeLocal_[2]*rSizeLocal_[2]*rSizeLocal_[1]], rSizeLocal_[1],rSizeLocal_[2],rSizeLocal_[2]);


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
                //fftwf_execute(bPlan_i_);
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


#endif
