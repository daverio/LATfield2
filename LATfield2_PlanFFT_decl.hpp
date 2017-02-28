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


#ifdef SINGLE
#define FFTW_CPLX fftwf_complex
#define FFTW_PLAN fftwf_plan
#else
#define FFTW_CPLX fftw_complex
#define FFTW_PLAN fftw_plan
#endif
/* these are not the actual compiler variables, these are
   only the declaration that they exist somewhere at all. */

extern  const int FFT_FORWARD;
extern  const int FFT_BACKWARD;
extern  const int FFT_IN_PLACE;
extern  const int FFT_OUT_OF_PLACE;

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


  void execute(int fft_type);

private:
  void PrintPlans();

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

  Real * rData_; //pointer to start of data (halo skip)
  FFTW_CPLX * cData_; //pointer to start of data (halo skip)
  FFTW_CPLX * kData_; //pointer to start of data (halo skip)
  FFTW_CPLX * temp_;
  FFTW_CPLX * temp1_;//needed if field got more than 1 component


  FFTW_PLAN fPlan_i_;
  FFTW_PLAN fPlan_j_;
  FFTW_PLAN fPlan_k_;
  FFTW_PLAN fPlan_k_real_;

  FFTW_PLAN bPlan_i_;
  FFTW_PLAN bPlan_j_;
  FFTW_PLAN bPlan_j_real_;
  FFTW_PLAN bPlan_k_;
};

#endif
