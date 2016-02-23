#ifndef LATFIELD2D_PLANFFT_HPP
#define LATFIELD2D_PLANFFT_HPP

/*! \file LATfield2_PlanFFT.hpp
 \brief FFT wrapper
 LATfield2_PlanFFT.hpp contain the class PlanFFT definition.
 */

#include "LATfield2_PlanFFT_decl.hpp"




const int FFT_FORWARD = 1;
const int FFT_BACKWARD = -1;
const int FFT_IN_PLACE = 16;
const int FFT_OUT_OF_PLACE = -16;


#ifndef DOXYGEN_SHOULD_SKIP_THIS


temporaryMemFFT::temporaryMemFFT()
{
    allocated_=0;
#pragma acc enter data create(this)
    
}
temporaryMemFFT::~temporaryMemFFT()
{
  
    freeTemp();
    
    freeTempReal();
  
#pragma acc exit data delete(this)
    
}

temporaryMemFFT::temporaryMemFFT(long size)
{

  setTemp(size);
#pragma acc enter data copyin(this)
#pragma acc enter data create(temp1_[0:alocSize][0:2])
#pragma acc enter data create(temp2_[0:alocSize][0:2])
}

void temporaryMemFFT::freeTemp() {
  if(allocated_!=0){
    
#pragma acc exit data delete(temp1_)
#pragma acc exit data delete(temp2_)
    
#ifndef FFT3D_ACC
#ifdef SINGLE
    fftwf_free(temp1_);
    fftwf_free(temp2_);
#endif
#ifndef SINGLE
    fftw_free(temp1_);
    fftw_free(temp2_);
#endif
#else
    free(temp1_);
    free(temp2_);
#endif
  }
  allocated_ = 0;
}


void temporaryMemFFT::freeTempReal() {
  if(allocatedReal_!=0)
  {
#pragma acc exit data delete(tempReal_)
    
    free(tempReal_);
  }
  allocatedReal_ = 0;
}

int temporaryMemFFT::setTempReal(long size)
{
    if(size>allocatedReal_)
    {
      
        long alocSize;
        alocSize = 2*size;
        
#ifdef SINGLE
        tempReal_ = (float*)malloc(alocSize * sizeof(float));
#else
        tempReal_ = (double*)malloc(alocSize * sizeof(double));
#endif
        
#pragma acc enter data create(tempReal_[0:alocSize])
        
    }
}
int temporaryMemFFT::setTemp(long size)
{
    if(size>allocated_)
    {
        freeTemp();
      
        long alocSize;
        alocSize = 2*size;
        
#ifndef FFT3D_ACC
#ifdef SINGLE
        temp1_ = (fftwf_complex *)fftwf_malloc(alocSize*sizeof(fftwf_complex));
        temp2_ = (fftwf_complex *)fftwf_malloc(alocSize*sizeof(fftwf_complex));
#endif
        
#ifndef SINGLE
        temp1_ = (fftw_complex *)fftw_malloc(alocSize*sizeof(fftw_complex));
        temp2_ = (fftw_complex *)fftw_malloc(alocSize*sizeof(fftw_complex));
#endif
#else
#ifdef SINGLE
        temp1_ = (fftwf_complex *)malloc(alocSize*sizeof(fftwf_complex));
        temp2_ = (fftwf_complex *)malloc(alocSize*sizeof(fftwf_complex));
#endif
        
#ifndef SINGLE
        temp1_ = (fftw_complex *)malloc(alocSize*sizeof(fftw_complex));
        temp2_ = (fftw_complex *)malloc(alocSize*sizeof(fftw_complex));
#endif
#endif
        /// maybe add temp2_ if segfault
#pragma acc enter data create(temp1_[0:alocSize][0:2])
#pragma acc enter data create(temp2_[0:alocSize][0:2])
    }

    allocated_=size;

    return 1;
}

#endif
//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;




#endif
