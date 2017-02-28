#ifndef LATFIELD2_PLANFFT_TEMPMEM_HPP
#define LATFIELD2_PLANFFT_TEMPMEM_HPP

#include "LATfield2_PlanFFT_tempMem_decl.hpp"


temporaryMemFFT::temporaryMemFFT()
{
    allocated_=0;

}
temporaryMemFFT::~temporaryMemFFT()
{
    freeTemp();
    freeTempReal();
}

temporaryMemFFT::temporaryMemFFT(long size)
{
  setTemp(size);
}

void temporaryMemFFT::freeTemp() {
  if(allocated_!=0){


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
    free(tempReal_);
  }
  allocatedReal_ = 0;
}

int temporaryMemFFT::setTempReal(long size)
{
    if(size>allocatedReal_)
    {
        long alocSize;
        alocSize = size;

#ifdef SINGLE
        tempReal_ = (float*)malloc(alocSize * sizeof(float));
#else
        tempReal_ = (double*)malloc(alocSize * sizeof(double));
#endif
    }
}
int temporaryMemFFT::setTemp(long size)
{
    if(size>allocated_)
    {
        freeTemp();
        
        long alocSize;
        alocSize = size;

#ifdef SINGLE
        temp1_ = (fftwf_complex *)fftwf_malloc(alocSize*sizeof(fftwf_complex));
        temp2_ = (fftwf_complex *)fftwf_malloc(alocSize*sizeof(fftwf_complex));
#endif

#ifndef SINGLE
        temp1_ = (fftw_complex *)fftw_malloc(alocSize*sizeof(fftw_complex));
        temp2_ = (fftw_complex *)fftw_malloc(alocSize*sizeof(fftw_complex));
#endif
    }

    allocated_=size;

    return 1;
}



#endif
