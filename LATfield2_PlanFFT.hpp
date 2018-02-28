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
	temp1_=NULL;
	temp2_=NULL;
	allocated_=0;
}
temporaryMemFFT::~temporaryMemFFT()
{

	if(allocated_>0)
	{
#ifdef SINGLE
		//if(temp1_!=NULL)fftwf_free(temp1_);
		//if(temp2_!=NULL)fftwf_free(temp2_);
#endif
#ifndef SINGLE
		//if(temp1_!=NULL)fftw_free(temp1_);
		//if(temp2_!=NULL)fftw_free(temp2_);
#endif
		//allocated_ = 0;
		//temp1_=NULL;
		//temp2_=NULL;
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
		if(allocated_)
		{
			if(temp1_!=NULL)fftwf_free(temp1_);
			if(temp2_!=NULL)fftwf_free(temp2_);
		}
		temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));

		//debug cout<<"("<< parallel.grid_rank()[0]<< ","<<parallel.grid_rank()[1] <<"): called temporary resize. Old size: " <<allocated_<<" , new size: "<< size<<endl;

#endif

#ifndef SINGLE
		if(allocated_)
		{
			if(temp1_!=NULL)fftw_free(temp1_);
			if(temp2_!=NULL)fftw_free(temp2_);
		}
		temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif
		allocated_ = size ;


	}
	return 1;
}


#endif
//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;




#endif
