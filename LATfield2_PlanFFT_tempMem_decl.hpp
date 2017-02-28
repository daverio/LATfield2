#ifndef LATFIELD2_PLANFFT_TEMPMEM_DECL_HPP
#define LATFIELD2_PLANFFT_TEMPMEM_DECL_HPP


class temporaryMemFFT
	{
	public:
	    temporaryMemFFT();
	    ~temporaryMemFFT();
	    temporaryMemFFT(long size);

	    int setTemp(long size);
	    int setTempReal(long size);
#ifdef SINGLE
	    fftwf_complex* temp1(){return temp1_;}
	    fftwf_complex* temp2(){return temp2_;}
	    float * tempReal(){return tempReal_;}
#endif

#ifndef SINGLE
	    fftw_complex * temp1(){return temp1_;}
	    fftw_complex * temp2(){return temp2_;}
	    double * tempReal(){return tempReal_;}
#endif

	private:
#ifdef SINGLE
	    fftwf_complex * temp1_;
	    fftwf_complex * temp2_;
	    float * tempReal_;
#endif

#ifndef SINGLE
	    fftw_complex * temp1_;
	    fftw_complex * temp2_;
	    double * tempReal_;
#endif
	    long allocated_; //number of variable stored (bit = allocated*sizeof(fftw(f)_complex))
	    long allocatedReal_;
};


#endif
