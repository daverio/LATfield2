#ifndef LATFIELD2D_PLANFFT_HPP
#define LATFIELD2D_PLANFFT_HPP

/*! \file LATfield2_PlanFFT.hpp
 \brief FFT wrapper
 LATfield2_PlanFFT.hpp contain the class PlanFFT definition.
 */


const int FFT_FORWARD = 1;
const int FFT_BACKWARD = -1;
const int FFT_IN_PLACE = 16;
const int FFT_OUT_OF_PLACE = -16;

temporaryMemFFT tempMemory;

#include "LATfield2_PlanFFT_tempMem.hpp"
#include "LATfield2_PlanFFT_transpose.hpp"
#include "LATfield2_PlanFFT_decl.hpp"


#endif
