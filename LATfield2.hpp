#ifndef LATFIELD2_HPP
#define LATFIELD2_HPP

/*! \file LATfield2.hpp
 \brief LATfield2 header
 \author David Daverio,Neil Bevis

 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <list>
#include <iomanip>

#ifdef FFT3D
#include "fftw3.h"
#endif

#ifndef HDF5
#define HDF5
#endif

#ifdef HDF5
  #include "hdf5.h"



  #ifdef SINGLE
    #define REAL_TYPE H5T_NATIVE_FLOAT
  #else
    #define REAL_TYPE H5T_NATIVE_DOUBLE
  #endif

  #ifdef BIG_ENDIAN_ORDER
    #define DATA_ORDER H5T_ORDER_BE
    #define LONG_TYPE_FILE H5T_STD_I64BE
    #define INT_TYPE_FILE H5T_STD_I32BE
    #define FLOAT_TYPE_FILE H5T_IEEE_F32BE
    #define DOUBLE_TYPE_FILE H5T_IEEE_F64BE
    #define BOOL_TYPE_FILE H5T_STD_I8BE
    #ifdef SINGLE
      #define REAL_TYPE_FILE H5T_IEEE_F32BE
    #else
      #define REAL_TYPE_FILE H5T_IEEE_F64BE
    #endif
  #else
    #define  DATA_ORDER H5T_ORDER_LE
    #define LONG_TYPE_FILE H5T_STD_I64LE
    #define INT_TYPE_FILE H5T_STD_I32LE
    #define FLOAT_TYPE_FILE H5T_IEEE_F32LE
    #define DOUBLE_TYPE_FILE H5T_IEEE_F64LE
    #define BOOL_TYPE_FILE H5T_STD_I8LE
    #ifdef SINGLE
      #define REAL_TYPE_FILE H5T_IEEE_F32LE
    #else
      #define REAL_TYPE_FILE H5T_IEEE_F64LE
    #endif
  #endif
#endif


#if defined(B_SSE2)
  #define ALIGNEMENT 16
  #define NUM_FLOATS  4
  #define NUM_DOUBLES 2
  #define VECT_DOUBLE Vec2d
  #define VECT_FLOAT  Vec4f
#elif defined(B_AVX)
  #define ALIGNEMENT 32
  #define NUM_FLOATS  8
  #define NUM_DOUBLES 4
#elif defined(B_AVX512)
  #define ALIGNEMENT 64
  #define NUM_FLOATS 16
  #define NUM_DOUBLES 8
#elif defined(NO_VECTORIZATION)
  #define ALIGNEMENT 1 /* not used */
  #define NUM_FLOATS 1
#else
	#error “Choose the build target, please.”
#endif

using namespace std;

namespace LATfield2
{

#include "LATfield2_Lattice_decl.hpp"
#include "LATfield2_vector.hpp"
}

#ifdef EXTERNAL_IO
#include "LATfield2_IO_server.hpp"
IOserver ioserver;
#endif

#include "LATfield2_parallel2d.hpp"
Parallel2d parallel;


#include "LATfield2_SettingsFile.hpp"



#ifdef HDF5
#include "LATfield2_save_hdf5.h"
#endif


namespace LATfield2
{
        #include "int2string.hpp"
        #include "Imag.hpp"
        #include "LATfield2_Lattice.hpp"
        #include "LATfield2_Site.hpp"
        #include "LATfield2_Field.hpp"
        #ifdef FFT3D
            #include "LATfield2_PlanFFT.hpp"
        #endif
        #include "particles/LATfield2_Particles.hpp"

//macros
        #include  "looping_macro.hpp"
        #include "timer.hpp"

}



#endif
