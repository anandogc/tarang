/* Tarang-2
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-2 .
 *
 * Tarang-2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file	basis_basicfn.h
 * 
 * @brief basic functions declarations that are common to all basis functions.
 *
 *	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *	with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$.
 *
 *	ring(n,m) = \f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *	with typical \f$ R^{sh} = 0,4,8,..., R_{max},\infty \f$ and
 *		\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals.  \f$\Theta_{max} \f$ could be
 *		\f$ \pi \f$ or \f$ \pi/2 \f$ (see function Get_max_anis_theta).
 *
 * @version 4.0 Parallel version
 * @author  M. K. Verma, A. G. Chatterjee
 * @date	Sept 2008
 * @bug		Precision=8 for double. Put condition for float.
 */ 

#ifndef _DEF_VARS_H
#define _DEF_VARS_H

#include <mpi.h>

#include <blitz/array.h>

#include <complex>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#include <sstream>		//Used in the function to_string(const T& t)
#include <sys/stat.h>	//Used for creating data_dir/out (if does not exist) by Global and subdirectories within it by BasicIO.



// RANDOM NOS
#include <random/uniform.h>
#include <random/normal.h>
#include <random/exponential.h>
#include <random/discrete-uniform.h>
#include <random/beta.h>
#include <random/gamma.h>
#include <random/chisquare.h>
#include <random/F.h>

#include <h5file.h>

using namespace ranlib;
// Random ends


using namespace blitz ;


 
//*********************************************************************************************
// DEF vars  


#if defined(REAL_FLOAT)

#define Real 								float
#define MPI_Real							MPI_FLOAT
#define H5Real 								"float"		 // for HDF5
#define H5Complex							"cfloat"     // for HDF5

#define ZERO 0.0f
#define ONE 1.0f
#define TWO 2.0f
#define THREE 3.0f
#define FOUR 4.0f
#define FIVE 5.0f
#define SIX 6.0f
#define SEVEN 7.0f
#define EIGHT 8.0f
#define NINE 9.0f
#define TEN 10.0f
#define ELEVEN 11.0f

const float MY_PI=3.141592653589793238f;

const int MY_PRECISION = 6;

#elif defined(REAL_DOUBLE)

#define Real								double
#define MPI_Real							MPI_DOUBLE
#define H5Real 								"double"		 // for HDF5
#define H5Complex							"cdouble"     // for HDF5

#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0
#define THREE 3.0
#define FOUR 4.0
#define FIVE 5.0
#define SIX 6.0
#define SEVEN 7.0
#define EIGHT 8.0
#define NINE 9.0
#define TEN 10.0
#define ELEVEN 11.0

const double MY_PI=3.141592653589793238;

const int MY_PRECISION = 12;

#endif

#ifndef Complex
#define Complex  complex<Real>
#endif

// Switches for FFTW plans 

#if defined(ESTIMATE)
#define FFTW_PLAN_FLAG FFTW_ESTIMATE

#elif defined(MEASURE)
#define FFTW_PLAN_FLAG FFTW_MEASURE

#elif defined(PATIENT)
#define FFTW_PLAN_FLAG FFTW_PATIENT

#elif defined(EXHAUSTIVE)
#define FFTW_PLAN_FLAG FFTW_EXHAUSTIVE
#endif


//----------------------------------------------------------------------------------------
  
//Convert any data type to string
template <class T>
inline std::string To_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

  
#endif

//***********************************  End of extern_vars.h  ********************************




