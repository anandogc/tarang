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

#ifndef _SPECTRAL_DEF_VARS_H
#define _SPECTRAL_DEF_VARS_H

#include <mpi.h>
#include <fftw3-mpi.h>

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

using namespace ranlib;
// Random ends


using namespace blitz ;


 
//*********************************************************************************************
// DEF vars  

// Defining FLOAT_DP: switch for setting float

#if defined(FLOAT_DP)

#define DP 								float
#define FFTW_PLAN_DP					fftwf_plan
#define FFTW_MPI_PLAN_DFT_R2C_3D_DP  	fftwf_mpi_plan_dft_r2c_3d
#define FFTW_MPI_PLAN_DFT_C2R_3D_DP  	fftwf_mpi_plan_dft_c2r_3d
#define FFTW_MPI_PLAN_DFT_R2C_2D_DP  	fftwf_mpi_plan_dft_r2c_2d
#define FFTW_MPI_PLAN_DFT_C2R_2D_DP  	fftwf_mpi_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_2D_DP			fftwf_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D_DP			fftwf_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_1D_DP				fftwf_plan_dft_1d
#define FFTW_PLAN_DFT_R2C_1D_DP			fftwf_plan_dft_r2c_1d
#define FFTW_PLAN_DFT_C2R_1D_DP			fftwf_plan_dft_c2r_1d
#define FFTW_PLAN_R2R_1D_DP				fftwf_plan_r2r_1d
#define FFTW_EXECUTE_DFT_R2C_DP			fftwf_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R_DP			fftwf_execute_dft_c2r
#define FFTW_EXECUTE_DFT_DP				fftwf_execute_dft
#define FFTW_EXECUTE_R2R_DP				fftwf_execute_r2r
//#define FFTW_MPI_INIT_DP				fftwf_mpi_init
#define FFTW_COMPLEX_DP					fftwf_complex
#define FFTW_DESTROY_PLAN_DP			fftwf_destroy_plan
#define MPI_DP							MPI_FLOAT
#define FFTW_MPI_LOCAL_SIZE_3D_TRANSPOSED_DP fftwf_mpi_local_size_3d_transposed
#define FFTW_MPI_LOCAL_SIZE_2D_TRANSPOSED_DP fftwf_mpi_local_size_2d_transposed

#define FFTW_MPI_PLAN_DFT_3D_DP         fftwf_mpi_plan_dft_3d   // for GP
#define FFTW_MPI_EXECUTE_DFT_DP         fftwf_mpi_execute_dft

// Defining DOUBLE_DP: switch for setting double

#elif defined(DOUBLE_DP)

#define DP 						double
#define FFTW_PLAN_DP					fftw_plan
#define FFTW_MPI_PLAN_DFT_R2C_3D_DP  	fftw_mpi_plan_dft_r2c_3d
#define FFTW_MPI_PLAN_DFT_C2R_3D_DP  	fftw_mpi_plan_dft_c2r_3d
#define FFTW_MPI_PLAN_DFT_R2C_2D_DP  	fftw_mpi_plan_dft_r2c_2d
#define FFTW_MPI_PLAN_DFT_C2R_2D_DP  	fftw_mpi_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_2D_DP			fftw_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D_DP			fftw_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_1D_DP				fftw_plan_dft_1d
#define FFTW_PLAN_DFT_R2C_1D_DP			fftw_plan_dft_r2c_1d
#define FFTW_PLAN_DFT_C2R_1D_DP			fftw_plan_dft_c2r_1d
#define FFTW_PLAN_R2R_1D_DP				fftw_plan_r2r_1d
#define FFTW_EXECUTE_DFT_R2C_DP			fftw_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R_DP			fftw_execute_dft_c2r
#define FFTW_EXECUTE_DFT_DP				fftw_execute_dft
#define FFTW_EXECUTE_R2R_DP				fftw_execute_r2r
//#define FFTW_MPI_INIT_DP				fftw_mpi_init
#define FFTW_COMPLEX_DP					fftw_complex
#define FFTW_DESTROY_PLAN_DP			fftw_destroy_plan
#define MPI_DP							MPI_DOUBLE
#define FFTW_MPI_LOCAL_SIZE_3D_TRANSPOSED_DP fftw_mpi_local_size_3d_transposed
#define FFTW_MPI_LOCAL_SIZE_2D_TRANSPOSED_DP fftw_mpi_local_size_2d_transposed

#define FFTW_MPI_PLAN_DFT_3D_DP         fftw_mpi_plan_dft_3d  // for GP
#define FFTW_MPI_EXECUTE_DFT_DP         fftw_mpi_execute_dft


#endif

#ifndef complx
#define complx  complex<DP>
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


template<class T, int N_rank>
bool operator==(TinyVector<T,N_rank> A, TinyVector<T,N_rank> B){
	for (int i=0; i<N_rank; i++)
		if (A(i)!=B(i))
			return false;
	return true;
}

    
#endif

//***********************************  End of extern_vars.h  ********************************




