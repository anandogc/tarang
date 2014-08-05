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
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;
using namespace blitz;


//*********************************************************************************************

// Defining FLOAT_DP: switch for setting float
#if defined(FLOAT_DP)

#define DP 								float
#define MPI_DP							MPI_FLOAT
#define FFTW_COMPLEX					fftwf_complex

#define FFTW_MPI_INIT                   fftwf_mpi_init

#define FFTW_PLAN						fftwf_plan
#define FFTW_MPI_PLAN_DFT_3D			fftwf_mpi_plan_dft_3d   // for GP
#define FFTW_MPI_PLAN_DFT_R2C_2D		fftwf_mpi_plan_dft_r2c_2d
#define FFTW_MPI_PLAN_DFT_C2R_2D		fftwf_mpi_plan_dft_c2r_2d
#define FFTW_MPI_PLAN_DFT_R2C_3D		fftwf_mpi_plan_dft_r2c_3d
#define FFTW_MPI_PLAN_DFT_C2R_3D		fftwf_mpi_plan_dft_c2r_3d
#define FFTW_PLAN_MANY_R2R				fftwf_plan_many_r2r
#define FFTW_PLAN_MANY_DFT				fftwf_plan_many_dft
#define FFTW_PLAN_MANY_DFT_R2C			fftwf_plan_many_dft_r2c
#define FFTW_PLAN_MANY_DFT_C2R			fftwf_plan_many_dft_c2r

#define FFTW_EXECUTE					fftwf_execute
#define FFTW_EXECUTE_R2R				fftwf_execute_r2r
#define FFTW_EXECUTE_DFT				fftwf_execute_dft
#define FFTW_EXECUTE_DFT_R2C			fftwf_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R			fftwf_execute_dft_c2r
#define FFTW_MPI_EXECUTE_DFT			fftwf_mpi_execute_dft
#define FFTW_MPI_EXECUTE_DFT_R2C        fftwf_mpi_execute_dft_r2c        
#define FFTW_MPI_EXECUTE_DFT_C2R        fftwf_mpi_execute_dft_c2r        


// Defining DOUBLE_DP: switch for setting double
#elif defined(DOUBLE_DP)

#define DP 								double
#define MPI_DP							MPI_DOUBLE
#define FFTW_COMPLEX					fftw_complex

#define FFTW_MPI_INIT                   fftw_mpi_init

#define FFTW_PLAN						fftw_plan
#define FFTW_MPI_PLAN_DFT_3D			fftw_mpi_plan_dft_3d   // for GP
#define FFTW_MPI_PLAN_DFT_R2C_2D		fftw_mpi_plan_dft_r2c_2d
#define FFTW_MPI_PLAN_DFT_C2R_2D		fftw_mpi_plan_dft_c2r_2d
#define FFTW_MPI_PLAN_DFT_R2C_3D		fftw_mpi_plan_dft_r2c_3d
#define FFTW_MPI_PLAN_DFT_C2R_3D		fftw_mpi_plan_dft_c2r_3d
#define FFTW_PLAN_MANY_R2R				fftw_plan_many_r2r
#define FFTW_PLAN_MANY_DFT				fftw_plan_many_dft
#define FFTW_PLAN_MANY_DFT_R2C			fftw_plan_many_dft_r2c
#define FFTW_PLAN_MANY_DFT_C2R			fftw_plan_many_dft_c2r

#define FFTW_EXECUTE					fftw_execute
#define FFTW_EXECUTE_R2R				fftw_execute_r2r
#define FFTW_EXECUTE_DFT				fftw_execute_dft
#define FFTW_EXECUTE_DFT_R2C			fftw_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R			fftw_execute_dft_c2r
#define FFTW_MPI_EXECUTE_DFT			fftw_mpi_execute_dft
#define FFTW_MPI_EXECUTE_DFT_R2C        fftw_mpi_execute_dft_r2c        
#define FFTW_MPI_EXECUTE_DFT_C2R        fftw_mpi_execute_dft_c2r  

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


#define IN_ALL_PROCS(l, expr) \
for (int l=0; l<numprocs; l++) { \
	if (my_id==l)\
		{expr;}\
	MPI_Barrier(MPI_COMM_WORLD);\
}

    
#endif

//***********************************  End of extern_vars.h  ********************************




