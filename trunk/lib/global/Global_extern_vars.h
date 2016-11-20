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

#ifndef _EXTERN_VAR_GLOBAL_H
#define _EXTERN_VAR_GLOBAL_H



#include "Global.h"
extern Global global;

//Aliases to Global.h

//program
extern string basis_type;
extern bool fftw_switch;                 // for fftw_original
extern string sincostr_switch;             // SINX, COSX, SCC, CSC, CCS, SSC, CSS, SCS
extern string sincostr_switch_Vx, sincostr_switch_Vy, sincostr_switch_Vz;
extern string sincostr_switch_F;
extern string sincostr_switch_Visqr;
extern string sincostr_switch_VxVy, sincostr_switch_VxVz, sincostr_switch_VyVz;
extern string sincostr_switch_FVx, sincostr_switch_FVy, sincostr_switch_FVz;
extern string sincostr_switch_divergence;


//field
extern int N[4];
extern vector<Real> kfactor;
extern int Nx;
extern int Ny;
extern int Nz;

//mpi
extern int num_p_cols;
extern int num_p_rows;

extern int num_x_procs;
extern int num_y_procs;
extern int num_z_procs;

extern int num_x_procs_real;
extern int num_y_procs_real;
extern int num_z_procs_real;

extern int my_id;
extern int numprocs;
extern int master_id;
extern bool master;

extern int my_hor_pcoord;
extern int my_vert_pcoord;

extern int my_x_pcoord;
extern int my_y_pcoord;
extern int my_z_pcoord;

extern int my_x_pcoord_real;
extern int my_y_pcoord_real;
extern int my_z_pcoord_real;



//fft
extern int maxlx;
extern int maxly;
extern int maxlz;

extern int lx_start;
extern int ly_start;
extern int lz_start;

extern int maxrx;
extern int maxry;
extern int maxrz;

extern int rx_start;
extern int ry_start;
extern int rz_start;


extern ptrdiff_t local_Nx;
extern ptrdiff_t local_Nx_start;
extern ptrdiff_t local_Nx_hor;
extern ptrdiff_t local_Nx_vert;
extern ptrdiff_t local_Nx_real;



extern ptrdiff_t local_Ny;
extern ptrdiff_t local_Ny_start;
extern ptrdiff_t local_Ny_hor;
extern ptrdiff_t local_Ny_vert;
extern ptrdiff_t local_Ny_real;



extern ptrdiff_t local_Nz;
extern ptrdiff_t local_Nz_start;
extern ptrdiff_t local_Nz_hor;
extern ptrdiff_t local_Nz_vert;
extern ptrdiff_t local_Nz_real;


extern TinyVector<int, 3> shape_complex_array;
extern TinyVector<int, 3> shape_real_array;

/*
extern int complex_arraydim_1;
extern int complex_arraydim_2;
extern int complex_arraydim_3;
extern int real_arraydim_1;
extern int real_arraydim_2;
extern int real_arraydim_3;*/

//Constants
extern Complex  I;

/// minusI = -sqrt(-1).			
extern Complex  minusI;

/// minus2I = -2*sqrt(-1).			
extern Complex  minus2I;

extern Real  MYEPS;
extern Real  MYEPS2;

extern int MY_MAX_INT;
// cut off while reading diagnostic_procedure() array from input file and similar ops

/// Infinite radius.. All the modes outside -- for flux and shelltr calc
extern Real INF_RADIUS; 
extern Real INF_TIME;

#endif

//***********************************  End of extern_vars.h  ********************************




