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


#include "def_vars.h"

//Aliases to Global.h



//program
string basis_type;
bool fftw_switch;           // for fftw_original
string sincostr_switch;     // SINX, COSX, SCC, CSC, CCS, SSC, CSS, SCS
string sincostr_switch_Vx;
string sincostr_switch_Vy;
string sincostr_switch_Vz;
string sincostr_switch_F;
string sincostr_switch_Visqr;
string sincostr_switch_VxVy;
string sincostr_switch_VxVz;
string sincostr_switch_VyVz;
string sincostr_switch_FVx;
string sincostr_switch_FVy;
string sincostr_switch_FVz;
string sincostr_switch_divergence;

//field
int N[4];
vector<Real> kfactor;

//mpi
int num_p_cols;
int num_p_rows;
int num_x_procs;
int num_y_procs;
int num_z_procs;

int num_x_procs_real;
int num_y_procs_real;
int num_z_procs_real;

int my_id;
int numprocs;
int master_id;
bool master;


//fft

int maxlx;
int maxly;
int maxlz;

int lx_start;
int ly_start;
int lz_start;

int maxrx;
int maxry;
int maxrz;

int rx_start;
int ry_start;
int rz_start;


int Nx;
ptrdiff_t local_Nx;
ptrdiff_t local_Nx_start;
ptrdiff_t local_Nx_hor;
ptrdiff_t local_Nx_vert;
ptrdiff_t local_Nx_real;


int Ny;
ptrdiff_t local_Ny;
ptrdiff_t local_Ny_start;
ptrdiff_t local_Ny_hor;
ptrdiff_t local_Ny_vert;
ptrdiff_t local_Ny_real;


int Nz;
ptrdiff_t local_Nz;
ptrdiff_t local_Nz_start;
ptrdiff_t local_Nz_hor;
ptrdiff_t local_Nz_vert;
ptrdiff_t local_Nz_real;


int my_hor_pcoord;
int my_vert_pcoord;
int my_x_pcoord;
int my_y_pcoord;
int my_z_pcoord;

int my_x_pcoord_real;
int my_y_pcoord_real;
int my_z_pcoord_real;

TinyVector<int, 3> shape_complex_array;
TinyVector<int, 3> shape_real_array;



/*
int complex_arraydim_1;
int complex_arraydim_2;
int complex_arraydim_3;
int real_arraydim_1;
int real_arraydim_2;
int real_arraydim_3;*/

//Constants
Complex I;		

/// minusI = -sqrt(-1).			
Complex  minusI;

/// minus2I = -2*sqrt(-1).			
Complex  minus2I;

Real  MYEPS;
Real  MYEPS2;

int MY_MAX_INT;
// cut off while reading diagnostic_procedure() array from input file and similar ops

/// Infinite radius.. All the modes outside -- for flux and shelltr calc
Real INF_RADIUS; 
Real INF_TIME;

