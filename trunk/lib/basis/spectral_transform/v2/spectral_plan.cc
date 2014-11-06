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


/*! \file  spectral_plan.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  A. G. Chatterjee, M. K. Verma
 * @version 2
 * @date March 2014
 * @bug	No known bugs
 */

#include "spectral_plan.h"

SpectralPlan::SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Nz)
   :basis(basis),
    my_id(my_id), 
	numprocs(numprocs),
	Nx(Nx),
	Ny(1),
	Nz(Nz),
	time_per_step(0)
{}

SpectralPlan::SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Ny, int Nz)
   :basis(basis),
    my_id(my_id), 
	numprocs(numprocs),
	Nx(Nx),
	Ny(Ny),
	Nz(Nz),
	time_per_step(0)
{}

SpectralPlan::SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Ny, int Nz, int num_p_cols)
   :basis(basis),
    my_id(my_id), 
	numprocs(numprocs),
	Nx(Nx),
	Ny(Ny),
	Nz(Nz),
	num_p_cols(num_p_cols),
	time_per_step(0)
{}


//3D FFF
void SpectralPlan::Forward_transform(Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Inverse_transform(Array<complx,3> A, Array<DP,3> Ar){}

//3D SFF, SSF, SSS
void SpectralPlan::Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar){}

//3D Transpose
void SpectralPlan::Transpose(Array<DP,3> Ar, Array<complx,3> A){}
void SpectralPlan::Transpose(Array<complx,3> A, Array<DP,3> Ar){}

//2D FFF
void SpectralPlan::Forward_transform(Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Inverse_transform(Array<complx,2> A, Array<DP,2> Ar){}

//2D SFF, SSF, SSS
void SpectralPlan::Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar){}

//2D Transpose
void SpectralPlan::Transpose(Array<DP,2> Ar, Array<complx,2> A){}
void SpectralPlan::Transpose(Array<complx,2> A, Array<DP,2> Ar){}


//******************************** End of field_basic.h  **************************************



