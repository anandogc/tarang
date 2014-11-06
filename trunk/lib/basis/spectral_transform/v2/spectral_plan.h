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


#ifndef _SPECTRAL_PLAN_H
#define _SPECTRAL_PLAN_H

#include "spectral_def_vars.h"


//*********************************************************************************************	

class SpectralPlan
{
protected:
	MPI_Request *request;
	MPI_Status *status;

public:
	string basis;

	int my_id;
	int numprocs;

	int my_id_col;
	int my_id_row;

	int num_p_cols;
	int num_p_rows;

	int Nx;
	int Ny;
	int Nz;

	int local_Nx;
	int local_Ny;
	int local_Nz;

	int local_Nx_col;
	int local_Ny_col;
	int local_Nz_col;

	int local_Nx_row;
	int local_Ny_row;
	int local_Nz_row;

	int local_Nx_start;
	int local_Ny_start;
	int local_Nz_start;

	int local_Nx_start_col;
	int local_Ny_start_col;
	int local_Nz_start_col;

	int local_Nx_start_row;
	int local_Ny_start_row;
	int local_Nz_start_row;

	double time_per_step;
	double *timers;
	int num_timers;


	int *PAPI_events;
	long long counters[3];
	long long dummy_counters[3];
	int papi_err;


	SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Nz);
	SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Ny, int Nz);
	SpectralPlan(string basis, int my_id, int numprocs, int Nx, int Ny, int Nz, int num_p_col);

	virtual void Init_array() = 0;
	virtual void Evaluate_time_per_step(double num_iter) = 0;
	virtual void Evaluate_time_per_step(string sincostr_option, double num_iter) = 0;

	//3D FFF
	virtual void Forward_transform(Array<DP,3> Ar, Array<complx,3> A);
	virtual void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar);

	//3D SFF, SSF, SSS
	virtual void Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A);
	virtual void Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar);

	//3D Transpose
	virtual void Transpose(Array<DP,3> Ar, Array<complx,3> A);
	virtual void Transpose(Array<complx,3> A, Array<DP,3> Ar);

	//2D FFF
	virtual void Forward_transform(Array<DP,2> Ar, Array<complx,2> A);
	virtual void Inverse_transform(Array<complx,2> A, Array<DP,2> Ar);

	//2D SFF, SSF, SSS
	virtual void Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A);
	virtual void Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar);
	
	//2D Transpose
	virtual void Transpose(Array<DP,2> Ar, Array<complx,2> A);
	virtual void Transpose(Array<complx,2> A, Array<DP,2> Ar);

};

#endif

//******************************** End of field_basic.h  **************************************



