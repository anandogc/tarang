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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _SPECTRAL_TRANSFORM_H
#define _SPECTRAL_TRANSFORM_H

#include "spectral_def_vars.h"
#include "spectral_plan.h"

//FFFW
#include "plans/FFFW_slab/2d/fffw_slab_in_order_2d.h"
#include "plans/FFFW_slab/3d/fffw_slab_transposed_order_3d.h"

//FFF:SLAB
#include "plans/FFF_slab/2d/fff_slab_alltoall_2d.h"

#include "plans/FFF_slab/3d/fff_slab_isend_recv_overlap_isend_forward_3d.h"
#include "plans/FFF_slab/3d/fff_slab_isend_recv_overlap_isend_both_3d.h"
#include "plans/FFF_slab/3d/fff_slab_alltoall_3d.h"

//FFF:PENCIL
//#include "plans/FFF_pencil/fff_pencil_alltoall_3d.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_plane_by_plane.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_all_at_once.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_buffered.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_threaded.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_morton.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_compressed_gzip.h"
#include "plans/FFF_pencil/fff_pencil_alltoall_3d_compressed_lz4.h"

//SFF:SLAB
#include "plans/SFF_slab/2d/sf_slab_alltoall_2d.h"

#include "plans/SFF_slab/3d/sff_slab_isend_recv_3d.h"
#include "plans/SFF_slab/3d/sff_slab_alltoall_3d.h"

//SFF:PENCIL
#include "plans/SFF_pencil/sff_pencil_alltoall_3d.h"

//SSF
#include "plans/SSF_slab/3d/ssf_slab_alltoall_3d.h"
#include "plans/SSF_pencil/ssf_pencil_alltoall_3d.h"


//SSS
#include "plans/SSS_slab/2d/ss_slab_alltoall_2d.h"
#include "plans/SSS_slab/3d/sss_slab_alltoall_3d.h"
#include "plans/SSS_pencil/sss_pencil_alltoall_3d.h"


//*********************************************************************************************	

class SpectralTransform
{
public:
	SpectralPlan *plan;
	
	string basis;

	int select_plan_id;

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


	SpectralTransform();

	void Init(string basis, int Nx, int Nz);
	void Init(string basis, int Nx, int Ny, int Nz);
	void Init(string basis, int Nx, int Ny, int Nz, int num_p_hor);

	void Forward_transform(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar);

	void Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar);

	void Transpose(Array<DP,3> Ar, Array<complx,3> A);
	void Transpose(Array<complx,3> A, Array<DP,3> Ar);

	void Forward_transform(Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform(Array<complx,2> A, Array<DP,2> Ar);

	void Forward_transform(string sincostr_option, Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform(string sincostr_option, Array<complx,2> A, Array<DP,2> Ar);

	void Transpose(Array<DP,2> Ar, Array<complx,2> A);
	void Transpose(Array<complx,2> A, Array<DP,2> Ar);

};

#endif


//******************************** End of field_basic.h  **************************************


