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
#include "plans/FFF_slab/fff_slab_isend_recv_overlap_isend_forward_3d.h"
#include "plans/FFF_slab/fff_slab_isend_recv_overlap_isend_both_3d.h"
#include "plans/FFF_slab/fff_slab_alltoall_3d.h"

#include "plans/SFF_slab/sff_slab_isend_recv_3d.h"
#include "plans/SFF_slab/sff_slab_alltoall_3d.h"


//*********************************************************************************************	

class SpectralTransform
{
	SpectralPlan *plan;
	
public:
	
	size_t Nx;
	size_t Ny;
	size_t Nz;
	
	size_t local_Nx;
	size_t local_Ny;
	size_t local_Nz;
	
	size_t local_Nx_start;
	size_t local_Ny_start;
	size_t local_Nz_start;

	void Init(string basis, string docomposition, int plan_id, int N0, int N1, int N2, int num_p_hor=1);

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


