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



/*! \file  Rvf.cc
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @sa	Rvf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 

#include "Rvf.h"

/*
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
 */

//*********************************************************************************************

RVF::RVF() 
{ 
  	  
	V1r.resize(local_N2, Nrv[1], Nrv[3]/2+1);  
	V2r.resize(local_N2, Nrv[1], Nrv[3]/2+1);  
	V3r.resize(local_N2, Nrv[1], Nrv[3]/2+1);
	
		// If Transpose switch on, FT(*Vtr) = *V
	
#ifdef TRANSPOSE
	V1r = new Array<complx,3>(local_N2, N[1], N[3]/2+1);  
	V2r = new Array<complx,3>(local_N2, N[1], N[3]/2+1);  
	V3r = new Array<complx,3>(local_N2, N[1], N[3]/2+1);  
	
#else 
	V1r = new Array<complx,3>(local_N1, N[2], N[3]/2+1);			
	V2r = new Array<complx,3>(local_N1, N[2], N[3]/2+1);  
	V3r = new Array<complx,3>(local_N1, N[2], N[3]/2+1); 
#endif	

	
	// Structure function
	
	if ( global.structure_function_switch || global.structure_function.planar_switch )
	{
		RV_structurefn_r_max = Max_radius_inside_real_space()/global.structure_fn.rmax_div_st_fn_r_max;
		
		int q_range = global.structure_fn.qmax - global.structure_fn.qmin;
		
		RV_St = new Array<DP,3>(global.structure_fn.array_size+1, q_range+1, 2);
		
		RV_St_count = new Array<DP,1>(global.structure_fn.array_size+1);
		
		// Contains averaged over all procs in master proc.
		RV_St_final = new Array<DP,3>(global.structure_fn.array_size+1, q_range+1, 2);
		
		RV_St_count_final = new Array<DP,1>(global.structure_fn.array_size+1);
		

		// To optimize, we take max distance to be N[i]/RV_rmax_div_st_fn_r_max-1,
	/*	if (RV_structure_fn_approx_switch == 0)
			dist_farthest_proc = numprocs-1;
		
		else if (RV_structure_fn_approx_switch == 1) */
		//	dist_farthest_proc = global.mpi.numprocs/RV_rmax_div_st_fn_r_max-1; 
	}	
	
/*	if (global.structure_function.planar_switch) {
	
		if (global.field.anisotropy_dirn == 1)
			RV_structurefn_pll_ind_max = N[1];
		
		else if (global.field.anisotropy_dirn == 2)
			RV_structurefn_pll_ind_max = N[2];
		
		else if (global.field.anisotropy_dirn == 3)
			RV_structurefn_pll_ind_max = N[3];
	} */
	
	
	
}

/**********************************************************************************************

				Forward transform: Inplace 
   
**********************************************************************************************/

void RVF::RV_Forward_transform()
{
	Forward_transform_array(*V1r, 1);			// SFT
	
	Forward_transform_array(*V2r, 0);			// CFT
	
	Forward_transform_array(*V3r, 0);			// CFT	
}

/**********************************************************************************************

			Forward_transform(*Vir) = Vi 
				Vir unchanged.
   
**********************************************************************************************/


void RVF::RV_Forward_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3 
)
{
	global.temp_array.Ar = *V1r;
	Forward_transform_array_transpose_order(global.temp_array.Ar, V1, 1);			// SFT
	
	global.temp_array.Ar = *V2r;
	Forward_transform_array_transpose_order(global.temp_array.Ar, V2, 0);			// CFT
	
	global.temp_array.Ar = *V3r;
	Forward_transform_array_transpose_order(global.temp_array.Ar, V3, 0);			// CFT	
}


/**********************************************************************************************

		Inplace Inverse Fourier transform.
   
**********************************************************************************************/

void RVF::RV_Inverse_transform()
{

	Inverse_transform_array(V1r, 1);			// ISFT

	Inverse_transform_array(V2r, 0);			// ICFT

	Inverse_transform_array(V3r, 0);			// ICFT

}


/**********************************************************************************************

			Inverse_transform(Vi) = *Vir 
				Keeping Vi unchanged....
				temp = N1 x N2 x N3
   
**********************************************************************************************/

void RVF::RV_Inverse_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3
)
{
	global.temp_array.X = V1;
	Inverse_transform_array_transpose_order(global.temp_array.X, *V1r, 1);	// ISFT
	
	global.temp_array.X = V2;
	Inverse_transform_array_transpose_order(global.temp_array.X, *V2r, 0);	// ICFT
	
	global.temp_array.X = V3;
	Inverse_transform_array_transpose_order(global.temp_array.X, *V3r, 0);	// ICFT

}
	
/**********************************************************************************************
/**********************************************************************************************



/**********************************************************************************************

			 Outputs Vi in real space
   
**********************************************************************************************/

void RVF::RV_Write_real_field()
{
	io.Write_array(V1r);
	io.Write_array(V2r);
	io.Write_array(V3r);
}


void RVF::RV_Read_real_field()
{
	io.Read_array(V1r);
	io.Read_array(V2r);
	io.Read_array(V3r);
}


void RVF::RV_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3)
{
	
	if (global.structure_function.approx_switch) {
		Qi1_start = Pi1;
		Qi2_start = Pi2;
		Qi3_start = Pi3;
		
		Qi1_end = min(((int) (Pi1 + N[1]/global.structure_function.rmax_div_st_fn_r_max)),local_N1-1);
		Qi2_end = min(((int) (Pi2 + N[2]/global.structure_function.rmax_div_st_fn_r_max)), N[2]-1);
		Qi3_end = min(((int) (Pi3 + N[3]/(2*global.structure_function.rmax_div_st_fn_r_max))), N[3]/2-1);
		// To optimize, we take max distance to be min of the above two.
	}
	
	else{
		Qi1_start = Qi2_start = Qi3_start = 0;
		
		Qi1_end = local_N1-1;
		Qi2_end = Nrv[2]-1;
		Qi3_end = Nrv[3]/2-1;		
	}
}

//**************************  End of RVF class definitions ************************************






