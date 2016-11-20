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


/*! \file  Rsf.cc
 * 
 * @brief  Class declaration of Rsf, a Real Scalar Field 
 *
 * @sa	Rsf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug No known bugs
 */
 

#include "Rsf.h"


//*********************************************************************************************

//  Class constructor; Allocation for F, Ek etc. Initialize the arrays.

RSF::RSF()
{
	
	// Structure function
	
	
	if ( global.structure_fn.box_switch || global.structure_fn.planar_switch )
	{
		RS_structurefn_r_max = Max_radius_inside_real_space()
								/global.structure_fn.rmax_div_st_fn_r_max;
		
		int q_range = global.structure_fn.qmax - global.structure_fn.qmin;
		
		RS_St = new Array<Real,2>(global.structure_fn.array_size+1, q_range+1);
		
		RS_St_count = new Array<Real,1>(global.structure_fn.array_size+1);
		
		// Contains averaged over all procs in master proc.
		RS_St_final = new Array<Real,2>(global.structure_fn.array_size+1, q_range+1);
		
		RS_St_count_final = new Array<Real,1>(global.structure_fn.array_size+1);
		
		// To optimize, we take max distance to be N[i]/RS_rmax_div_st_fn_r_max,
		/*	if (RS_structure_fn_approx_switch == 0)
			dist_farthest_proc = numprocs-1;
		
		else if (RS_structure_fn_approx_switch == 1) */
		//dist_farthest_proc = numprocs/RS_rmax_div_st_fn_r_max-1;
	}
	
	
	/*if (RS_planar_structure_fn_switch == 1)
	{
		
		if (globalvar_anisotropy_switch == 1)
			RS_structure_fn_rpll_max = Nrs[1];
		
		else if (globalvar_anisotropy_switch == 2)
			RS_structure_fn_rpll_max = Nrs[2];
		
		else if (globalvar_anisotropy_switch == 3)
			RS_structure_fn_rpll_max = Nrs[3];
	}*/
	

#ifdef TRANSPOSE
	Fr = new Array<Complex,3>(local_N2, Nrs[1], Nrs[3]/2+1); 
#else    
    Fr = new Array<Complex,3>(local_N1, Nrs[2], Nrs[3]/2+1);  
#endif	

    *Fr = 0.0;   // initialize
	
}

/**********************************************************************************************

					Inplace Forward Fourier transform
   
**********************************************************************************************/

void RSF::RS_Forward_transform()
{
	Forward_transform_array(*Fr,1);			// SFT
}



/**********************************************************************************************

					Forward_transform_transpose_order(*Fr) = F 
						*Fr unchanged
   
**********************************************************************************************/

void RSF::RS_Forward_transform_transpose_order(Array<Complex,3> F)
{
	global.temp_array.Xr = *Fr;
	
	Forward_transform_array_transpose_order(global.temp_array.Xr, F, 1);		// SFT
}



/**********************************************************************************************

					Inplace Inverse Fourier transform
   
**********************************************************************************************/

void RSF::RS_Inverse_transform()
{
	Inverse_transform_array(*Fr, 1);			// ISFT
}

/**********************************************************************************************

					Inverse_transform_transpose_order(F) = *Fr 
						F unchanged
   
**********************************************************************************************/

void RSF::RS_Inverse_transform_transpose_order(Array<Complex,3> F)
{
	global.temp_array.X = F;
	
	Inverse_transform_array_transpose_order(global.temp_array.X, *Fr, 1);		// ISFT
}


/**********************************************************************************************



/**********************************************************************************************

    Outputs F in real space
   
**********************************************************************************************/

void RSF::RS_Write_real_field()
{
	io.Write_array(Fr);
}


void RSF::RS_Read_real_field()
{
	io.Read_array(Fr);
}


void RSF::RS_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3)
{
	// To optimize, we take max distance to be min(N/5,
	if (RS_structure_fn_approx_switch == 0) {
		Qi1_start = Qi2_start = Qi3_start = 0;
		
		Qi1_end = local_N1-1;
		Qi2_end = Nrs[2]-1;
		Qi3_end = Nrs[3]/2-1;
	}
	
	else if (RS_structure_fn_approx_switch == 1) {
		Qi1_start = Pi1;
		Qi2_start = Pi2;
		Qi3_start = Pi3;
		
		Qi1_end = min(((int) (Pi1 + Nrs[1]/global.structure_function.rmax_div_st_fn_r_max)), local_N1-1);
		Qi2_end = min(((int) (Pi2 + Nrs[2]/global.structure_function.rmax_div_st_fn_r_max)), Nrs[2]-1);
		Qi3_end = min(((int) (Pi3 + Nrs[3]/(2*global.structure_function.rmax_div_st_fn_r_max))), Nrs[3]/2-1);
		// To optimize, we take max distance to be min of the above two.
	}
}

//************************ END of RSF class Definitions ***************************************



