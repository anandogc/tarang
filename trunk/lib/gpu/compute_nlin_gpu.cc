
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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 


#include "gpu.h"

//****************************************************************************************					
 
void Compute_nlin_fluid_gpu()
{	
	Tnow_gpu += Tdt_gpu;
	
	// Diag-nlin computation
	Real_space_multiply_gpu(U_V1r_gpu, U_V1r_gpu, Xr_gpu);
	Forward_transform(Xr_gpu, U_nlin1_gpu);
	Xderiv(U_nlin1_gpu, U_nlin1_gpu);
	
	Real_space_multiply_gpu(U_V2r_gpu, U_V2r_gpu, Xr_gpu);
	Forward_transform(Xr_gpu, U_nlin2_gpu);
	Xderiv(U_nlin2_gpu, U_nlin2_gpu);
	
	Real_space_multiply_gpu(U_V3r_gpu, U_V3r_gpu, Xr_gpu);
	Forward_transform(Xr_gpu, U_nlin3_gpu);
	Xderiv(U_nlin3_gpu, U_nlin3_gpu);
	
	// Off-diag nlin computations
	
	Xr_gpu = U_V1r_gpu;
	Real_space_multiply_gpu(U_V1r_gpu, U_V2r_gpu, U_V1r_gpu);			// U_V1r= U_V1rU_V2r
	Real_space_multiply_gpu(U_V2r_gpu, U_V3r_gpu, U_V2r_gpu);			// U_V2r= U_V2rU_V3r
	Real_space_multiply_gpu(U_V3r_gpu, Xr_gpu, U_V3r_gpu);			// U_V3r= U_V3rU_V1r
	
	// Forward transform(U_V1rU_V2r) and derivatives.
	Forward_transform(U_V1r_gpu, X_gpu);
	
	Add_Yderiv(X_gpu, U_nlin1_gpu);
	Add_Xderiv(X_gpu, U_nlin2_gpu);
	
	// Forward transform(U_V2rU_V3r) and derivatives.
	Forward_transform(U_V1r_gpu, X_gpu);
	
	Add_Zderiv(X_gpu, U_nlin2_gpu);
	Add_Yderiv(X_gpu, U_nlin3_gpu);
	
	// Forward transform(U_V1rU_V3r) and derivatives.
	Forward_transform(U_V3r_gpu, X_gpu);
	
	Add_Zderiv(X_gpu, U_nlin1_gpu);
	Add_Xderiv(X_gpu, U_nlin3_gpu);
}


//********************************** End of Ifluid_main.cc ************************************



