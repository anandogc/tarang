
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



// GPU arrays
Array<Complex,3> V1_gpu(shape_complex_array);
Array<Complex,3> V2_gpu(shape_complex_array);
Array<Complex,3> V3_gpu(shape_complex_array);

Array<Complex,3> nlin1_gpu(shape_complex_array);
Array<Complex,3> nlin2_gpu(shape_complex_array);
Array<Complex,3> nlin3_gpu(shape_complex_array);


Array<Complex,3> Pressure_gup(shape_complex_array);

Array<Complex,3> V1r_gpu(shape_real_array);
Array<Complex,3> V2r_gpu(shape_real_array);
Array<Complex,3> V3r_gpu(shape_real_array);



//****************************************************************************************
 
void Time_iterate_fluid_gpu()
{
	V1_gpu = U.cvf.V1;
	V2_gpu = U.cvf.V2;
	V3_gpu = U.cvf.V3;
	
allocate	V1r_gpu, V2r_gpu,V3r_gpu, nlin1_gpu, 
	
//	fluidIO_incompress.Output_all_inloop(U, P);  // for initial cond

	Tnow_gpu = global.time.init;
	Tfinal_gpu = global.time.final;
	Tdt_gpu = global.time.dt;
	
	do 	{
		if (integration_scheme == "EULER")
			Euler_fluid_gpu();
		
		Compute_total_energy();
		print total_energy;
			
	}
	while	(Tnow_gpu < Tfinal_gpu);
	
	U.cvf.V1 = V1_gpu;
	U.cvf.V2 = V2_gpu;
	U.cvf.V3 = V3_gpu;

} 


//********************************** End of Ifluid_main.cc ************************************	



