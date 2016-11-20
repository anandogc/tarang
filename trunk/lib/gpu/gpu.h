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

/*! \file Ifluid_main.h
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 2 October 2008
 * 
 * @bug  No known bug
 */ 

//********************************************************************************************* 

extern Array<Complex,3> V1_gpu(shape_complex_array);
extern Array<Complex,3> V2_gpu(shape_complex_array);
extern Array<Complex,3> V3_gpu(shape_complex_array);

extern Array<Complex,3> nlin1_gpu(shape_complex_array);
extern Array<Complex,3> nlin2_gpu(shape_complex_array);
extern Array<Complex,3> nlin3_gpu(shape_complex_array);

extern Array<Complex,3> Pressure_gup(shape_complex_array);
Array<Complex,3>  X_gpu(shape_complex_array);


extern Array<Real,3> V1r_gpu(shape_real_array);
extern Array<Real,3> V2r_gpu(shape_real_array);
extern Array<Real,3> V3r_gpu(shape_real_array);

extern Array<Real,3>  Xr_gpu(shape_real_array);


//******************************** End of Ifluid_main.h ***************************************

																																																																																																							
