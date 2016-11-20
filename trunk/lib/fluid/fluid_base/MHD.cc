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

/*! \file elsasser.cc
 * 
 * CONVERTS FIELD VARS AND NLINS FROM U_B TO ELSASSER VARS AND VICE VERSA
 *		
 *	void UB_to_Elsasser_field(IncVF& W) -- U(Zp) = U + W; W(Zm) = U - W;
 *	void Elsasser_to_UB_field(IncVF& W) -- U = (U(Zp) + W(Zm))/2; W =  (U(Zp) - W(Zm))/2; 
 *	
 *	Same as above for force.
 *	
 *	void UB_to_Elsasser_nlin(IncVF& W) 
 *		-- U.nlin(Zp) = U.nlin + W.nlin; W(Zm) = U.nlin - W.nlin;
 *	void Elsasser_to_UB(IncVF& W) 
 *		-- U.nlin = (U(Zp).nlin + W(Zm).nlin)/2; W =  (U.nlin(Zp) - W.nlin(Zm))/2; 
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug   No known bugs
 */

#include "MHD.h"

//*********************************************************************************************

void MHD::UB_to_Elsasser_field(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.cvf.V1;   U.cvf.V1 = U.cvf.V1 + W.cvf.V1;   W.cvf.V1 = global.temp_array.X - W.cvf.V1; 
	global.temp_array.X = U.cvf.V2;   U.cvf.V2 = U.cvf.V2 + W.cvf.V2;   W.cvf.V2 = global.temp_array.X - W.cvf.V2; 
	global.temp_array.X = U.cvf.V3;   U.cvf.V3 = U.cvf.V3 + W.cvf.V3;   W.cvf.V3 = global.temp_array.X - W.cvf.V3; 
}

void MHD::Elsasser_to_UB_field(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.cvf.V1; 
	U.cvf.V1 = Complex(0.5,0)*(U.cvf.V1 + W.cvf.V1);   
	W.cvf.V1 =  Complex(0.5,0)*(global.temp_array.X - W.cvf.V1);
	
	global.temp_array.X = U.cvf.V2; 
	U.cvf.V2 = Complex(0.5,0)*(U.cvf.V2 + W.cvf.V2);   
	W.cvf.V2 =  Complex(0.5,0)*(global.temp_array.X - W.cvf.V2);
	
	global.temp_array.X = U.cvf.V3; 
	U.cvf.V3 = Complex(0.5,0)*(U.cvf.V3 + W.cvf.V3);   
	W.cvf.V3 =  Complex(0.5,0)*(global.temp_array.X - W.cvf.V3);
}


void MHD::UB_to_Elsasser_real_field(FluidVF& U, FluidVF& W)
{

    global.temp_array.Xr = U.rvf.V1r;   
    U.rvf.V1r = U.rvf.V1r + W.rvf.V1r;   W.rvf.V1r = global.temp_array.Xr - W.rvf.V1r; 
    
    global.temp_array.Xr = U.rvf.V2r;   
    U.rvf.V2r = U.rvf.V2r + W.rvf.V2r;   W.rvf.V2r = global.temp_array.Xr - W.rvf.V2r; 
    
    global.temp_array.Xr = U.rvf.V3r;   
    U.rvf.V3r = U.rvf.V3r + W.rvf.V3r;   W.rvf.V3r = global.temp_array.Xr - W.rvf.V3r;
	

}


void MHD::Elsasser_to_UB_real_field(FluidVF& U, FluidVF& W)
{

    global.temp_array.Xr = U.rvf.V1r;   
    U.rvf.V1r = (U.rvf.V1r + W.rvf.V1r)/2;	W.rvf.V1r = (global.temp_array.Xr - W.rvf.V1r)/2;
    
    global.temp_array.Xr = U.rvf.V2r;   
    U.rvf.V2r = (U.rvf.V2r + W.rvf.V2r)/2;	W.rvf.V2r = (global.temp_array.Xr - W.rvf.V2r)/2;
    
    global.temp_array.Xr = U.rvf.V3r;   
    U.rvf.V3r = (U.rvf.V3r + W.rvf.V3r)/2;	W.rvf.V3r = (global.temp_array.Xr - W.rvf.V3r)/2;
}

// Force field
//

void MHD::UB_to_Elsasser_force(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.Force1;   
	U.Force1 = U.Force1 + W.Force1;   W.Force1 = global.temp_array.X - W.Force1;
	
	global.temp_array.X = U.Force2;   
	U.Force2 = U.Force2 + W.Force2;   W.Force2 = global.temp_array.X - W.Force2; 
	
	global.temp_array.X = U.Force3;   
	U.Force3 = U.Force3 + W.Force3;   W.Force3 = global.temp_array.X - W.Force3; 
}

void MHD::Elsasser_to_UB_force(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.Force1;   
	U.Force1 = Complex(0.5,0)*(U.Force1 + W.Force1);   
	W.Force1 = Complex(0.5,0)*(global.temp_array.X - W.Force1);
	
	global.temp_array.X = U.Force2;   
	U.Force2 = Complex(0.5,0)*(U.Force2 + W.Force2);   
	W.Force2 = Complex(0.5,0)*(global.temp_array.X - W.Force2);
	
	global.temp_array.X = U.Force3;   
	U.Force3 = Complex(0.5,0)*(U.Force3 + W.Force3);   
	W.Force3 = Complex(0.5,0)*(global.temp_array.X - W.Force3);
}

//
// nlin
//

void MHD::UB_to_Elsasser_nlin(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.nlin1;   
	U.nlin1 = U.nlin1 + W.nlin1;   W.nlin1 = global.temp_array.X - W.nlin1;
	
	global.temp_array.X = U.nlin2;   
	U.nlin2 = U.nlin2 + W.nlin2;   W.nlin2 = global.temp_array.X - W.nlin2; 
	
	global.temp_array.X = U.nlin3;   
	U.nlin3 = U.nlin3 + W.nlin3;   W.nlin3 = global.temp_array.X - W.nlin3; 
}

void MHD::Elsasser_to_UB_nlin(FluidVF& U, FluidVF& W)
{
	global.temp_array.X = U.nlin1;   
	U.nlin1 = Complex(0.5,0)*(U.nlin1 + W.nlin1);   
	W.nlin1 = Complex(0.5,0)*(global.temp_array.X - W.nlin1);
	
	global.temp_array.X = U.nlin2;   
	U.nlin2 = Complex(0.5,0)*(U.nlin2 + W.nlin2);   
	W.nlin2 = Complex(0.5,0)*(global.temp_array.X - W.nlin2);
	
	global.temp_array.X = U.nlin3;   
	U.nlin3 = Complex(0.5,0)*(U.nlin3 + W.nlin3);   
	W.nlin3 = Complex(0.5,0)*(global.temp_array.X - W.nlin3);
}

//******************************  End of Elasser.cc  ******************************************

	
