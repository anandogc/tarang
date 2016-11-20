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

/*! \file  Cvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "PlainFluidVF.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
   
**********************************************************************************************/

PlainFluidVF::PlainFluidVF
(
	string field_name
) : cvf(field_name), rvf(field_name) 
{}


// Vir = inv_transform(Vi)	
void PlainFluidVF::Inverse_transform()
{
    rvf.Inverse_transform(cvf); 
}
	
	
// Vi = Forward_transform(Vir)		
void PlainFluidVF::Forward_transform()
{
    cvf.Forward_transform(rvf);
}	

// rvf.V1 = Inv FT(cvf.V1)
void PlainFluidVF::Inverse_transform_first_component()
{
    global.program.sincostr_switch = sincostr_switch_Vx;
    universal->Inverse_transform(cvf.V1, rvf.V1r);	
}

// cvf.V1 = FT(rvf.V1)
void PlainFluidVF::Forward_transform_first_component()
{
    global.program.sincostr_switch = sincostr_switch_Vx;
    universal->Forward_transform(rvf.V1r, cvf.V1);
}


//************************ END of CVF class Definitions ***************************************




