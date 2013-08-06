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
#include "Csf.h"


//*********************************************************************************************

//  Class constructor; Allocation for F, Ek etc. Initialize the arrays.

RSF::RSF(string field_name)
{
	this->field_name=field_name;
    
    Fr.resize(shape_real_array);
	
}

/**********************************************************************************************

					Inplace Forward Fourier transform
   
**********************************************************************************************/

/*void RSF::Forward_transform()
{
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Forward_transform_array(Fr);			
}*/



/**********************************************************************************************

					Forward_transform_transpose_order(*Fr) = F 
						*Fr unchanged
   
**********************************************************************************************/

void RSF::Forward_transform(CSF csf)
{
	global.program.sincostr_switch = sincostr_switch_F;
	
	universal->Forward_transform(Fr, csf.F);
}



/**********************************************************************************************

					Inplace Inverse Fourier transform
   
**********************************************************************************************/

/*void RSF::Inverse_transform()
{
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Inverse_transform_array(Fr);			
}*/

/**********************************************************************************************

					Inverse_transform_transpose_order(F) = *Fr 
						F unchanged
   
**********************************************************************************************/

void RSF::Inverse_transform(CSF csf)
{
	global.program.sincostr_switch = sincostr_switch_F;
	
	universal->Inverse_transform(csf.F, Fr);
}






/*******************************************************************************

    Outputs F in real space
PS: In Pencil:SFF & SSF basis, z axis along the mid axis.  So complex data
 is jumbled up.  Therefore, complex<Ar> is copied to real<Ar> for writing.
 For reading, do the reverse.
   
*******************************************************************************/

// field_kind = Tr
void RSF::Write_real_field()
{
 /*   if (global.program.decomposition == "SLAB")
        BasicIO::Write_real_array(Fr, field_name+".Fr");
    
    else if (global.program.decomposition == "PENCIL") {
        if  ((global.program.basis_type == "FFF") || (global.program.basis_type == "SSS"))
            BasicIO::Write_real_array(Fr, field_name+".Fr");
        else {
            universal->Copy_realarray_complex_to_real(Fr, global.temp_array.Ar_copy);
            BasicIO::Write_real_array(global.temp_array.Ar_copy, field_name+".Fr");
        }
    }*/
}


void RSF::Read_real_field()
{
 /*   if (global.program.decomposition == "SLAB")
        BasicIO::Read_real_array(Fr, field_name+".Fr");
    
    else if (global.program.decomposition == "PENCIL") {
       if  ((global.program.basis_type == "FFF") || (global.program.basis_type == "SSS"))
           BasicIO::Read_real_array(Fr, field_name+".Fr");
       else {
           BasicIO::Read_real_array(global.temp_array.Ar_copy, field_name+".Fr");
           universal->Copy_realarray_real_to_complex(global.temp_array.Ar_copy, Fr);
       }
    }*/
}


//************************ END of RSF class Definitions ***************************************



