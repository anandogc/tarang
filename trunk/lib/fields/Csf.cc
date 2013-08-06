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



/*! \file  Csf.cc
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 * @sa Csf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "Csf.h"
#include "Rsf.h"

//*********************************************************************************************

CSF::CSF(string field_name)
{
	this->field_name=field_name;
	total_energy = total_k2energy = 0.0;
}


/**********************************************************************************************

		 COPY CSF

**********************************************************************************************/


void CSF::Copy_to(CSF& to)
{
	to.F = F;
}


void CSF::Copy_to(PlainCSF& to)
{
	to.F = F;
}


void CSF::Copy_from(CSF& from)
{
	F = from.F;
}

void CSF::Copy_from(PlainCSF& from)
{
	F = from.F;
}

/**********************************************************************************************

		Forward transform(*F) -> *F
			-- Here temp_r is temporary array
		
		Forward_transform_transpose_order(*Ftr) -> *F 
		  -- Here Ftr unchanged.
	
***********************************************************************************************/

/* void CSF::Forward_transform()
{
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Forward_transform_array(F);		
} */

//
//

void CSF::Forward_transform(RSF rsf)
{
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Forward_transform(rsf.Fr, F);
}


/**********************************************************************************************

		Inverse transform(*F) -> *F
			-- Here temp_r is temporary array
		
		Inverse_transform_transpose_order(*F) -> *Ftr    

***********************************************************************************************/

/*void CSF::Inverse_transform()
{
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Inverse_transform_array(F);			
}*/


void CSF::Inverse_transform(RSF rsf)
{	
	global.program.sincostr_switch = sincostr_switch_F;
	universal->Inverse_transform(F, rsf.Fr);
}


/**********************************************************************************************

		  F(k) = F(k)/k^2

**********************************************************************************************/

void CSF::Divide_ksqr()
{
	universal->Array_divide_ksqr(F);
}


/**********************************************************************************************

		 Computes the total energy sum[|F(k)|^2/2]
			 dissipation sum[k^2|F(k)|^2]
			 Entropy

**********************************************************************************************/

void CSF::Compute_total_energy()
{
		
	total_energy = universal->Get_total_energy(F);
}



void CSF::Compute_total_k2energy()
{
	total_k2energy = 2 * universal->Get_total_Sn(F, 2);
}


void CSF::Compute_total_kn_energy(int n, DP &result)
{
	result = 2 * universal->Get_total_Sn(F, n);
}
//
//

void CSF::Compute_entropy()
{
	entropy = universal->Get_entropy_scalar(F);
}


/**********************************************************************************************

		 Modal energy

**********************************************************************************************/

DP CSF::Modal_energy(int i1, int i2, int i3)
{
	return universal->Modal_energy(i1, i2, i3, F);	
}



/**********************************************************************************************

		 Dealias CSF

**********************************************************************************************/

void CSF::Dealias()
{ 
	if (global.program.alias_option == "DEALIAS")
		universal->Dealias(F);
}


/**********************************************************************************************

		 Input and Output array

**********************************************************************************************/
// field_kind = T
void CSF::Write_complex_field()
{
 /*   if (global.program.decomposition == "SLAB")
        BasicIO::Write_complex_array(F, field_name+".F");
    
    else if (global.program.decomposition == "PENCIL") {
        if (global.program.basis_type == "SSF")
             BasicIO::Write_complex_array(F, field_name+".F");
       
        else {
            universal->Copy_array_complex_to_real(F, global.temp_array.Acopy);
            BasicIO::Write_complex_array(global.temp_array.Acopy, field_name+".F");
        }
    } */
}

void CSF::Write_reduced_complex_field()
{
 /*   if (global.program.decomposition == "SLAB")
        BasicIO::Write_reduced_complex_array(F, field_name+".F");
    
    else if (global.program.decomposition == "PENCIL") {
        if (global.program.basis_type == "SSF")
            BasicIO::Write_reduced_complex_array(F, field_name+".F");
        
        else {
            universal->Copy_array_complex_to_real(F, global.temp_array.Acopy);
            BasicIO::Write_reduced_complex_array(global.temp_array.Acopy, field_name+".F");
        }
    } */
}


void CSF::Read_complex_field()
{
/*    if (global.program.decomposition == "SLAB")
        BasicIO::Read_complex_array(F, field_name+".F");
    
    else if (global.program.decomposition == "PENCIL") {
        if (global.program.basis_type == "SSF")
            BasicIO::Read_complex_array(F, field_name+".F");
        
        else {
            BasicIO::Read_complex_array(global.temp_array.Acopy, field_name+".F");
            universal->Copy_array_real_to_complex(global.temp_array.Acopy, F);
        }
    } */
}

void CSF::Read_reduced_complex_field()
{
 /*   F = 0.0;
    
    if (global.program.decomposition == "SLAB")
        BasicIO::Read_reduced_complex_array(F, field_name+".F");
    
    else if (global.program.decomposition == "PENCIL") {
        if (global.program.basis_type == "SSF")
            BasicIO::Read_reduced_complex_array(F, field_name+".F");
        
        else {
            BasicIO::Read_reduced_complex_array(global.temp_array.Acopy, field_name+".F");
            universal->Copy_array_real_to_complex(global.temp_array.Acopy, F);
        }
    }
    
    universal->Adjust_array_after_reading_reduced(F);
  */
}



//******************************** End of CSF.cc **********************************************


