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

/**********************************************************************************************

					Inverse_transform_transpose_order(F) = *Fr 
						F unchanged
   
**********************************************************************************************/

void RSF::Inverse_transform(CSF csf)
{
	global.program.sincostr_switch = sincostr_switch_F;
	
	universal->Inverse_transform(csf.F, Fr);
}


// field_kind = Tr

void RSF::Read_real_field()
{
	Fr = 0;
	universal->Read(Fr, universal->H5_real, field_name+".Fr");
}

void RSF::Write_real_field()
{
	string folder_name="real_" + To_string(global.time.now);

	universal->Write(Fr, universal->H5_real , "w", folder_name, field_name+".Fr");
}

void RSF::Write_real_field_slice(unsigned int slice_file_counter)
{
	string file_name;

	for (vector<h5::Plan>::size_type i=0; i<universal->H5_slices.size(); i++) {
		file_name = "slice_" + To_string(i) + "_" + To_string(slice_file_counter);

		universal->Write(Fr, universal->H5_slices[i], "w", "slices", file_name, "Fr");
	}
}



//************************ END of RSF class Definitions ***************************************
