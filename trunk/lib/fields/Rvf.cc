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
#include "Cvf.h"

//*********************************************************************************************

RVF::RVF(string field_name)
{
	this->field_name=field_name;
	
    V1r.resize(shape_real_array);
    V2r.resize(shape_real_array);
    V3r.resize(shape_real_array);

}

/**********************************************************************************************

				Forward transform: Inplace 
   
**********************************************************************************************/

/*void RVF::Forward_transform()
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Forward_transform_array(V1r);			
	
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Forward_transform_array(V2r);
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Forward_transform_array(V3r);			
} */

/**********************************************************************************************

			Forward_transform(*Vir) = Vi 
				Vir unchanged.
   
**********************************************************************************************/


void RVF::Forward_transform(CVF cvf)
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Forward_transform(V1r, cvf.V1);
	
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Forward_transform(V2r, cvf.V2);
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Forward_transform(V3r, cvf.V3);
}


/**********************************************************************************************

		Inplace Inverse Fourier transform.
   
**********************************************************************************************/

/* void RVF::Inverse_transform()
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Inverse_transform_array(V1r);			

	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Inverse_transform_array(V2r);	
	}

	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Inverse_transform_array(V3r);		
}
 */


/**********************************************************************************************

			Inverse_transform(Vi) = *Vir 
				Keeping Vi unchanged....
				temp = N1 x N2 x N3
   
**********************************************************************************************/

void RVF::Inverse_transform(CVF cvf)
{

	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Inverse_transform(cvf.V1, V1r);

	if (!global.program.two_dimension) { // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Inverse_transform(cvf.V2, V2r);

	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Inverse_transform(cvf.V3, V3r);
}


/*******************************************************************************

			 Outputs Vi in real space
 
 PS: In Pencil:SFF & SSF basis, z axis along the mid axis.  So complex data
 is jumbled up.  Therefore, complex<Ar> is copied to real<Ar> for writing.
 For reading, do the reverse.
   
********************************************************************************/
// field_kind = Ur or Wr
void RVF::Write_real_field()
{
    BasicIO::Write(V1r, field_name+".V1r", BasicIO::real);
	
/*	if (my_id==0)
		cout << "my_id , V1r : " << my_id << V1r << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_id==1)
		cout << "my_id , V1r : " << my_id << V1r << endl;
	MPI_Barrier(MPI_COMM_WORLD);
 */
	
    if (!global.program.two_dimension)
        BasicIO::Write(V2r, field_name+".V2r", BasicIO::real);
	
	BasicIO::Write(V3r, field_name+".V3r", BasicIO::real);
}


void RVF::Read_real_field()
{
    BasicIO::Read(V1r, field_name+".V1r", BasicIO::real);
	
    if (!global.program.two_dimension)
        BasicIO::Read(V2r, field_name+".V2r", BasicIO::real);
	
	BasicIO::Read(V3r, field_name+".V3r", BasicIO::real);
}

//**************************  End of RVF class definitions ************************************




