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


/*! \file four_ET.cc
 *
 * @sa four_ET.h
 *
 * @author M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug no known bug
 */

#include "BasicIO.h"
#include "Global_extern_vars.h"
#include <universal.h>

/*Assumption:
	0<= ap.Z < rank
*/
template<typename Plans, int rank>
void BasicIO::Set_H5_plans(Array_properties<rank> ap, Plans* plans)
{


	// if (master) cerr <<" Full: " << endl;
	plans->H5_full.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, h5::Select::all(3), ap.shape_full_complex_array, h5::Select::all(3), ap.datatype_complex_space);

	//Kz0_full
	TinyVector<hsize_t,rank> Kz0_filespace_shape = ap.shape_full_complex_array;
	Kz0_filespace_shape(2) = 1;


	// if (master) cerr <<" kz0 Full: " << endl;
	plans->H5_kz0_full.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, h5::Select(Range::all(), Range::all(), Range(0)), Kz0_filespace_shape, h5::Select::all(3), ap.datatype_complex_space);

	//Real
	// if (master) cerr <<" Real: " << endl;
	TinyVector<hsize_t,rank> Real_filespace_shape = shape(Nx, Ny, Nz);
	plans->H5_real.set_plan(ap.id_real_array, ap.numprocs_real_array, ap.shape_full_real_array, h5::Select(Range::all(), Range::all(), Range(0, Nz-1)), Real_filespace_shape, h5::Select::all(3), ap.datatype_real_space);

	// if (master) cerr <<" Slices: " << endl;
/*	plans->H5_slices.resize(global.io.slice_save.size());

	for (vector<H5Plan>::size_type i=0; i<plans->H5_slices.size(); i++) {
		plans->H5_slices[i].set_plan(MPI_COMMUNICATOR, ap.id_real_array, ap.numprocs_real_array, ap.datatype_real_space, ap.shape_full_real_array, global.io.slice_save[i]);
	}*/
/*
	 = Set_plan(rank, ap.id_real_array.data(), ap.numprocs_real_array.data(), dataspace_filter, memspace_filter, ap.datatype_real_space);
	if ( (plans->H5_real.dataspace == 0) || (plans->H5_real.memspace == 0) ){
		if (master)
			cerr << "BasicIO::Set_H5_plans - Error in setting HDF5 real IO plan." << endl;
		exit(1);
	}
*/

	//in_reduced
	if (sum(ap.shape_N_in_reduced)>0){
		if (all(ap.shape_full_complex_array>=ap.shape_N_in_reduced)) {

			h5::Expression e_mem_reduced;
			h5::Expression e_file_reduced;

			e_mem_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0,ap.shape_N_in_reduced(2)-1))
			              + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1), Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1))
						  + h5::Select(Range(Nx-ap.shape_N_in_reduced(0)/2,toEnd), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0,ap.shape_N_in_reduced(2)-1))
			              + h5::Select(Range(Nx-ap.shape_N_in_reduced(0)/2,toEnd), Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1), Range(0,ap.shape_N_in_reduced(1)/2-1)  , Range(0,ap.shape_N_in_reduced(2)-1))
			               + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1), Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1))
						   + h5::Select(Range(ap.shape_N_in_reduced(0)/2,toEnd), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0,ap.shape_N_in_reduced(2)-1))
			               + h5::Select(Range(ap.shape_N_in_reduced(0)/2,toEnd), Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1));


			// if (master) cerr <<" in_reduced:. " << endl;
			plans->H5_in_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, ap.shape_N_in_reduced, e_file_reduced, ap.datatype_complex_space);

			// if (master) cerr <<"In reduced kz0: " << endl;

			e_mem_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)     , Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0))
			              + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)     , Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0))
 			              + h5::Select(Range(Nx-ap.shape_N_in_reduced(0)/2,toEnd), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0))
			              + h5::Select(Range(Nx-ap.shape_N_in_reduced(0)/2,toEnd), Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)  , Range(0,ap.shape_N_in_reduced(1)/2-1)  , Range(0))
			               + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)  , Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0))
						   + h5::Select(Range(ap.shape_N_in_reduced(0)/2,toEnd), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0))
			               + h5::Select(Range(ap.shape_N_in_reduced(0)/2,toEnd), Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0));


			TinyVector<hsize_t,rank> N_in_reduced_kz0_shape = ap.shape_N_in_reduced;
			N_in_reduced_kz0_shape(2) = 1;

			//in_kz0_reduced
			// if (master) cerr <<" in_kz0_reduced: " << endl;
			plans->H5_in_kz0_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, N_in_reduced_kz0_shape, e_file_reduced, ap.datatype_complex_space);

			/*e_mem_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)-1), Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0,ap.shape_N_in_reduced(2)-1))
			              + h5::Select(Range(0,ap.shape_N_in_reduced(0)-1), Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)-1), Range(0,ap.shape_N_in_reduced(1)/2-1)  , Range(0,ap.shape_N_in_reduced(2)-1))
			               + h5::Select(Range(0,ap.shape_N_in_reduced(0)-1), Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0,ap.shape_N_in_reduced(2)-1));


			// if (master) cerr <<" in_reduced:. " << endl;
			plans->H5_in_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, ap.shape_N_in_reduced, e_file_reduced, ap.datatype_complex_space);

			// if (master) cerr <<"In reduced kz0: " << endl;

			e_mem_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)     , Range(0,ap.shape_N_in_reduced(1)/2-1)     , Range(0))
			              + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)     , Range(Ny-ap.shape_N_in_reduced(1)/2,toEnd), Range(0));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)  , Range(0,ap.shape_N_in_reduced(1)/2-1)  , Range(0))
			               + h5::Select(Range(0,ap.shape_N_in_reduced(0)/2-1)  , Range(ap.shape_N_in_reduced(1)/2,toEnd), Range(0));


			TinyVector<hsize_t,rank> N_in_reduced_kz0_shape = ap.shape_N_in_reduced;
			N_in_reduced_kz0_shape(2) = 1;

			//in_kz0_reduced
			// if (master) cerr <<" in_kz0_reduced: " << endl;
			plans->H5_in_kz0_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, N_in_reduced_kz0_shape, e_file_reduced, ap.datatype_complex_space);*/
		}
		else {
			if (master)
				cerr << "BasicIO::Set_H5_plans: For each axis, the condition 'field.N >= io.N_in_reduced' must hold." << endl;
			exit(1);
		}
	}

	//out_reduced
	if (sum(ap.shape_N_out_reduced)>0){
		if (all(ap.shape_full_complex_array>=ap.shape_N_out_reduced)) {

			h5::Expression e_mem_reduced;
			h5::Expression e_file_reduced;
			if (master) cerr <<"Out reduced: " << endl;

			e_mem_reduced = h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(0,ap.shape_N_out_reduced(1)/2-1)     , Range(0,ap.shape_N_out_reduced(2)/2))
			              + h5::Select(Range(Nx-ap.shape_N_out_reduced(0)/2,toEnd), Range(0,ap.shape_N_out_reduced(1)/2-1)     , Range(0,ap.shape_N_out_reduced(2)/2))
			              + h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(Ny-ap.shape_N_out_reduced(1)/2,toEnd), Range(0,ap.shape_N_out_reduced(2)/2))
			              + h5::Select(Range(Nx-ap.shape_N_out_reduced(0)/2,toEnd), Range(Ny-ap.shape_N_out_reduced(1)/2,toEnd), Range(0,ap.shape_N_out_reduced(2)/2));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(0,ap.shape_N_out_reduced(1)/2-1)  , Range(0,ap.shape_N_out_reduced(2)/2))
			               + h5::Select(Range(ap.shape_N_out_reduced(0)/2,toEnd), Range(0,ap.shape_N_out_reduced(1)/2-1)  , Range(0,ap.shape_N_out_reduced(2)/2))
			               + h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(ap.shape_N_out_reduced(1)/2,toEnd), Range(0,ap.shape_N_out_reduced(2)/2))
			               + h5::Select(Range(ap.shape_N_out_reduced(0)/2,toEnd)   , Range(ap.shape_N_out_reduced(1)/2,toEnd), Range(0,ap.shape_N_out_reduced(2)/2));


			// if (master) cerr <<" out_reduced:. " << endl;
			plans->H5_out_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, ap.shape_N_out_reduced, e_file_reduced, ap.datatype_complex_space);


			if (master) cerr <<"Out reduced kz0: " << endl;

			e_mem_reduced = h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(0,ap.shape_N_out_reduced(1)/2-1)     , Range(0))
			              + h5::Select(Range(Nx-ap.shape_N_out_reduced(0)/2,toEnd), Range(0,ap.shape_N_out_reduced(1)/2-1)     , Range(0))
			              + h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)     , Range(Ny-ap.shape_N_out_reduced(1)/2,toEnd), Range(0))
			              + h5::Select(Range(Nx-ap.shape_N_out_reduced(0)/2,toEnd), Range(Ny-ap.shape_N_out_reduced(1)/2,toEnd), Range(0));

			e_file_reduced = h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)  , Range(0,ap.shape_N_out_reduced(1)/2-1)  , Range(0))
			               + h5::Select(Range(ap.shape_N_out_reduced(0)/2,toEnd), Range(0,ap.shape_N_out_reduced(1)/2-1)  , Range(0))
			               + h5::Select(Range(0,ap.shape_N_out_reduced(0)/2-1)  , Range(ap.shape_N_out_reduced(1)/2,toEnd), Range(0))
			               + h5::Select(Range(ap.shape_N_out_reduced(0)/2,toEnd), Range(ap.shape_N_out_reduced(1)/2,toEnd), Range(0));


			TinyVector<hsize_t,rank> N_out_reduced_kz0_shape = ap.shape_N_out_reduced;
			N_out_reduced_kz0_shape(2) = 1;

			//out_kz0_reduced
			// if (master) cerr <<" out_kz0_reduced: " << endl;
			plans->H5_out_kz0_reduced.set_plan(ap.id_complex_array, ap.numprocs_complex_array, ap.shape_full_complex_array, e_mem_reduced, N_out_reduced_kz0_shape, e_file_reduced, ap.datatype_complex_space);

		}
		else {
			if (master)
				cerr << "BasicIO::Set_H5_plans: For each axis, the condition 'field.N >= io.N_out_reduced' must hold." << endl;
			exit(1);
		}
	}

}
