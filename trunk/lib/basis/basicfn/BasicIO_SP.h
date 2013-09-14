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
	0<= array_properties.Z < rank
*/
template<typename Planner, int rank>
void BasicIO::Set_H5_plans(Array_properties<rank> array_properties, Planner* planner)
{

	Array<int,1> dataspace_filter[rank];
	Array<int,1> memspace_filter[rank];

	//Full
	for (int i=0; i<rank; i++){
		dataspace_filter[i].resize(array_properties.shape_full_complex_array(i));
		dataspace_filter[i]=1;

		memspace_filter[i].resize(array_properties.shape_full_complex_array(i));
		memspace_filter[i]=1;

	}

	//if (master) cerr <<" Full: " << endl;
	planner->H5_full = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
	if ( (planner->H5_full.dataspace == 0) || (planner->H5_full.memspace == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 full array IO plan." << endl;
		exit(1);
	}


	//Kz0_full
	int kz0_Z_length=0;
	if (H5Tequal(array_properties.datatype_complex_space, H5T_COMPLX)>0)
		kz0_Z_length = 1;
	else if (H5Tequal(array_properties.datatype_complex_space, H5T_DP)>0)
		kz0_Z_length = 2;
	else{
		if (master)
			cerr << "BasicIO: Unrecognized datatype in 'Array_properties.datatype_complex_space'" << endl;
		exit(0);
	}

	dataspace_filter[array_properties.Z].resize(kz0_Z_length);
	dataspace_filter[array_properties.Z]=1;
	memspace_filter[array_properties.Z]=0;
	memspace_filter[array_properties.Z](Range(0,kz0_Z_length-1))=1;

	//if (master) cerr <<" kz0 Full: " << endl;
	planner->H5_kz0_full = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
	if ( (planner->H5_kz0_full.dataspace == 0) || (planner->H5_kz0_full.memspace == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 kz0 IO plan." << endl;
		exit(1);
	}



	//Real
	for (int i=0; i<rank; i++){
		dataspace_filter[i].resize(array_properties.shape_full_real_array(i));
		dataspace_filter[i]=1;
		memspace_filter[i].resize(array_properties.shape_full_real_array(i));
		memspace_filter[i]=1;
	}

	//if (master) cerr <<" Real: " << endl;
	planner->H5_real = Set_plan(rank, array_properties.id_real_array.data(), array_properties.numprocs_real_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_real_space);
	if ( (planner->H5_real.dataspace == 0) || (planner->H5_real.memspace == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 real IO plan." << endl;
		exit(1);
	}


	//in_reduced
	if (sum(array_properties.shape_N_in_reduced)>0){
		if (all(array_properties.shape_full_complex_array>=array_properties.shape_N_in_reduced)) {
			for (int i=0; i<rank; i++){
				dataspace_filter[i].resize(array_properties.shape_N_in_reduced(i));
				dataspace_filter[i]=0;

				memspace_filter[i].resize(array_properties.shape_full_complex_array(i));
				memspace_filter[i]=0;

				if (array_properties.Fourier_directions(i) && i!=array_properties.Z){ //Fourier axis
					dataspace_filter[i](Range(0,array_properties.shape_N_in_reduced(i)/2-1))=1;
					dataspace_filter[i](Range(array_properties.shape_N_in_reduced(i)/2+1, toEnd))=1;

					memspace_filter[i](Range(0,array_properties.shape_N_in_reduced(i)/2-1))=1;
					memspace_filter[i](Range(array_properties.shape_full_complex_array(i)-array_properties.shape_N_in_reduced(i)/2+1, toEnd))=1;
				}
				else{
					dataspace_filter[i]=1;

					memspace_filter[i](Range(0,array_properties.shape_N_in_reduced(i)-1))=1;			
				}
			}

			//if (master) cerr <<" in_reduced:. " << endl;
			planner->H5_in_reduced = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
			if ( (planner->H5_in_reduced.dataspace == 0) || (planner->H5_in_reduced.memspace == 0) ) {
				if (master)
					cerr << "BasicIO: Error in setting HDF5 in_reduced IO plan." << endl;
				exit(1);
			}


			//in_kz0_reduced
			dataspace_filter[array_properties.Z].resize(kz0_Z_length);
			dataspace_filter[array_properties.Z]=1;

			memspace_filter[array_properties.Z].resize(array_properties.shape_full_complex_array(array_properties.Z));
			memspace_filter[array_properties.Z]=0;
			memspace_filter[array_properties.Z](Range(0,kz0_Z_length-1))=1;

			//if (master) cerr <<" in__kz0_reduced: " << endl;
			planner->H5_in_kz0_reduced = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
			if ( (planner->H5_in_kz0_reduced.dataspace == 0) || (planner->H5_in_kz0_reduced.memspace == 0) ){
				if (master)
					cerr << "BasicIO: Error in setting HDF5 in_kz0_reduced IO plan." << endl;
				exit(1);
			}
		}
		else {
			if (master)
				cerr << "BasicIO::Set_H5_plans: For each axis, the condition 'field.N >= io.N_in_reduced' must hold." << endl;
			exit(1);
		}
	}


	//out_reduced
	if (sum(array_properties.shape_N_out_reduced)>0){
		if (all(array_properties.shape_full_complex_array>=array_properties.shape_N_in_reduced)) {
			for (int i=0; i<rank; i++){
				dataspace_filter[i].resize(array_properties.shape_N_out_reduced(i));
				dataspace_filter[i]=0;

				memspace_filter[i].resize(array_properties.shape_full_complex_array(i));
				memspace_filter[i]=0;

				if (array_properties.Fourier_directions(i) && i!=array_properties.Z){ //Fourier axis
					dataspace_filter[i](Range(0,array_properties.shape_N_out_reduced(i)/2-1))=1;
					dataspace_filter[i](Range(array_properties.shape_N_out_reduced(i)/2+1, toEnd))=1;

					memspace_filter[i](Range(0,array_properties.shape_N_out_reduced(i)/2-1))=1;
					memspace_filter[i](Range(array_properties.shape_full_complex_array(i)-array_properties.shape_N_out_reduced(i)/2+1, toEnd))=1;
				}
				else{
					dataspace_filter[i]=1;
					memspace_filter[i](Range(0,array_properties.shape_N_out_reduced(i)-1))=1;			
				}
			}

			//if (master) cerr <<" out_reduced:. " << endl;
			planner->H5_out_reduced = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
			if ( (planner->H5_out_reduced.dataspace == 0) || (planner->H5_out_reduced.memspace == 0) ){
				if (master)
					cerr << "BasicIO: Error in setting HDF5 out_reduced IO plan." << endl;
				exit(1);
			}

			//out_kz0_reduced
			dataspace_filter[array_properties.Z].resize(kz0_Z_length);
			dataspace_filter[array_properties.Z]=1;

			memspace_filter[array_properties.Z].resize(array_properties.shape_full_complex_array(array_properties.Z));
			memspace_filter[array_properties.Z]=0;
			memspace_filter[array_properties.Z](Range(0,kz0_Z_length-1))=1;

			//if (master) cerr <<" out__kz0_reduced:. " << endl;
			planner->H5_out_kz0_reduced = Set_plan(rank, array_properties.id_complex_array.data(), array_properties.numprocs_complex_array.data(), dataspace_filter, memspace_filter, array_properties.datatype_complex_space);
			if ( (planner->H5_out_kz0_reduced.dataspace == 0) || (planner->H5_out_kz0_reduced.memspace == 0) ){
				if (master)
					cerr << "BasicIO: Error in setting HDF5 out_kz0_reduced IO plan." << endl;
				exit(1);
			}

		}
		else {
			if (master)
				cerr << "BasicIO::Set_H5_plans: For each axis, the condition 'field.N >= io.N_out_reduced' must hold." << endl;
			exit(1);
		}
	}
}
