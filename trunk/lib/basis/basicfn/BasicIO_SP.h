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

	Array<int,1> filespace_filter[rank];
	Array<int,1> memoryspace_filter[rank];

	//Full
	for (int i=0; i<rank; i++){
		filespace_filter[i].resize(ap.shape_full_complex_array(i));
		filespace_filter[i]=1;

		memoryspace_filter[i].resize(ap.shape_full_complex_array(i));
		memoryspace_filter[i]=1;

	}

	//if (master) cerr <<" Full: " << endl;
	plans->H5_full.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
	if ( (plans->H5_full.filespace() == 0) || (plans->H5_full.memoryspace() == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 full array IO plan." << endl;
		exit(1);
	}


	//Kz0_full
	int kz0_Z_length=0;
	if (H5Tequal(ap.datatype_complex_space, h5::Dtype(H5Complex))>0)
		kz0_Z_length = 1;
	else if (H5Tequal(ap.datatype_complex_space, h5::Dtype(H5Real))>0)
		kz0_Z_length = 2;
	else{
		if (master)
			cerr << "BasicIO: Unrecognized datatype in 'ap.datatype_complex_space'" << endl;
		exit(0);
	}

	filespace_filter[ap.Z].resize(kz0_Z_length);
	filespace_filter[ap.Z]=1;
	memoryspace_filter[ap.Z]=0;
	memoryspace_filter[ap.Z](Range(0,kz0_Z_length-1))=1;

	//if (master) cerr <<" kz0 Full: " << endl;
	plans->H5_kz0_full.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
	if ( (plans->H5_kz0_full.filespace() == 0) || (plans->H5_kz0_full.memoryspace() == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 kz0 IO plan." << endl;
		exit(1);
	}


	//Real
	for (int i=0; i<rank; i++) {
		filespace_filter[i].resize(ap.shape_full_real_array(i));
		filespace_filter[i]=1;
		// filespace_filter[i](ap.shape_full_real_array(i))=1;
		memoryspace_filter[i].resize(ap.shape_full_real_array(i));
		memoryspace_filter[i]=1;
	}

	//if (master) cerr <<" Real: " << endl;
	plans->H5_real.set_plan(rank, ap.id_real_array.data(), ap.numprocs_real_array.data(), filespace_filter, memoryspace_filter, ap.datatype_real_space);
	if ( (plans->H5_real.filespace() == 0) || (plans->H5_real.memoryspace() == 0) ){
		if (master)
			cerr << "BasicIO: Error in setting HDF5 real IO plan." << endl;
		exit(1);
	}


	//in_reduced
	if (sum(ap.shape_N_in_reduced)>0){
		if (all(ap.shape_full_complex_array>=ap.shape_N_in_reduced)) {
			for (int i=0; i<rank; i++){
				filespace_filter[i].resize(ap.shape_N_in_reduced(i));
				filespace_filter[i]=0;

				memoryspace_filter[i].resize(ap.shape_full_complex_array(i));
				memoryspace_filter[i]=0;

				if (ap.Fourier_directions(i) && i!=ap.Z){ //Fourier axis
					filespace_filter[i](Range(0,ap.shape_N_in_reduced(i)/2-1))=1;
					filespace_filter[i](Range(ap.shape_N_in_reduced(i)/2+1, toEnd))=1;

					memoryspace_filter[i](Range(0,ap.shape_N_in_reduced(i)/2-1))=1;
					memoryspace_filter[i](Range(ap.shape_full_complex_array(i)-ap.shape_N_in_reduced(i)/2+1, toEnd))=1;
				}
				else{
					filespace_filter[i]=1;

					memoryspace_filter[i](Range(0,ap.shape_N_in_reduced(i)-1))=1;			
				}
			}

			//if (master) cerr <<" in_reduced:. " << endl;
			plans->H5_in_reduced.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
			if ( (plans->H5_in_reduced.filespace() == 0) || (plans->H5_in_reduced.memoryspace() == 0) ) {
				if (master)
					cerr << "BasicIO: Error in setting HDF5 in_reduced IO plan." << endl;
				exit(1);
			}


			//in_kz0_reduced
			filespace_filter[ap.Z].resize(kz0_Z_length);
			filespace_filter[ap.Z]=1;

			memoryspace_filter[ap.Z].resize(ap.shape_full_complex_array(ap.Z));
			memoryspace_filter[ap.Z]=0;
			memoryspace_filter[ap.Z](Range(0,kz0_Z_length-1))=1;

			//if (master) cerr <<" in__kz0_reduced: " << endl;
			plans->H5_in_kz0_reduced.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
			if ( (plans->H5_in_kz0_reduced.filespace() == 0) || (plans->H5_in_kz0_reduced.memoryspace() == 0) ){
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
	if (sum(ap.shape_N_out_reduced)>0) {
		if (all(ap.shape_full_complex_array>=ap.shape_N_in_reduced)) {
			for (int i=0; i<rank; i++){
				filespace_filter[i].resize(ap.shape_N_out_reduced(i));
				filespace_filter[i]=0;

				memoryspace_filter[i].resize(ap.shape_full_complex_array(i));
				memoryspace_filter[i]=0;

				if (ap.Fourier_directions(i) && i!=ap.Z){ //Fourier axis
					filespace_filter[i](Range(0,ap.shape_N_out_reduced(i)/2-1))=1;
					filespace_filter[i](Range(ap.shape_N_out_reduced(i)/2+1, toEnd))=1;

					memoryspace_filter[i](Range(0,ap.shape_N_out_reduced(i)/2-1))=1;
					memoryspace_filter[i](Range(ap.shape_full_complex_array(i)-ap.shape_N_out_reduced(i)/2+1, toEnd))=1;
				}
				else{
					filespace_filter[i]=1;
					memoryspace_filter[i](Range(0,ap.shape_N_out_reduced(i)-1))=1;			
				}
			}

			//if (master) cerr <<" out_reduced:. " << endl;
			plans->H5_out_reduced.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
			if ( (plans->H5_out_reduced.filespace() == 0) || (plans->H5_out_reduced.memoryspace() == 0) ){
				if (master)
					cerr << "BasicIO: Error in setting HDF5 out_reduced IO plan." << endl;
				exit(1);
			}

			//out_kz0_reduced
			filespace_filter[ap.Z].resize(kz0_Z_length);
			filespace_filter[ap.Z]=1;

			memoryspace_filter[ap.Z].resize(ap.shape_full_complex_array(ap.Z));
			memoryspace_filter[ap.Z]=0;
			memoryspace_filter[ap.Z](Range(0,kz0_Z_length-1))=1;

			//if (master) cerr <<" out__kz0_reduced:. " << endl;
			plans->H5_out_kz0_reduced.set_plan(rank, ap.id_complex_array.data(), ap.numprocs_complex_array.data(), filespace_filter, memoryspace_filter, ap.datatype_complex_space);
			if ( (plans->H5_out_kz0_reduced.filespace() == 0) || (plans->H5_out_kz0_reduced.memoryspace() == 0) ){
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
