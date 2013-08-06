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

int Get_num_high_segments(Array<int,1> filter, int id, int numprocs){
	int high_segments=0;

	int length_per_proc = filter.size() / numprocs;

	for (int n=0; n<id; n+=length_per_proc){
		for (int i=0; i<length_per_proc; i++)
			if (filter(n*length_per_proc+i)==1){
				high_segments++;
				break;
			}
	}
	
	return high_segments;
}

Filter_prop Get_prop(Array<int,1> filter){
	Filter_prop prop = {0,0,0};

	//Set first_rise
	for (int i=0; i<filter.size(); i++){
		if (filter(i)==1){
			prop.first_rise=i;
			break;
		}
	}

	//Set first_fall
	int first_one;
	for (int i=0; i<filter.size(); i++)
		if (filter(i)==1){
			first_one=i;
			break;
		}

	for (int i=first_one; i<filter.size(); i++)
		if (filter(i)==1)
			prop.first_fall=i;
		else
			break;

	if (prop.first_fall>prop.first_rise)
		prop.first_fall++;


	//Set zero_length
	if (filter(0)==1){
		for (int i=prop.first_fall; i<filter.size(); i++)
			if (filter(i)==0)
				prop.zero_length++;
			else
				break;
	}
	else{
		for (int i=0; i<filter.size(); i++)
			if (filter(i)==0)
				prop.zero_length++;
			else
				break;
	}

	return prop;
}

void BasicIO::Set_parameters()
{

	TinyVector<int, 3> shape_complex_array;
	TinyVector<int, 3> shape_real_array;
	TinyVector<int, 3> numprocs_all_dirn;
	TinyVector<int, 3> direction_z_complex_array_complement;

	shape_complex_array = global.field.shape_complex_array*shape(1,1,2);
	shape_real_array = global.field.shape_real_array;

	shape_full_complex_array*=shape(1,1,2);
	
	shape_in_reduced_array*=shape(1,1,2);
	shape_out_reduced_array*=shape(1,1,2);
	
	direction_z_complex_array_complement = (direction_z_complex_array^1);
	direction_z_complex_array*=shape(1,1,2);

	numprocs_all_dirn = shape_full_complex_array/shape_complex_array;


	// full array
	//reading
	// DATASPACE (On file)
	full.in.dataspace.dimension = shape_full_complex_array;
	full.in.dataspace.blockdim = shape_complex_array;
	full.in.dataspace.start = id_complex_array*shape_complex_array;
	full.in.dataspace.blockcount = 1,1,1;
	full.in.dataspace.blockstride = 1,1,1;
	
	// MEM SPACE (In Array)
	full.in.memspace.dimension = shape_complex_array;
	full.in.memspace.start = 0,0,0;
	full.in.memspace.blockcount = 1,1,1;
	full.in.memspace.blockstride = 1,1,1;
	full.in.memspace.blockdim = shape_complex_array;

	//writing
	full.out = full.in;


	// for kz=0 plane
	//reading
	kz0_full.in.dataspace.dimension = shape_full_complex_array*direction_z_complex_array_complement + direction_z_complex_array;
	kz0_full.in.dataspace.blockdim = shape_complex_array*direction_z_complex_array_complement + direction_z_complex_array;
	kz0_full.in.dataspace.start = id_complex_array*shape_complex_array;
	kz0_full.in.dataspace.blockcount = 1,1,1;
	kz0_full.in.dataspace.blockstride = 1,1,1;

	
	kz0_full.in.memspace.dimension = shape_complex_array;
	kz0_full.in.memspace.blockdim = shape_complex_array*direction_z_complex_array_complement + direction_z_complex_array;
	kz0_full.in.memspace.start = 0,0,0;
	kz0_full.in.memspace.blockcount = 1,1,1;
	kz0_full.in.memspace.blockstride = 1,1,1;
	
	// writing
	kz0_full.out = kz0_full.in;


	//*****************************
	// for real space data
	//*****************************
	
	real.in.dataspace.dimension = shape_full_real_array;
	real.in.dataspace.blockdim = shape_real_array;
	real.in.dataspace.start = id_real_array*shape_real_array;
	real.in.dataspace.blockcount = 1,1,1;
	real.in.dataspace.blockstride = 1,1,1;
	
	real.in.memspace.dimension = shape_real_array;
	real.in.memspace.blockdim = shape_real_array;
	real.in.memspace.start = 0,0,0;
	real.in.memspace.blockcount = 1,1,1;
	real.in.memspace.blockstride = 1,1,1;

	real.out = real.in;
	//********************

	//reduced
	
	Filter_prop prop;
	int num_high_segments;
	
	Array<int,1> filter_reduced[3];
	Array<int,1> my_filter_reduced[3];
	
	filter_reduced[0].resize(shape_full_complex_array(0));
	filter_reduced[1].resize(shape_full_complex_array(1));
	filter_reduced[2].resize(shape_full_complex_array(2));
	
	my_filter_reduced[0].resize(shape_complex_array(0));
	my_filter_reduced[1].resize(shape_complex_array(1));
	my_filter_reduced[2].resize(shape_complex_array(2));

	if (global.io.N_in_reduced.size()==3){
		reduced.in.dataspace.dimension = shape_in_reduced_array;
		reduced.in.memspace.dimension = shape_complex_array;
		reduced.out.dataspace.dimension = shape_out_reduced_array;
		reduced.out.memspace.dimension = shape_complex_array;


		for (int i=0; i<3; i++){
			filter_reduced[i]=0;

			// negative k directions
			if (Fourier_directions(i)*direction_z_complex_array_complement(i)){		//Fourier and not Z
				filter_reduced[i](Range(0,shape_in_reduced_array(i)/2-1))=1;
				filter_reduced[i](Range(shape_full_complex_array(i)-shape_in_reduced_array(i)/2, toEnd))=1;
			}
			else  {// positive k directions
				filter_reduced[i](Range(0,shape_in_reduced_array(i)-1))=1;
			}

			my_filter_reduced[i]=filter_reduced[i](Range(id_complex_array(i)*shape_complex_array(i),(id_complex_array(i)+1)*shape_complex_array(i)-1));
		}

	

		for (int i=0; i<3; i++){
			prop = Get_prop(my_filter_reduced[i]);
			num_high_segments = Get_num_high_segments(filter_reduced[i], id_complex_array(i), numprocs_all_dirn(i));

			reduced.in.memspace.blockdim(i) = prop.first_fall - prop.first_rise;
			reduced.in.memspace.start(i) = prop.first_rise;
			

			if (Fourier_directions(i)*direction_z_complex_array_complement(i)) {
				if (numprocs_all_dirn(i)==1){
					reduced.in.memspace.blockcount(i)=2;
					reduced.in.memspace.blockstride(i) = prop.first_fall - prop.first_rise + prop.zero_length;
				}
				else {
					if (prop.first_fall - prop.first_rise == 0) {
						reduced.in.memspace.blockcount(i)=0;
						reduced.in.memspace.blockstride(i) = 1;
					}
					else {
						reduced.in.memspace.blockcount(i)=1;
						reduced.in.memspace.blockstride(i)=1;
					}
					
				}
			}
			else {
				if (prop.first_fall - prop.first_rise == 0) {
					reduced.in.memspace.blockcount(i)=0;
					reduced.in.memspace.blockstride(i) = 1;
				}
				else {
					reduced.in.memspace.blockcount(i)=1;
					reduced.in.memspace.blockstride(i)=1;
				}
			}
				
			reduced.in.dataspace.blockdim(i) = prop.first_fall - prop.first_rise;
			reduced.in.dataspace.start(i) = num_high_segments*(reduced.in.dataspace.blockdim(i));
			reduced.in.dataspace.blockcount(i) = reduced.in.memspace.blockcount(i);
			
			if (Fourier_directions(i)*direction_z_complex_array_complement(i) && (numprocs_all_dirn(i)==1))
				reduced.in.dataspace.blockstride(i) = shape_in_reduced_array(i)/2;
			else
				reduced.in.dataspace.blockstride(i) = 1;
		}
		
		kz0_reduced.in = reduced.in;
		
		kz0_reduced.in.dataspace.dimension(1) = 1;
		kz0_reduced.in.dataspace.blockdim(1) = 1;
		
		kz0_reduced.in.memspace.blockdim(1) = 1;
		
	}
	
	
	if (global.io.N_out_reduced.size()==3){
		
		//Reduced out
		for (int i=0; i<3; i++){
			filter_reduced[i]=0;
			
			// negative k directions
			if (Fourier_directions(i)*direction_z_complex_array_complement(i)){		//Fourier and not Z
				filter_reduced[i](Range(0,shape_out_reduced_array(i)/2-1))=1;
				filter_reduced[i](Range(shape_full_complex_array(i)-shape_out_reduced_array(i)/2, toEnd))=1;

			}
			else  {// positive k directions
				filter_reduced[i](Range(0,shape_out_reduced_array(i)-1))=1;
			}
			
			my_filter_reduced[i]=filter_reduced[i](Range(id_complex_array(i)*shape_complex_array(i),(id_complex_array(i)+1)*shape_complex_array(i)-1));

		}
		
		
		for (int i=0; i<3; i++){
			prop = Get_prop(my_filter_reduced[i]);
			num_high_segments = Get_num_high_segments(filter_reduced[i], id_complex_array(i), numprocs_all_dirn(i));
			
			reduced.out.memspace.blockdim(i) = prop.first_fall - prop.first_rise;
			reduced.out.memspace.start(i) = prop.first_rise;
			
			
			if (Fourier_directions(i)*direction_z_complex_array_complement(i)) {
				if (numprocs_all_dirn(i)==1){
					reduced.out.memspace.blockcount(i)=2;
					reduced.out.memspace.blockstride(i) = prop.first_fall - prop.first_rise + prop.zero_length;
				}
				else {
					if (prop.first_fall - prop.first_rise == 0) {
						reduced.out.memspace.blockcount(i)=0;
						reduced.out.memspace.blockstride(i) = 1;
					}
					else {
						reduced.out.memspace.blockcount(i)=1;
						reduced.out.memspace.blockstride(i)=1;
					}
					
				}
			}
			else {
				if (prop.first_fall - prop.first_rise == 0) {
					reduced.out.memspace.blockcount(i)=0;
					reduced.out.memspace.blockstride(i) = 1;
				}
				else {
					reduced.out.memspace.blockcount(i)=1;
					reduced.out.memspace.blockstride(i)=1;
				}
			}
			

			reduced.out.dataspace.blockdim(i) = prop.first_fall - prop.first_rise;
			reduced.out.dataspace.start(i) = num_high_segments*(reduced.out.dataspace.blockdim(i));
			reduced.out.dataspace.blockcount(i) = reduced.out.memspace.blockcount(i);
			
			if (Fourier_directions(i)*direction_z_complex_array_complement(i) && (numprocs_all_dirn(i)==1))
				reduced.out.dataspace.blockstride(i) = shape_out_reduced_array(i)/2;
			else
				reduced.out.dataspace.blockstride(i) = 1;
		}


		kz0_reduced.out = reduced.out;
		
		kz0_reduced.out.dataspace.dimension(1) = 1;
		kz0_reduced.out.dataspace.blockdim(1) = 1;
		
		kz0_reduced.out.memspace.blockdim(1) = 1;
	}
}

