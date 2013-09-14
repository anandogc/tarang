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


/*! \file array_basic.cc
 *
 * @sa field_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */


#include "BasicIO.h"

//initialize prefixes
string BasicIO::prefix_prinfinity[1];
string BasicIO::prefix_fluid[3];
string BasicIO::prefix_fluid_T[4];
string BasicIO::prefix_mhd[6];
string BasicIO::prefix_mhd_T[7];
map<int,string*> BasicIO::dataset_name_old;

hsize_t BasicIO::start_3d[3];
hsize_t BasicIO::stride_3d[3];
hsize_t BasicIO::count_3d[3];

//**********************
//Tarang-1 to Tarang-2 field converter
//**********************

void BasicIO::Tarang1_Tarang2_convertor_init(){
	transition_table["CV1"]="U.V1";
	transition_table["CV2"]="U.V2";
	transition_table["CV3"]="U.V3";
	
	transition_table["CW1"]="B.V1";
	transition_table["CW2"]="B.V2";
	transition_table["CW3"]="B.V3";
	
	transition_table["T"]="T.F";
	
	//initialize prefixes
	prefix_prinfinity[0]="T";
	
	prefix_fluid[0]="CV1";
	prefix_fluid[1]="CV2";
	prefix_fluid[2]="CV3";
	
	prefix_fluid_T[0]="CV1";
	prefix_fluid_T[1]="CV2";
	prefix_fluid_T[2]="CV3";
	prefix_fluid_T[3]="T";
	
	prefix_mhd[0]="CV1";
	prefix_mhd[1]="CV2";
	prefix_mhd[2]="CV3";
	prefix_mhd[3]="CW1";
	prefix_mhd[4]="CW2";
	prefix_mhd[5]="CW3";
	
	prefix_mhd_T[0]="CV1";
	prefix_mhd_T[1]="CV2";
	prefix_mhd_T[2]="CV3";
	prefix_mhd_T[3]="CW1";
	prefix_mhd_T[4]="CW2";
	prefix_mhd_T[5]="CW3";
	prefix_mhd_T[6]="T";
	
	
	//assign them to map
	dataset_name_old[1]=prefix_prinfinity;
	dataset_name_old[3]=prefix_fluid;
	dataset_name_old[4]=prefix_fluid_T;
	dataset_name_old[6]=prefix_mhd;
	dataset_name_old[7]=prefix_mhd_T;
}

void BasicIO::CV_Input_HDF5(int number_of_fields, ...)
{
	va_list vl;
	va_start( vl, number_of_fields);
	
	file_name = data_in_folder+"/input.h5";
	file_identifier = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, acc_template);
	
	Array<complx, 3>* A3;
	
	for (int i=0; i<number_of_fields; i++){
	    A3 = reinterpret_cast<Array<complx, 3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    Read_DP_Array(dataset_name_old[number_of_fields][i], dataArray, local_Nx, Ny, Nz);
	}
	
	H5Fclose(file_identifier);
}

	//**************************************************************************************************
void BasicIO::CV_Input_plane_HDF5(int number_of_fields, ...)
{
	va_list vl;
	va_start( vl, number_of_fields);
	
	file_name = data_in_folder+"/input.h5";
	file_identifier = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, acc_template);
	
	Array<complx,2> temp_array(local_Nx, Ny);
	
	Array<complx, 3>* A3;
	
	for (int i=0; i<number_of_fields; i++){
		if (i==2 || i==5){
			dataArray = reinterpret_cast<DP *>(temp_array.data());
			Read_DP_Array(dataset_name_old[number_of_fields][i], dataArray, local_Nx, Ny, 2);
			A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
			(*A3)(Range::all(),Range::all(),0)=temp_array;
		}
		else
			{
			A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
			dataArray = reinterpret_cast<DP *>((*A3).data());
			Read_DP_Array(dataset_name_old[number_of_fields][i], dataArray, local_Nx, Ny, 2*(Nz/2+1));
			}
	}
	
	H5Fclose(file_identifier);
}


	//**************************************************************************************************

void BasicIO::CV_Output_HDF5(int number_of_fields, ...)
{
	
	va_list vl;
	va_start( vl, number_of_fields);
	
	Array<complx, 3>* A3;
	
	for (int i=0; i<number_of_fields; i++){
	    A3= reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    if (global.mpi.master) 
			cout << dataset_name_old[number_of_fields][i] << "  -->  " << transition_table[dataset_name_old[number_of_fields][i]] << endl;
	    Write_DP_Array(transition_table[dataset_name_old[number_of_fields][i]], dataArray, int(local_Nx), Ny, 2*(Nz/2+1));
	}
}

	//**************************************************************************************************

void BasicIO::CV_Output_plane_HDF5(int number_of_fields, ...)
{
	
	va_list vl;
	va_start( vl, number_of_fields);
	
	Array<complx,2> temp_array(local_Nx, Ny);
	Array<complx, 3>* A3;
	
	for (int i=0; i<number_of_fields; i++){
		if (i==2 || i==5){
			A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
			temp_array = (*A3)(Range::all(),Range::all(),0);
			dataArray = reinterpret_cast<DP *>(temp_array.data());
			if (global.mpi.master)
				cout << dataset_name_old[number_of_fields][i] << "  -->  " << transition_table[dataset_name_old[number_of_fields][i]] << endl;
			Write_DP_Array(transition_table[dataset_name_old[number_of_fields][i]]+"kz0", dataArray, int(local_Nx), Ny, 2);
		}
		else
			{
			A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
			dataArray = reinterpret_cast<DP *>((*A3).data());
			if (global.mpi.master)
				cout << dataset_name_old[number_of_fields][i] << "  -->  " << transition_table[dataset_name_old[number_of_fields][i]] << endl;
			Write_DP_Array(transition_table[dataset_name_old[number_of_fields][i]], dataArray, int(local_Nx), Ny, 2*(Nz/2+1));
			}
	}
}

void BasicIO::Read_DP_Array(string dataset_name, DP *dataArray, int dim1, int dim2, int dim3)
{
	hsize_t dimens_3d[3], start_3d[3], stride_3d[3], count_3d[3];
	
    dimens_3d[0] = (dim1)*numprocs;
    dimens_3d[1] = (dim2);
    dimens_3d[2] = (dim3);
	
    dataspace = H5Screate_simple(rank, dimens_3d, NULL);
	
    start_3d[0] = (dim1)*my_id;
    start_3d[1] = 0;
    start_3d[2] = 0;
	
    stride_3d[0] = 1;
    stride_3d[1] = 1;
    stride_3d[2] = 1;
	
    count_3d[0] = (dim1);
    count_3d[1] = (dim2);
    count_3d[2] = (dim3);
	
    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d, stride_3d, count_3d, NULL);
	
	
    dimens_3d[0] = (dim1);
    dimens_3d[1] = (dim2);
    dimens_3d[2] = (dim3);
	
    memspace = H5Screate_simple(rank, dimens_3d, NULL);
	
	
    start_3d[0] = 0;
    start_3d[1] = 0;
    start_3d[2] = 0;
	
    stride_3d[0] = 1;
    stride_3d[1] = 1;
    stride_3d[2] = 1;
	
    count_3d[0] = (dim1);
    count_3d[1] = (dim2);
    count_3d[2] = (dim3);
	
    err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
							  start_3d, NULL, count_3d, NULL);
	
	dataset = H5Dopen2(file_identifier, dataset_name.c_str(), H5P_DEFAULT);
	
	err = H5Dread(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, dataArray);
	
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}


	//**************************************************************************************************

void BasicIO::Write_DP_Array(string dataset_name, DP *dataArray, int dim1, int dim2, int dim3)
{
	
	hsize_t dimens_3d[3], start_3d[3], stride_3d[3], count_3d[3];
	
	
	dimens_3d[0] = (dim1)*numprocs;
    dimens_3d[1] = (dim2);
    dimens_3d[2] = (dim3);
	
	dataspace = H5Screate_simple(rank, dimens_3d, NULL);
	
	start_3d[0] = (dim1)*my_id;
	start_3d[1] = 0;
	start_3d[2] = 0;
	
	stride_3d[0] = 1;
	stride_3d[1] = 1;
	stride_3d[2] = 1;
	
	count_3d[0] = (dim1);
	count_3d[1] = (dim2);
	count_3d[2] = (dim3);
	
	err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d, stride_3d, count_3d, NULL);
	
	
	dimens_3d[0] = (dim1);
	dimens_3d[1] = (dim2);
	dimens_3d[2] = (dim3);
	
	memspace = H5Screate_simple(rank, dimens_3d, NULL);
	
	start_3d[0] = 0;
	start_3d[1] = 0;
	start_3d[2] = 0;
	
	stride_3d[0] = 1;
	stride_3d[1] = 1;
	stride_3d[2] = 1;
	
	count_3d[0] = (dim1);
	count_3d[1] = (dim2);
	count_3d[2] = (dim3);
	
	err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
							  start_3d, NULL, count_3d, NULL);
	
		//Open the file for writing
	file_name = data_out_folder + "/" + dataset_name + ".h5";
	
	file_identifier = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
	
	dataset = H5Dcreate2(file_identifier, dataset_name.c_str(), H5_DATATYPE,
						 dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	err = H5Dwrite(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, dataArray);
	
	H5Dclose(dataset);
	H5Fclose(file_identifier);
	H5Sclose(memspace);
	H5Sclose(dataspace);
	
}





