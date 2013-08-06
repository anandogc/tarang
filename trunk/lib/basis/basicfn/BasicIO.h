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


/*! \file  IO.h
 * 
 *	@brief  Class declaration of input/output function.
 * 
 *
 *	@author  A. G. Chatterjee
 *	@version 2.0 MPI
 *	@date feb 2012
 *
 * @bug  No known bugs
 */

//*********************************************************************************************
#ifndef _H_BASIC_IO
#define _H_BASIC_IO

#include <hdf5.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "def_vars.h"
#include "Global_extern_vars.h"
// #include "basis.h"

struct H5_dataset_meta {
	string name;
	vector<hsize_t> dimensions;
	H5T_class_t datatype;
};

struct H5_space{
	TinyVector<hsize_t, 3> dimension;
	TinyVector<hsize_t, 3> start;
	TinyVector<hsize_t, 3> blockstride;
	TinyVector<hsize_t, 3> blockcount;
	TinyVector<hsize_t, 3> blockdim;

};


struct H5_space_group{
	H5_space dataspace;
	H5_space memspace;

};


struct H5_IO {
	H5_space_group in;
	H5_space_group out;

};

struct Filter_prop {
	int first_rise;
	int first_fall;
	int zero_length;
};

ostream& operator<< (ostream &out, H5_space &h5_space);
bool operator== (H5_space &h5_space1, H5_space &h5_space2);

ostream& operator<< (ostream &out, H5_space_group &h5_space_group);
bool operator== (H5_space_group &h5_space_group1, H5_space_group &h5_space_group2);

ostream& operator<< (ostream &out, H5_IO &h5_io);
bool operator== (H5_IO &h5_io1, H5_IO &h5_io2);


ostream& operator<< (ostream &out, Filter_prop &prop);

class BasicIO
{
	public:

	static string data_in_folder, data_out_folder;
	
	static TinyVector<int, 3> id_complex_array;
	static TinyVector<int, 3> id_real_array;
	
	
	static TinyVector<int, 3> shape_complex_array;
	static TinyVector<int, 3> shape_real_array;
	
	static TinyVector<int, 3> shape_full_complex_array;
	static TinyVector<int, 3> shape_full_real_array;
	
	static TinyVector<int, 3> direction_z_complex_array;
	static TinyVector<int, 3> direction_z_real_array;
	
	static TinyVector<int, 3> shape_in_reduced_array;
	static TinyVector<int, 3> shape_out_reduced_array;
	
	static TinyVector<int, 3> Fourier_directions;

    
	static H5_IO full;
	static H5_IO kz0_full;
	static H5_IO reduced;
	static H5_IO kz0_reduced;
	static H5_IO real;  
    

	private:
		static string LOG_filename;
		static ofstream LOG_file;
		static time_t rawtime;

		static map<string, string> transition_table;
	
		static string prefix_prinfinity[1];
		static string prefix_fluid[3];
		static string prefix_fluid_T[4];
		static string prefix_mhd[6];
		static string prefix_mhd_T[7];
		static map<int,string*> dataset_name_old;
	
		static hid_t err;
		static int rank;
		static hsize_t start_3d[3];
		static hsize_t stride_3d[3];
		static hsize_t count_3d[3];
		
		static hid_t dataspace, memspace, dataset;
		static hid_t acc_template;

		// define an info object to store MPI-IO information 
		static MPI_Info FILE_INFO_TEMPLATE;
		static hid_t file_identifier;
		static int ierr;

		static int dim1, dim2, dim3;
		static DP *dataArray;

		static string folder_name, file_name;
		static bool frequent;

	public:

		static void Initialize();
		static void Finalize();
		static void Log(string message);
		static void Begin_frequent();
		static void End_frequent();
	
		static vector<H5_dataset_meta> Get_meta(string filename);
    
    	static void Set_parameters();
    

        static void Read(Array<complx,3> A, string dataset_name, H5_IO config);
        static void Read(Array<DP,3> A, string dataset_name, H5_IO config);
        static void Read_basic(DP* data, string dataset_name, H5_IO config);

        static void Write(Array<complx,3> A, string dataset_name, H5_IO config);
        static void Write(Array<DP,3> A, string dataset_name, H5_IO config);
        static void Write_basic(DP* data, string dataset_name, H5_IO config);
    
	
		//Tarang-1 to Tarang-2 field converter
        static void CV_Input_HDF5(int number_of_fields, ...);
        static void CV_Input_plane_HDF5(int number_of_fields, ...);
        static void CV_Output_HDF5(int number_of_fields, ...);
        static void CV_Output_plane_HDF5(int number_of_fields, ...);
        static void Read_DP_Array(string dataset_name, DP *dataArray, int dim1, int dim2, int dim3);
        static void Write_DP_Array(string dataset_name, DP *dataArray, int dim1, int dim2, int dim3);
};


#endif
//========================= Class declaration of BasicIO ends ============================== 

