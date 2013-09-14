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

class Universal;


class BasicIO
{
	public:

	static hid_t H5T_FLOAT;
	static hid_t H5T_DOUBLE;

	static hid_t H5T_COMPLEX_FLOAT;
	static hid_t H5T_COMPLEX_DOUBLE;

	struct H5_dataset_meta {
		string name;
		vector<hsize_t> dimensions;
		H5T_class_t datatype;
	};

	struct H5_plan{
		hid_t dataspace;
		hid_t memspace;

		hid_t datatype;
	};

	template<int rank>
	struct Array_properties {
		TinyVector<int,rank> shape_full_complex_array;
		TinyVector<int,rank> shape_full_real_array;

		TinyVector<int,rank> id_complex_array;
		TinyVector<int,rank> id_real_array;

		TinyVector<int,rank> numprocs_complex_array;
		TinyVector<int,rank> numprocs_real_array;

		TinyVector<int,rank> shape_N_in_reduced;
		TinyVector<int,rank> shape_N_out_reduced;

		TinyVector<int,rank> Fourier_directions;
		int Z;

		hid_t datatype_complex_space;
		hid_t datatype_real_space;

		Array_properties():shape_N_in_reduced(0), shape_N_out_reduced(0){};
	};

	static string data_in_folder, data_out_folder;

	private:
		static string LOG_filename;
		static ofstream LOG_file;
		static time_t rawtime;

		//Required for Tarant-1 to Tarang-2 convertor
		static map<string, string> transition_table;
	
		static string prefix_prinfinity[1];
		static string prefix_fluid[3];
		static string prefix_fluid_T[4];
		static string prefix_mhd[6];
		static string prefix_mhd_T[7];
		static map<int,string*> dataset_name_old;
		//Convertor properties ends
	
		static herr_t err;
		static hid_t acc_template;

		// define an info object to store MPI-IO information 
		static MPI_Info FILE_INFO_TEMPLATE;
		static hid_t dataset;
		static hid_t file_identifier;

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
    
		template<typename Planner, int rank>
    	static void Set_H5_plans(Array_properties<rank> array_properties, Planner* planner);

    	static H5_plan Set_plan(int rank, int* my_id, int* numprocs, Array<int,1>* dataspace_filter, Array<int,1>* memspace_filter, hid_t datatype);
    

        static void Read(void* data, string dataset_name, H5_plan plan);
        static void Write(const void* data, string dataset_name, H5_plan plan, string folder_suffix="");
    
	
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

