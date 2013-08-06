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

#include "BasicIO.h"
//*********************************************************************************************

#define TERMINATE(function, message) \
{	\
	cerr << __FILE__ << ": In function '" << function << "': " << __LINE__ << message << endl;\
	BasicIO::Finalize();\
	MPI_Finalize();\
}


//public:
string BasicIO::data_in_folder="/in/";     
string BasicIO::data_out_folder="/out/";

TinyVector<int, 3> BasicIO::id_complex_array;
TinyVector<int, 3> BasicIO::id_real_array;


TinyVector<int, 3> BasicIO::shape_complex_array;
TinyVector<int, 3> BasicIO::shape_real_array;

TinyVector<int, 3> BasicIO::shape_full_complex_array;
TinyVector<int, 3> BasicIO::shape_full_real_array;

TinyVector<int, 3> BasicIO::direction_z_complex_array;
TinyVector<int, 3> BasicIO::direction_z_real_array;

TinyVector<int, 3> BasicIO::shape_in_reduced_array;
TinyVector<int, 3> BasicIO::shape_out_reduced_array;

TinyVector<int, 3> BasicIO::Fourier_directions;

// for reading Fourier space data (both full and reduced)

H5_IO BasicIO::full;
H5_IO BasicIO::kz0_full;
H5_IO BasicIO::reduced;
H5_IO BasicIO::kz0_reduced;
H5_IO BasicIO::real;

	//private:
hid_t  BasicIO::err=0;

map<string, string> BasicIO::transition_table;

hid_t BasicIO::dataspace=0;
hid_t BasicIO::memspace=0;
hid_t BasicIO::dataset=0;
hid_t BasicIO::acc_template=0;

// define an info object to store MPI-IO information 
MPI_Info BasicIO::FILE_INFO_TEMPLATE;
hid_t BasicIO::file_identifier;
int BasicIO::ierr;

int BasicIO::rank=3;
int BasicIO::dim1;
int BasicIO::dim2;
int BasicIO::dim3;
DP *BasicIO::dataArray;

string BasicIO::folder_name="";
string BasicIO::file_name="";

bool BasicIO::frequent=false;

string BasicIO::LOG_filename="README.log";
ofstream BasicIO::LOG_file;
time_t BasicIO::rawtime;

//*********************************************************************************************

void BasicIO::Initialize()
{
	// set the file access template for parallel IO access 
	acc_template = H5Pcreate(H5P_FILE_ACCESS);

	// create an MPI_INFO object -- on some platforms it is useful to
	// pass some information onto the underlying MPI_File_open call

	ierr = MPI_Info_create(&FILE_INFO_TEMPLATE);
	ierr = H5Pset_sieve_buf_size(acc_template, SIEVE_BUF_SIZE);
	ierr = H5Pset_alignment(acc_template, 524288, 262144);

	//
	// visit http://www.mpi-forum.org/docs/mpi21-report/node267.htm for information on access_style, ..

	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"access_style", (char*)"write_once");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"collective_buffering", (char*)"true");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_block_size",  (char*)"1048576");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_buffer_size", (char*)"4194304");

	// tell the HDF5 library that we want to use MPI-IO to do the writing 
	ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, FILE_INFO_TEMPLATE);

	dim1=0;
	dim2=0;
	dim3=0;
	
	//If data_out_folder does not exist then create it
	string out_path = global.io.data_dir + data_out_folder;
	mkdir(out_path.c_str(), S_IRWXU|S_IRWXG);

	//Open the log file
	if (global.mpi.master) LOG_file.open((global.io.data_dir + data_out_folder + "/" + LOG_filename).c_str());

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

	Set_parameters();
}


void BasicIO::Finalize()
{
	ierr = H5Pclose(acc_template);
	ierr = MPI_Info_free(&FILE_INFO_TEMPLATE);
	if (global.mpi.master) LOG_file.close();
	
}

//*********************************************************************************************

void BasicIO::Log(string message){
	time ( &rawtime );
 	string t_now = ctime(&rawtime);
	t_now = t_now.substr( 0, t_now.find_last_of('\n'));
	if (global.mpi.master) LOG_file << "[ " << t_now << " ] " << message << endl;
}

void BasicIO::Begin_frequent(){
	frequent=true;
	Log("Field Frequent writing started.  Time = "+To_string(global.time.now));
}
void BasicIO::End_frequent(){
	frequent=false;
	Log("Field Frequent writing finisheded.");
}

//*********************************************************************************************

vector<H5_dataset_meta> BasicIO::Get_meta(string filename){
	vector<H5_dataset_meta> file_meta_data;
	H5_dataset_meta dataset_meta_data;
	int max_rank=3;
	int max_object_name_length = 3;
	
	herr_t   status;
	hsize_t num_obj;
	int object_type;
	hid_t group_id, type_id, dataset_id;
	char object_name[max_object_name_length+1];	// +1 for null character
	hid_t           space;    /* Handles */
	hsize_t         dims[max_rank],
	rank;
	
	/*
	 *  Open a file, open the root and scan the root.
	 */
	file_identifier = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	
	group_id = H5Gopen(file_identifier,"/", H5P_DEFAULT);
	err = H5Gget_num_objs(group_id, &num_obj);
	
	for (int i = 0; i < num_obj; i++) {
		H5Gget_objname_by_idx(group_id, (hsize_t)i, object_name, (size_t)(max_object_name_length+1));
		object_type =  H5Gget_objtype_by_idx(group_id, (size_t)i );
		
		if (object_type == H5G_DATASET){  //Member is a dataset
			dataset_meta_data.name = object_name;
			
			dataset_id = H5Dopen(group_id, object_name, H5P_DEFAULT);
			space = H5Dget_space(dataset_id);
			
			rank = H5Sget_simple_extent_dims (space, dims, NULL);
			if (rank==3)
				dataset_meta_data.dimensions.assign(dims, dims + 3);
			else
				cerr << "Not a tarang-1 output field file." << endl;
			
			type_id = H5Dget_type(dataset_id);
			dataset_meta_data.datatype = H5Tget_class(type_id);
			
			file_meta_data.push_back(dataset_meta_data);
			
			H5Dclose(dataset_id);
		}
		else{
			cerr << "Not a tarang-1 output field file." << endl;
			break;
		}
	}
	
	status = H5Fclose(file_identifier);
	
	return file_meta_data;
}


//******************************************************************************


void BasicIO::Read(Array<complx,3> A, string dataset_name, H5_IO config)
{
	Read_basic(reinterpret_cast<DP*>(A.data()), dataset_name, config);
	
	shape_complex_array = global.field.shape_complex_array;
	
	//cout << "shape complex array: " << shape_complex_array << endl;
	
	if (config == reduced || config == kz0_reduced) {
		
		cout << "basic IO set 0" << endl;
		TinyVector<int, 3> direction_z_complex_array_complement(3);
		direction_z_complex_array_complement = (direction_z_complex_array^1);
		
		if (Fourier_directions(0)*direction_z_complex_array_complement(0)){		//Fourier and not Z
			Array<int,1> filter_y(shape_full_complex_array(0));
			Array<int,1> my_filter_y(shape_complex_array(0));
			
			filter_y=0;
			
			filter_y(Ny-shape_in_reduced_array(0)/2)=1;
			my_filter_y=filter_y(Range(id_complex_array(0)*shape_complex_array(0), (id_complex_array(0)+1)*shape_complex_array(0)-1));
			
			//cout << "filters: " << filter_y << " " << my_filter_y << endl;
			int index = first(my_filter_y==1);
			//cout << "myid, index: " << my_id << " " << id_complex_array(0) << " " << index << endl;
			if ( (index >= 0) && (index < A.extent(0)) )
			   A(index, Range::all(), Range::all()) = 0;
			
			cout << "Y: my_id, index " << my_id << " " << index << endl;
		}

		if (Fourier_directions(2)*direction_z_complex_array_complement(2)){		//Fourier and not Z
			Array<int,1> filter_x(shape_full_complex_array(2));
			Array<int,1> my_filter_x(shape_complex_array(2));
			
			filter_x=0;
			
			filter_x(Nx-shape_in_reduced_array(2)/4)=1;
			my_filter_x=filter_x(Range(id_complex_array(2)*shape_complex_array(2), (id_complex_array(2)+1)*shape_complex_array(2)-1));
			
			int index = first(my_filter_x==1);
			if ( (index >= 0) && (index < A.extent(2)) )
				A(Range::all(), Range::all(), index) = 0;
			
			cout << "X: my_id, index " << my_id << " " << index << endl;
		}
	}
}

void BasicIO::Read(Array<DP,3> A, string dataset_name, H5_IO config)
{
	Read_basic(A.data(), dataset_name, config);
}

//******************************************************************************

void BasicIO::Write(Array<complx,3> A, string dataset_name, H5_IO config)
{
	Write_basic(reinterpret_cast<DP*>(A.data()), dataset_name, config);
}

void BasicIO::Write(Array<DP,3> A, string dataset_name, H5_IO config)
{
	Write_basic(A.data(), dataset_name, config);
}


//******************************************************************************
void BasicIO::Read_basic(DP *data, string dataset_name, H5_IO config)
{	

	/*for (int i=0; i<numprocs; i++){
		if (my_id == i)
			cerr << "Read_basic: my_id = " << my_id << " " << dataset_name << "\n    in:\n" << config.in << endl;
		usleep(100);
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
    
	dataspace = H5Screate_simple(rank, config.in.dataspace.dimension.data(), NULL);
    

	err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, config.in.dataspace.start.data(), config.in.dataspace.blockstride.data(), config.in.dataspace.blockcount.data(), config.in.dataspace.blockdim.data());
    
	MPI_Barrier(MPI_COMM_WORLD);
	
	memspace = H5Screate_simple(rank, config.in.memspace.dimension.data(), NULL);
    
	err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, config.in.memspace.start.data(), config.in.memspace.blockstride.data(), config.in.memspace.blockcount.data(), config.in.memspace.blockdim.data());
    
	MPI_Barrier(MPI_COMM_WORLD);
	
	folder_name = global.io.data_dir + "/" + data_in_folder;
	file_name=folder_name + "/" + dataset_name+".h5";
	file_identifier = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, acc_template);
    
    if (file_identifier<0){
    	cerr << "BasicIO:: Unable to open " << file_name << " for reading." << endl;
    	return;
    }
    
	// now create the dataset -- if it is the first time through, otherwise
	//     open the exisiting dataset 
	dataset = H5Dopen2(file_identifier, dataset_name.c_str(), H5P_DEFAULT);

    
	// and finally read data from disk -- in this call, we need to
	//     include both the memory space and the data space.
	err = H5Dread(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, data);
	
	MPI_Barrier(MPI_COMM_WORLD);
    
	H5Dclose(dataset);
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Fclose(file_identifier);
	
}

//*********************************************************************************************

void BasicIO::Write_basic(DP* data, string dataset_name, H5_IO config)
{
	/*for (int i=0; i<numprocs; i++){
		if (my_id == i)
			cerr << "Write_basic: my_id = " << my_id << " " << dataset_name << "\n    out:\n" << config.out << endl;
		usleep(100);
		MPI_Barrier(MPI_COMM_WORLD);
	}*/

	dataspace = H5Screate_simple(rank, config.out.dataspace.dimension.data(), NULL);

    
	// figure out the offset into the dataspace for the current processor 
    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, config.out.dataspace.start.data(), config.out.dataspace.blockstride.data(), config.out.dataspace.blockcount.data(), config.out.dataspace.blockdim.data());
    
	memspace = H5Screate_simple(rank, config.out.memspace.dimension.data(), NULL);
    
	err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, config.out.memspace.start.data(), config.out.memspace.blockstride.data(), config.out.memspace.blockcount.data(), config.out.memspace.blockdim.data());
 
    
    //Create the time.now folder in data_dir/out
	folder_name = global.io.data_dir + data_out_folder; 

	if (frequent)
		folder_name += "/frequent";

	else if (config==reduced || config==kz0_reduced)
		folder_name += "/reduced_time_" + To_string(global.time.now);
	
	else if (config==real)
		folder_name += "/real_time_" + To_string(global.time.now);

	else
		folder_name += "/time_" + To_string(global.time.now);
	
	
	mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG);
    
    
	//Open the file for writing
	file_name= folder_name + "/" + dataset_name + ".h5";
	file_identifier = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    if (file_identifier<0){
    	cerr << "BasicIO:: Unable to open " << file_name << " for writing." << endl;
    	return;
    }

	// now create the dataset -- if it is the first time through, otherwise
	//     open the exisiting dataset
	dataset = H5Dcreate2(file_identifier, dataset_name.c_str(), H5_DATATYPE,
                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
	// and finally write the data to disk -- in this call, we need to
	//     include both the memory space and the data space.
	err = H5Dwrite(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, data);
    
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	H5Fclose(file_identifier);
}

template<class T, int N_rank>
bool operator==(TinyVector<T,N_rank> A, TinyVector<T,N_rank> B){
	for (int i=0; i<N_rank; i++)
		if (A(i)!=B(i))
			return false;
	return true;
}

ostream& operator<< (ostream &out, H5_space &h5_space)
{
	out << "            dimension = " << h5_space.dimension << endl;
	out << "            blockdim = " << h5_space.blockdim << endl;
	out << "            start = " << h5_space.start << endl;
	out << "            blockstride = " << h5_space.blockstride << endl;
	out << "            blockcount = " << h5_space.blockcount << endl;

   return out;
}

bool operator== (H5_space &h5_space1, H5_space &h5_space2)
{
	return  (h5_space1.dimension == h5_space2.dimension) &&
	        (h5_space1.blockdim == h5_space2.blockdim) &&
	       (h5_space1.start == h5_space2.start) &&
	       (h5_space1.blockstride == h5_space2.blockstride) &&
	       (h5_space1.blockcount == h5_space2.blockcount);
}


ostream& operator<< (ostream &out, H5_space_group &H5_space_group)
{
	out << "        dataspace:\n" << H5_space_group.dataspace << endl;
	out << "        memspace:\n" << H5_space_group.memspace << endl;
	
   return out;
}

bool operator== (H5_space_group &h5_space_group1, H5_space_group &h5_space_group2){
	return (h5_space_group1.dataspace == h5_space_group2.dataspace);
		   (h5_space_group1.memspace == h5_space_group2.memspace);
}


ostream& operator<< (ostream &out, H5_IO &h5_io)
{
   out << "    in:\n" << h5_io.in << endl;
   out << "    out:\n" << h5_io.in << endl;

   return out;
}

bool operator== (H5_IO &h5_io1, H5_IO &h5_io2){
	return (h5_io1.in == h5_io2.in) &&
	      (h5_io1.out == h5_io2.out);
}

ostream& operator<< (ostream &out, Filter_prop &prop)
{
   out << "    first_rise: " << prop.first_rise << endl;
   out << "    first_fall: " << prop.first_fall << endl;
   out << "    zero_length: " << prop.zero_length << endl;

   return out;
}

//========================= Class declaration of BasicIO ends ============================== 

