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
#include "Global_extern_vars.h"
#include <universal.h>
//*********************************************************************************************

hid_t BasicIO::H5T_FLOAT=-1;
hid_t BasicIO::H5T_DOUBLE=-1;

hid_t BasicIO::H5T_COMPLEX_FLOAT=-1;
hid_t BasicIO::H5T_COMPLEX_DOUBLE=-1;

//public:
string BasicIO::data_in_folder="";
string BasicIO::data_out_folder="";


//private:

hid_t BasicIO::acc_template=0;

// define an info object to store MPI-IO information 
MPI_Info BasicIO::FILE_INFO_TEMPLATE;

bool BasicIO::frequent=false;

string BasicIO::LOG_filename="README.log";
ofstream BasicIO::LOG_file;

//*********************************************************************************************

void BasicIO::Initialize(string data_base_folder)
{
	data_in_folder = data_base_folder + "/in";
	data_out_folder = data_base_folder + "/out";

	// set the file access template for parallel IO access 
	acc_template = H5Pcreate(H5P_FILE_ACCESS);

	// create an MPI_INFO object -- on some platforms it is useful to
	// pass some information onto the underlying MPI_File_open call

	MPI_Info_create(&FILE_INFO_TEMPLATE);
	H5Pset_sieve_buf_size(acc_template, SIEVE_BUF_SIZE);
	H5Pset_alignment(acc_template, 524288, 262144);

	//
	// visit http://www.mpi-forum.org/docs/mpi21-report/node267.htm for information on access_style, ..

	MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"access_style", (char*)"write_once");
	MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"collective_buffering", (char*)"true");
	MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_block_size",  (char*)"1048576");
	MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_buffer_size", (char*)"4194304");

	// tell the HDF5 library that we want to use MPI-IO to do the writing 
	H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, FILE_INFO_TEMPLATE);
	
	//If data_out_folder does not exist then create it
	mkdir(data_out_folder.c_str(), S_IRWXU|S_IRWXG);

	//Open the log file
	if (master && LOG_filename.size()>0)
		LOG_file.open((data_out_folder + "/" + LOG_filename).c_str());


	H5T_FLOAT = H5Tcopy(H5T_NATIVE_FLOAT);
	H5T_DOUBLE = H5Tcopy(H5T_NATIVE_DOUBLE);

	H5T_COMPLEX_FLOAT = H5Tcreate(H5T_COMPOUND, sizeof(complex<float>));
	H5Tinsert(H5T_COMPLEX_FLOAT, "real", 0, H5T_NATIVE_FLOAT);
	H5Tinsert(H5T_COMPLEX_FLOAT, "imag", sizeof(float), H5T_NATIVE_FLOAT);

	H5T_COMPLEX_DOUBLE = H5Tcreate(H5T_COMPOUND, sizeof(complex<double>));
	H5Tinsert(H5T_COMPLEX_DOUBLE, "real", 0, H5T_NATIVE_DOUBLE);
	H5Tinsert(H5T_COMPLEX_DOUBLE, "imag", sizeof(double), H5T_NATIVE_DOUBLE);
}


void BasicIO::Finalize()
{
	H5Pclose(acc_template);
	MPI_Info_free(&FILE_INFO_TEMPLATE);
	if (master) LOG_file.close();
}

//*********************************************************************************************

void BasicIO::Log(string message) {
	time_t rawtime;
	time ( &rawtime );
 	string t_now = ctime(&rawtime);
	t_now = t_now.substr( 0, t_now.find_last_of('\n'));
	if (master) LOG_file << "[ " << t_now << " ] " << message << endl;
}

void BasicIO::Begin_frequent() {
	frequent=true;
	Log("Field Frequent writing started.  Time = "+To_string(global.time.now));
}
void BasicIO::End_frequent() {
	frequent=false;
	Log("Field Frequent writing finisheded.");
}

//*********************************************************************************************

bool BasicIO::File_exists(const std::string& name)
{
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

vector<BasicIO::H5_dataset_meta> BasicIO::Get_meta(string filename, string path)
{
	if (!File_exists(filename)) {
		if (master)
			cerr << filename << " does not exists." << endl;
		exit(1);
	}

	vector<H5_dataset_meta> file_meta_data;
	H5_dataset_meta dataset_meta_data;
	int max_rank=3;
	int max_object_name_length = 100;
	
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
	hid_t file_identifier = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	
	group_id = H5Gopen(file_identifier, path.c_str(), H5P_DEFAULT);
	H5Gget_num_objs(group_id, &num_obj);
	
	for (int i = 0; i < num_obj; i++) {
		H5Gget_objname_by_idx(group_id, (hsize_t)i, object_name, (size_t)(max_object_name_length+1));
		object_type =  H5Gget_objtype_by_idx(group_id, (size_t)i );
		
		dataset_meta_data.name = object_name;

		dataset_id = H5Dopen(group_id, object_name, H5P_DEFAULT);
		space = H5Dget_space(dataset_id);
		
		rank = H5Sget_simple_extent_dims (space, dims, NULL);
		dataset_meta_data.dimensions.assign(dims, dims + 3);
		
		type_id = H5Dget_type(dataset_id);
		dataset_meta_data.datatype = H5Tget_class(type_id);
		
		file_meta_data.push_back(dataset_meta_data);
		
		H5Dclose(dataset_id);
	}
	
	status = H5Fclose(file_identifier);
	
	return file_meta_data;
}

int BasicIO::h5filter(const struct dirent *dir)
// post: returns 1/true if name of dir ends in .h5
{
	const char *s = dir->d_name;
	int len = strlen(s) - 3;	// index of start of . in .h5
	if ( (len >= 0) && (strncmp(s + len, ".h5", 3) == 0) )
		return 1;
	return 0;
}

vector<string> BasicIO::Get_file_list(string path){
	vector<string> name_list;

	struct dirent **namelist=NULL;
	int n=-1;

	n = scandir(path.c_str(), &namelist, h5filter, alphasort);
	if (n < 0)
	    cerr << "No such directory: " << path << endl;
	else {
	    for (int i = 0; i < n; i++) {
	        name_list.push_back(namelist[i]->d_name);
	        free(namelist[i]);
	        }
    }
	free(namelist);

	return name_list;
}


//******************************************************************************

herr_t Hdf5_reader_error_handler(int index,void *file_path)
{
	if (master){
		cerr << "Error reading '" << (char*)file_path << "'" << endl;
		H5Eprint(H5E_DEFAULT, stderr);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit (1);
}

//******************************************************************************
int BasicIO::Read(void *data, H5_plan plan, string file_name, string dataset_name)
{	
	string file_path = data_in_folder+"/"+file_name+".h5";

	if (!File_exists(file_path)){
		if (master)
			cerr << "File does not exist: '" << file_path << "'" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
	}

	H5Eset_auto(H5E_DEFAULT, Hdf5_reader_error_handler, (void*)file_path.c_str());
	hid_t file_identifier = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, acc_template);

    if (file_identifier<0){
    	if (master) cerr << "BasicIO:: Unable to open " << file_name << " for reading." << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
    }

    if (dataset_name.length()==0)
    	dataset_name = file_name;

	// now create the dataset -- if it is the first time through, otherwise
	//     open the exisiting dataset 
	hid_t dataset = H5Dopen2(file_identifier, dataset_name.c_str(), H5P_DEFAULT);

    
	// and finally read data from disk -- in this call, we need to
	//     include both the memory space and the data space.
	herr_t err = H5Dread(dataset, plan.datatype, plan.memspace, plan.dataspace, H5P_DEFAULT, data);
	
	H5Dclose(dataset);
	H5Fclose(file_identifier);

    if (err<0){
    	if (master)
    		cerr << "BasicIO:: Error in reading from file  '" << file_name << "'." << endl;
    	return 1;
    }

	return 0;
}

//*********************************************************************************************

herr_t Hdf5_writer_error_handler(int index,void *file_path)
{
	if (master){
		cerr << "Error writing '" << (char*)file_path << "'" << endl;
		H5Eprint(H5E_DEFAULT, stderr);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit (1);
}

int BasicIO::Write(const void* data, H5_plan plan, string folder_name, string file_name, string dataset_name)
{
   

   	if (frequent)
		folder_name = "frequent";

    //Create the time.now folder in data_dir/out
	string folder_path = data_out_folder + "/" + folder_name; 
	mkdir(folder_path.c_str(), S_IRWXU | S_IRWXG);
    
	
	//Open the file for writing
	string file_path = folder_path + "/" + file_name + ".h5";

	H5Eset_auto(H5E_DEFAULT, Hdf5_writer_error_handler, (void*)file_path.c_str());
	hid_t file_identifier = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);

    if (file_identifier<0){
    	cerr << "BasicIO:: Unable to open " << file_name << " for writing." << endl;
    	return 1;
    }


    if (dataset_name.length()==0)
    	dataset_name = file_name;

	// now create the dataset -- if it is the first time through, otherwise
	//     open the exisiting dataset

	hid_t dataset = H5Dcreate2(file_identifier, dataset_name.c_str(), plan.datatype,
                         plan.dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
	// and finally write the data to disk -- in this call, we need to
	//     include both the memory space and the data space.
	herr_t err = H5Dwrite(dataset, plan.datatype, plan.memspace, plan.dataspace, H5P_DEFAULT, data);
    
	H5Dclose(dataset);
	H5Fclose(file_identifier);

    if (err<0){
    	if (master)
    		cerr << "BasicIO:: Error in writing to file  '" << file_name << "'." << endl;
    	return 1;
    }

	return 0;
}


/*************
* Structures and Functions useful for:
* BasicIO::H5_plan BasicIO::Set_plan(int rank, int* my_id, int* numprocs, Array<int,1>* dataspace_filter, Array<int,1>* memspace_filter, hid_t datatype)
*
* These structures and functions are not required anywhere else.
*/

struct H5_filter_blocks
{
	size_t dataspace_start;
	size_t memspace_start;
	size_t blocklength;

	H5_filter_blocks(size_t dataspace_start, size_t memspace_start, size_t blocklength): dataspace_start(dataspace_start), memspace_start(memspace_start), blocklength(blocklength){}
};

 struct H5_hyperspace{
	Array<hsize_t, 1> dimension;
	Array<hsize_t, 1> start;
	Array<hsize_t, 1> blockdim;
	Array<hsize_t, 1> blockstride;
	Array<hsize_t, 1> blockcount;

	H5_hyperspace(size_t rank);
};

ostream &operator<<(ostream &out, H5_hyperspace hyperspace){
	out << "dimension: " << hyperspace.dimension << '\n'
	    << "start: " << hyperspace.start << '\n'
	    << "blockdim: " << hyperspace.blockdim << '\n'
	    << "blockstride: " << hyperspace.blockstride << '\n'
	    << "blockcount: " << hyperspace.blockcount << endl;

	return out;
}

H5_hyperspace::H5_hyperspace(size_t rank){
	dimension.resize(rank);
	start.resize(rank);
	blockdim.resize(rank);
	blockstride.resize(rank);
	blockcount.resize(rank);

	dimension=0;
	start=0;
	blockdim=0;
	blockstride=1;
	blockcount=0;
}


hsize_t Get_dataspace_start(Array<int,1> filter, hsize_t count){
	size_t sum=0;
	hsize_t index=0;
	for (index=0; ( (index<filter.size()) && (sum<count) ); index++)
		sum += filter(index);

	return index-1;
}

void Get_properties_memspace(Array<int,1> filter, hsize_t begin_at, int my_id, int numprocs, hsize_t &start, hsize_t &blocklength)
{
	hsize_t local_start = my_id*filter.size()/numprocs;
	hsize_t local_end = (my_id+1)*filter.size()/numprocs-1;
	hsize_t end;

	start = local_start + begin_at + first(filter(Range(local_start+begin_at, local_end))==1);
	end = start + first(filter(Range(start, local_end))==0);

	if (start<local_start || start >local_end){
		start=0;
		blocklength=0;
	}
	else if (end<local_start || end >local_end){
		end=local_end+1;
		blocklength=end-start;
		start-=local_start;
	}
	else{
		blocklength=end-start;
		start-=local_start;
	}
}



BasicIO::H5_plan BasicIO::Set_plan(int rank, int* my_id, int* numprocs, Array<int,1>* dataspace_filter, Array<int,1>* memspace_filter, hid_t datatype)
{	
	vector<H5_filter_blocks>* filter_blocks;
	filter_blocks = new vector<H5_filter_blocks>[rank];

	size_t local_start_index;
	//size_t local_end;

	hsize_t dataspace_start, memspace_start, blocklength;
	
	hsize_t dataspace_dimension[rank];
	hsize_t memspace_dimension[rank];

	for (int r=0; r<rank; r++){
		if (sum(dataspace_filter[r]) != sum(memspace_filter[r])){
			if (master)
				cerr << "BasicIO::Set_plan: dataspace_filter["<<r<<"] and memspace_filter["<<r<<"] has different counts." << endl;
			exit(1);
		}

		dataspace_dimension[r]=dataspace_filter[r].size();
		memspace_dimension[r]=memspace_filter[r].size()/numprocs[r];

		local_start_index = my_id[r]*memspace_filter[r].size()/numprocs[r];
		//local_end = (my_id[r]+1)*memspace_filter[r].size()/numprocs[r]-1; 

		int j=0;
		do {
			Get_properties_memspace(memspace_filter[r], j, my_id[r], numprocs[r], memspace_start, blocklength);
			j = memspace_start + blocklength;

			if (j>0){
				dataspace_start = Get_dataspace_start(dataspace_filter[r],sum(memspace_filter[r](Range(0, local_start_index + memspace_start))));

				filter_blocks[r].push_back(H5_filter_blocks(dataspace_start, memspace_start, blocklength));
			}
		}
		while(j>0);
	}


	/*for (int r=0; r<rank; r++){
		for (int j=0; j<filter_blocks[r].size(); j++){
			cout << "j = " << my_id[r] << " " << j << " " << filter_blocks[r][j].dataspace_start << " " << filter_blocks[r][j].memspace_start << " " << filter_blocks[r][j].blocklength << endl;
		}
		cout << endl;
	}*/

	size_t num_combinations=1;
	for (int r=0; r<rank; r++)
		num_combinations*=filter_blocks[r].size();

	vector<H5_hyperspace> dataspace;
	vector<H5_hyperspace> memspace;

	dataspace.reserve(num_combinations);
	memspace.reserve(num_combinations);

	for (size_t i=0; i<num_combinations; i++){
		dataspace.push_back(H5_hyperspace(rank));
		memspace.push_back(H5_hyperspace(rank));
	}

	H5_plan h5_plan;

	h5_plan.dataspace = H5Screate_simple(rank, dataspace_dimension, NULL);
	h5_plan.memspace = H5Screate_simple(rank, memspace_dimension, NULL);

	//Select no elements in the hyperspace
	H5Sselect_none(h5_plan.dataspace);
	H5Sselect_none(h5_plan.memspace);

	if (num_combinations>0){
		//Perform all the combination
		size_t repeat_rth_rank=num_combinations;
		size_t index;

		for (int r=0; r<rank; r++){
			index=-1;
			repeat_rth_rank = repeat_rth_rank/filter_blocks[r].size();

			for (size_t j=0; j<num_combinations;){
				index = (index+1)%filter_blocks[r].size();

				for (size_t k=0; k<repeat_rth_rank; j++, k++){
					dataspace[j].dimension(r)=dataspace_filter[r].size();
					dataspace[j].start(r)=filter_blocks[r][index].dataspace_start;
					dataspace[j].blockdim(r)=filter_blocks[r][index].blocklength;
					dataspace[j].blockstride(r)=1;
					dataspace[j].blockcount(r)=1;

					memspace[j].dimension(r)=memspace_filter[r].size()/numprocs[r];
					memspace[j].start(r)=filter_blocks[r][index].memspace_start;
					memspace[j].blockdim(r)=filter_blocks[r][index].blocklength;
					memspace[j].blockstride(r)=1;
					memspace[j].blockcount(r)=1;
				}
			}
		}

	    
		//Select required region in the hyperspace
	    for (int i=0; i<num_combinations; i++){
	    	/*if (my_id[2]==0){
				cout << "Dataspace: " << my_id[2] << " " << dataspace[i].dimension << " " << dataspace[i].start << " " << dataspace[i].blockstride << " " << dataspace[i].blockcount << " " << dataspace[i].blockdim << endl;

				cout << "Memspace: " << my_id[2] << " " << memspace[i].dimension << " " << memspace[i].start << " " << memspace[i].blockstride << " " << memspace[i].blockcount << " " << memspace[i].blockdim << endl;
			}*/
			

	    	 H5Sselect_hyperslab(h5_plan.dataspace, H5S_SELECT_OR, dataspace[i].start.data(), dataspace[i].blockstride.data(), dataspace[i].blockcount.data(), dataspace[i].blockdim.data());

	    	 H5Sselect_hyperslab(h5_plan.memspace, H5S_SELECT_OR, memspace[i].start.data(), memspace[i].blockstride.data(), memspace[i].blockcount.data(), memspace[i].blockdim.data());
	    }
	}

	h5_plan.datatype = datatype;
	return h5_plan;
}
//========================= Class declaration of BasicIO ends ============================== 

