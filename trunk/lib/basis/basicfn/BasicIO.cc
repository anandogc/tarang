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

//public:
string BasicIO::data_in_folder="";
string BasicIO::data_out_folder="";

MPI_Comm BasicIO::MPI_COMMUNICATOR;


//private:

bool BasicIO::frequent=false;

string BasicIO::LOG_filename="README.log";
ofstream BasicIO::LOG_file;

//*********************************************************************************************

void BasicIO::Initialize(string data_base_folder, MPI_Comm MPI_COMMUNICATOR)
{
	data_in_folder = data_base_folder + "/in";
	data_out_folder = data_base_folder + "/out";

	//If data_out_folder does not exist then create it
	mkdir(data_out_folder.c_str(), S_IRWXU|S_IRWXG);

	//Open the log file
	if (master && LOG_filename.size()>0)
		LOG_file.open((data_out_folder + "/" + LOG_filename).c_str());

	BasicIO::MPI_COMMUNICATOR = MPI_COMMUNICATOR;

}


void BasicIO::Finalize() {
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
	Log("Field Frequent writing finished.");
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


//Get list of datasets and their propertties.
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

//Used in convertors
//filters "*.h5" in a given directory
int BasicIO::h5filter(const struct dirent *dir)
// post: returns 1/true if name of dir ends in .h5
{
	const char *s = dir->d_name;
	int len = strlen(s) - 3;	// index of start of . in .h5
	if ( (len >= 0) && (strncmp(s + len, ".h5", 3) == 0) )
		return 1;
	return 0;
}

//Used in convertors
//Get list of files in a given directory using 'h5filter'
vector<string> BasicIO::Get_file_list(string path) {
	vector<string> name_list;

	struct dirent **namelist=NULL;
	int n=-1;

	n = scandir(path.c_str(), &namelist, h5filter, alphasort);
	if (n < 0)
		if (master)
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

template<typename T>
ostream& operator<<(ostream& output, const std::vector<T>& v) {
	output << "[";
	for (typename std::vector<T>::size_type i=0; i<v.size()-1; i++)
		output << v[i] << ", ";
	if (v.size()>0)
		output << v[v.size()-1];
	output << "]";

	return output;
}


//******************************************************************************
int BasicIO::Read(void *data, h5::Plan plan, string file_name, string dataset_name)
{

	string file_path = data_in_folder+"/"+file_name+".h5";

	if (!File_exists(file_path)){
		if (master)
			cerr << "File does not exist: '" << file_path << "'" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
	}


	if (dataset_name.length()==0)
		dataset_name = file_name;

	h5::File f(file_path, "r");
	f.mpioInit(MPI_COMMUNICATOR);

	std::vector<std::string> keys = f.keys();

	if (!f.contains(dataset_name) || !f[dataset_name].isDataset()) {
		if (master) {
			cerr << "Dataset '" << dataset_name << "' not found in file '" << file_path << "'" << endl;
			MPI_Abort(MPI_COMMUNICATOR, 1);
		}
	}

	h5::Dataset ds = f[dataset_name];
	//Match dimensions of dataset with dimensions in plan

	std::vector<hsize_t> dataset_dimensions = ds.shape();
	std::vector<hsize_t> planned_dataset_dimensions = plan.filespace_dimension();

	if ( !( dataset_dimensions.size()==planned_dataset_dimensions.size() and std::equal(dataset_dimensions.begin(), dataset_dimensions.end(), planned_dataset_dimensions.begin()) ) ){
		if (master) {

			cerr << "Dataset dimensions in input file '" << file_name << "' [";

			if (dataset_dimensions.size()>0)
				cerr << dataset_dimensions[0];

			for (std::vector<hsize_t>::size_type i=1; i<dataset_dimensions.size(); i++)
				cerr << "x" << dataset_dimensions[i];

			cerr << "] differ from specified dimensions [";

			if (planned_dataset_dimensions.size()>0)
				cerr << planned_dataset_dimensions[0];
			
			for (std::vector<hsize_t>::size_type i=1; i<planned_dataset_dimensions.size(); i++)
				cerr << "x" << planned_dataset_dimensions[i];

			cerr << "]" << endl;

			MPI_Abort(MPI_COMMUNICATOR, 1);
		}
	}

	ds.set_plan(plan);
	ds >> data;
	
	return 0;
}

//*********************************************************************************************

int BasicIO::Write(const void* data, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
{
	if (frequent)
		folder_name = "frequent";

	//Create the time.now folder in data_dir/out
	string folder_path = data_out_folder + "/" + folder_name; 
	mkdir(folder_path.c_str(), S_IRWXU | S_IRWXG);
	
	//Open the file for writing
	string file_path = folder_path + "/" + file_name + ".h5";

	if (dataset_name.length()==0)
		dataset_name = file_name;

	h5::File f;
	f.mpioInit(MPI_COMMUNICATOR);
	f.open(file_path, "w");

	h5::Dataset ds = f.create_dataset(dataset_name, plan);
	ds << (void*)data;

	return 0;
}

//========================= Class declaration of BasicIO ends ============================== 

