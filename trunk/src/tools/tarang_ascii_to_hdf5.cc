/* Tarang-2
 *
 * Copyright (C) 2008, 2009 Mahendra K. Verma
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

/*! \file tarang_ascii_to_hdf5.cc
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author A. G. Chatterjee, Abhishek Kumar, M. K. Verma
 * @date 2 April 2013
 * 
 * @bugs  No known bug
 */ 
 
#include "BasicIO.h"
#include "BasicIO_SP.h"

Global global;

int Convert_PRINFTY();
int Convert_FLUID();
int Convert_RBC();
int Convert_MHD();
int Convert_MRBC();

string file_name;
bool is_plane;
	
ifstream file_in;


struct IO_plans{
    BasicIO::H5_plan H5_real;
    BasicIO::H5_plan H5_full;
    BasicIO::H5_plan H5_kz0_full;
    BasicIO::H5_plan H5_in_reduced;
    BasicIO::H5_plan H5_out_reduced;
    BasicIO::H5_plan H5_in_kz0_reduced;
    BasicIO::H5_plan H5_out_kz0_reduced;
};

IO_plans tarang_2_3;

int main(int argc, char** argv)
{	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &global.mpi.my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &global.mpi.numprocs);
	
	if (global.mpi.numprocs!=1){
		if (global.mpi.my_id == 0)
			cerr << "This program is made to run on one processor!" << endl;
		exit(1);
	}

	global.mpi.master = (global.mpi.my_id == 0);
	
	if (argc < 6){
        cerr<< "Usage:\n    " << argv[0] << " /path/to/field_in.d PRINFTY/FLUID/RBC/MHD/MRBC dim1 dim2 dim3 [PLANE]" << endl;
        return 1;
	}

	global.Parse(argc, argv, false);
	BasicIO::Initialize();

	global.io.input_vx_vy_switch = (argc==7);

    
	global.io.output_vx_vy_switch=global.io.input_vx_vy_switch;

	global.field.N[1] = atoi(argv[3]);
	global.field.N[2] = atoi(argv[4]);
	global.field.N[3] = atoi(argv[5]);

	global.Process_basic_vars();

	global.program.kind = argv[2];

	if (Nx*Ny*Nz <= 0){
		cerr << "Invalid grid size." << endl;
		exit(1);
	}

	if (master) {
		cout << "N = [" << Nx << "," << Ny << "," << Nz << "]" << endl;
		cout << boolalpha << "input_vx_vy_switch = " << global.io.input_vx_vy_switch << endl;
		cout << boolalpha << "output_vx_vy_switch = " << global.io.output_vx_vy_switch << endl;
		
		cout << endl;
	}


	file_name = global.io.data_dir + "/tarang1/field_in.d";
	file_in.open(file_name.c_str());


	//Configure tarang 2.3 output
	
	BasicIO::Array_properties<3> array_properties;

	array_properties.shape_full_complex_array = Ny, Nz/2+1, Nx;
	array_properties.shape_full_real_array = Ny, Nz+2, Nx;

	array_properties.id_complex_array = my_id, 0, 0;
	array_properties.id_real_array = 0, 0, my_id;

	array_properties.numprocs_complex_array = numprocs, 1, 1;
	array_properties.numprocs_real_array = 1, 1, numprocs;

	array_properties.Fourier_directions = 1,1,1;
	array_properties.Z = 1;

	array_properties.datatype_complex_space = BasicIO::H5T_Complex;
	array_properties.datatype_real_space = BasicIO::H5T_Real;

	BasicIO::Set_H5_plans(array_properties, &tarang_2_3);


	int status;

	if (global.program.kind == "PRINFTY")
			status = Convert_PRINFTY();
	/*else if (simulation_type == "FLUID")
			status = Convert_FLUID();
	else if (simulation_type == "RBC")
			status = Convert_RBC();
	else if (simulation_type == "MHD")
			status = Convert_MHD();
	else if (simulation_type == "MRBC")
			status = Convert_MRBC();*/
	else{
			cerr << "Not a valid Simulation type." << endl;
			return(1);
		}

	BasicIO::Finalize();
	MPI_Finalize();
	
	return 0;
}


void Read_complex_array(Array<Complex, 3> &A, int Nz, string field_name)
{	
	cout << "Reading " << field_name << " ..." << flush;
	Complex temp_x;
	for (int i=0; i<Nx; i++)
		for (int j=0; j<Ny; j++)
			for (int k=0; k<Nz/2+1; k++) {
				if (file_in >> temp_x){
					A(i,j,k) = (Complex)temp_x;
				}
				else {
					cout << endl << "ERROR: INPUT DATA < Required size. (" << i << "x" << j << "x" << k << ") points read." << endl;
					exit(1);
				}

			}
	cout << "Done. (" << Nx << "x" << Ny << "x" << Nz/2+1 << ") points read."<<endl;
}


int Convert_PRINFTY()
{
	cout << "Allocating 1 Array ..." << flush;
	Array<Complex, 3> A(Nx, Ny, Nz/2+1);
	Array<Complex, 3> B(Ny, Nz/2+1, Nx);
	cout << "Done" << endl;

	
	Read_complex_array(A, Nz, "T");

	Complex* A_data = A.data();
	Complex* B_data = B.data();

	for (int x=0; x<Nx; x++){
		for (int y=0; y<local_Ny; y++){
			for (int z=0; z<Nz/2+1; z++){
				*(B_data+y*Nx*(Nz/2+1) + z*Nx + x) = *(A_data + x*(Nz/2+1)*local_Ny + y*(Nz/2+1) + z);
			}
		}
	}

	BasicIO::Write(B_data, tarang_2_3.H5_full, BasicIO::data_out_folder, "T.F");

	return 0;
}

/*int Convert_FLUID()
{
	cout << "Allocating 3 Arrays ..." << flush;
	Array<Complex, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV3(local_Nx, Ny, Nz/2+1);
	cout << "Done" << endl;
	
	Read_complex_array(CV1, Nz, "CV1");
	Read_complex_array(CV2, Nz, "CV2");
	
	if (is_plane)
		Read_complex_array(CV3, 2, "CV3");
	else
		Read_complex_array(CV3, Nz/2+1, "CV3");
		
	if (is_plane)
		BasicIO::CV_Output_plane_HDF5(3, &CV1, &CV2, &CV3);
	else
		BasicIO::CV_Output_HDF5(3, &CV1, &CV2, &CV3);
	
	return 0;
}

int Convert_RBC()
{
	cout << "Allocating 4 Arrays ..." << flush;
	Array<Complex, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<Complex, 3> T(local_Nx, Ny, Nz/2+1);
	cout << "Done" << endl;
	
	
	
	Read_complex_array(CV1, Nz, "CV1");
	Read_complex_array(CV2, Nz, "CV2");
	
	if (is_plane)
		Read_complex_array(CV3, 2, "CV3");
	else
		Read_complex_array(CV3, Nz/2+1, "CV3");
	
	Read_complex_array(T, Nz, "T");



	if (is_plane)
		BasicIO::CV_Output_plane_HDF5(4, &CV1, &CV2, &CV3, &T);
	else
		BasicIO::CV_Output_HDF5(4, &CV1, &CV2, &CV3, &T);

	return 0;
}

int Convert_MHD()
{
	cout << "Allocating 6 Arrays ..." << flush;
	Array<Complex, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<Complex, 3> CW1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CW2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CW3(local_Nx, Ny, Nz/2+1);
	cout << "Done" << endl;
	
	
	Read_complex_array(CV1, Nz, "CV1");
	Read_complex_array(CV2, Nz, "CV2");
	
	if (is_plane)
		Read_complex_array(CV3, 2, "CV3");
	else
		Read_complex_array(CV3, Nz/2+1, "CV3");
		
	Read_complex_array(CW1, Nz, "CW1");
	Read_complex_array(CW2, Nz, "CW2");
	
	if (is_plane)
		Read_complex_array(CW3, 2, "CW3");
	else
		Read_complex_array(CW3, Nz/2+1, "CW3");
		
	
	if (is_plane)
		BasicIO::CV_Output_plane_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
	else
		BasicIO::CV_Output_HDF5(6, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3);
	
	return 0;
}

int Convert_MRBC()
{
	cout << "Allocating 7 Arrays ..." << flush;
	Array<Complex, 3> CV1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CV3(local_Nx, Ny, Nz/2+1);

	Array<Complex, 3> CW1(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CW2(local_Nx, Ny, Nz/2+1);
	Array<Complex, 3> CW3(local_Nx, Ny, Nz/2+1);

	Array<Complex, 3> T(local_Nx, Ny, Nz/2+1);
	cout << "Done" << endl;




	Read_complex_array(CV1, Nz, "CV1");
	Read_complex_array(CV2, Nz, "CV2");
	
	if (is_plane)
		Read_complex_array(CV3, 2, "CV3");
	else
		Read_complex_array(CV3, Nz/2+1, "CV3");
		
	Read_complex_array(CW1, Nz, "CW1");
	Read_complex_array(CW2, Nz, "CW2");
	
	if (is_plane)
		Read_complex_array(CW3, 2, "CW3");
	else
		Read_complex_array(CW3, Nz/2+1, "CW3");
		
	Read_complex_array(T, Nz, "T");

	


	if (is_plane)
		BasicIO::CV_Output_plane_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
	else
		BasicIO::CV_Output_HDF5(7, &CV1, &CV2, &CV3, &CW1, &CW2, &CW3, &T);
	return 0;
}
*/
