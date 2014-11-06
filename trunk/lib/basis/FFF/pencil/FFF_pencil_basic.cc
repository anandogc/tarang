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

/*! \file four_basic.cc
 * 
 * @sa four_basic.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	22/08/2008
 * @bug		No known bugs
 */ 

#include "FFF_pencil.h"
#include "FFF_pencil_inline.h"




/**********************************************************************************************
 
 COmment
 
 ***********************************************************************************************/
void FFF_PENCIL::Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz)
{}

void FFF_PENCIL::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
	DP Kx = kx*kfactor[1];
	DP Ky = ky*kfactor[2];
	DP Kz = kz*kfactor[3];
	
	if (kz != 0) 
		Vz = (complx(0,Kx)*Vx + complx(0,Ky)*Vy)/complx(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = (complx(0,Kx)*Vx)/complx(0,-Ky); 
		}
		else {	// k = (kx,0,0); input fields are (Vy, Vz); Vx=0
			Vz = Vy;
			Vy = Vx;
			Vx = complx(0,0);
		}
	}
}

/**********************************************************************************************
 
Dealias
 
 step 1: A(Range(Ny/3+1,2*Ny/3-1),:,:) = 0;
 
 step 2: A(:,Range(Nz/3+1,Nz/2),:) = 0;
 
 step 3: A(:,Range(0,Nz/3),Range(Nx/3+1,2*Nx/3-1)) = 0;
 
 ***********************************************************************************************/

void FFF_PENCIL::Dealias(Array<complx,3> A)
{
	Assign_sub_array(Range(Nx/3+1,2*Nx/3-1), Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,Nz/2), A, complx(0,0));
}


// Data resides till outer_radius in k-space
bool FFF_PENCIL::Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
		return true;
	
	return false;
}

/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 
 // For kz =0 and kz = N[3]/2 planes
 // f(-kx, -ky, 0) = conj(f(kx, ky, 0))- do for kz=N[3]/2
 // Implement:  f(kx, ky, 0 or Nz/2) = conj(f(-kx, -ky, 0 or Nz/2))
 // ky=Ny/2: set A(kx,ky,Nz/2)=0; quite efficient without losing many data points
 
 ***********************************************************************************************/

void FFF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<complx,3> A)
{
	int array_index_minus_kx, array_index_minus_ky;
	
	// For a given (minuskx, minusky), locate (kx,ky) and then subst.
	// A(minuskx, minusky, 0) = conj(A(kx,ky,0))
	if (my_z_pcoord == 0) {
		for (int lx=Nx/2+1; lx<Nx; lx++)  {
			array_index_minus_kx = -Get_kx(lx);         // kx<0; minuskx = -Get_kx(lx) > 0

			for (int ly=0; ly<maxly; ly++) {
				if (Get_ky(ly) != Ny/2) {  // Do not apply for ky=Ny/2
			
				array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly);
				

				A(lx,ly,0) = conj(global.temp_array.plane_xy(array_index_minus_ky,array_index_minus_kx));
				}
				// for (ky=0,kz=0) line
				if (Get_ky(ly) < 0)
					A(0,ly,0) = conj(global.temp_array.plane_xy(array_index_minus_ky,0));
			}
		}
	}

	// for kz=Nz/2
	if (my_z_pcoord == num_z_procs-1)
			A(Range::all(),Range::all(),maxlz-1) = 0.0;
	
}

void FFF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<complx,3> A)
{
	// for kz=Nz/2
	if (my_z_pcoord == num_z_procs-1)
		 A(Range::all(),Range::all(),maxlz-1) = 0.0;
}

void FFF_PENCIL::Test_reality_condition_in_Array(Array<complx,3> A)
{
	int array_index_minus_kx, array_index_minus_ky;
	
	// For a given (minuskx, minusky), locate (kx,ky) and then subst.
	// A(minuskx, minusky, 0) = conj(A(kx,ky,0))
	if (my_z_pcoord == 0) {
		Get_XY_plane(A, global.temp_array.plane_xy, 0);
		
		for (int lx=Nx/2+1; lx<Nx; lx++)  {
			array_index_minus_kx = -Get_kx(lx);         // kx<0; minuskx = -Get_kx(lx) > 0

			for (int ly=0; ly<maxly; ly++) {
				if (Get_ky(ly) != Ny/2) {  // Do not apply for ky=Ny/2
				
				array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly);
				
					
					if (abs(A(lx,ly,0)-conj(global.temp_array.plane_xy(array_index_minus_ky,array_index_minus_kx))) > MYEPS2)
						cout << "Reality condition voilated for (kx,ky,kz)=(" <<array_index_minus_kx <<  "," << Get_ky(ly) << "," << 0 << ")" << endl;
				}
				// for (ky=0,kz=0) line
				if (Get_ky(ly) < 0)
					if (abs(A(0,ly,0)-conj(global.temp_array.plane_xy(array_index_minus_ky,array_index_minus_kx))) > MYEPS2)
						cout << "Reality condition voilated for (kx,ky,kz)=(" << 0 <<  "," << Get_ky(ly) << "," << 0 << ")" << endl;
			}
		}
	}
	
	// for kz=Nz/2
	//int last_index=maxlz-1;
	if (my_z_pcoord == num_z_procs-1) {
		for (int lx=0; lx<Nx; lx++)
			for (int ly=0; ly<maxly; ly++)
				if (abs(A(ly,Nz/2,lx)) > MYEPS)
					cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << Get_ky(ly) << "," << Nz/2 << ")" << endl;
	}
}


void FFF_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<complx,3> A, complx value)
{
	static Array<int,1> y_filter(Ny);
	static Array<int,1> z_filter(Nz/2+1);
	
	y_filter = 0;
	z_filter = 0;
	
	y_filter(y_range)=1;
	z_filter(z_range)=1;

	
	static Range y_apply, z_apply;
	
	
	y_apply = Range(first(y_filter(Range(my_y_pcoord*maxly,(my_y_pcoord+1)*maxly-1)) == 1),
					 last(y_filter(Range(my_y_pcoord*maxly,(my_y_pcoord+1)*maxly-1)) == 1));
	
	z_apply = Range(first(z_filter(Range(my_z_pcoord*maxlz,(my_z_pcoord+1)*maxlz-1)) == 1),
					 last(z_filter(Range(my_z_pcoord*maxlz,(my_z_pcoord+1)*maxlz-1)) == 1));
	
	
	if ( (y_apply(0)>=0) && (z_apply(0)>=0))
		A(x_range, y_apply, z_apply) = value;
}

void FFF_PENCIL::Get_XY_plane(Array<complx,3> A, Array<complx,2> plane_xy, int kz)
{
	if (num_y_procs == 1) {
		plane_xy= A(Range::all(),Range::all(),kz);
		cout << "one proc" << endl;
	}
	
	else if (num_y_procs > 1) {
		int lz = universal->Get_lz(kz);
		
		global.temp_array.plane_xy_inproc = A(Range::all(), Range::all(), lz);

		int data_size = 2*maxly*Nx;
		int full_data_size = 2*Ny*Nx;
	
		MPI_Gather(reinterpret_cast<DP*>(global.temp_array.plane_xy_inproc.data()), data_size, MPI_DP, reinterpret_cast<DP*>(plane_xy.data()), 1, MPI_Vector_resized_y_plane_block, 0, global.mpi.MPI_COMM_ROW);
		
		MPI_Bcast(reinterpret_cast<DP*>(global.temp_array.plane_xy.data()), full_data_size, MPI_DP, 0, global.mpi.MPI_COMM_ROW);

/*		if (c>77) {
			for (int i=0; i<num_y_procs; i++) {
				cout << "my_y_pcoord, i  = " << my_y_pcoord << " " << i << endl;
				if (my_y_pcoord==i) {
					cout << "my_y_pcoord = " << my_y_pcoord << endl;
					cout << global.temp_array.plane_xy_inproc << endl;
					cout << global.temp_array.plane_xy << endl;
				}
				MPI_Barrier(global.mpi.MPI_COMM_ROW);
			}
		}*/
	}
}


int FFF_PENCIL::Read(Array<complx,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	int err = BasicIO::Read(global.temp_array.Xr.data(), plan, file_name, dataset_name);
	spectralTransform.Transpose(global.temp_array.Xr, A);
	return err;
}

int FFF_PENCIL::Read(Array<DP,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name)
{
	return BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
}


int FFF_PENCIL::Write(Array<complx,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	spectralTransform.Transpose(A, global.temp_array.Xr);
	return BasicIO::Write(global.temp_array.Xr.data(), plan, folder_name, file_name, dataset_name);
}

int FFF_PENCIL::Write(Array<DP,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name)
{
	return BasicIO::Write(Ar.data(), plan, folder_name, file_name, dataset_name);  
}


//*****************************  End of four_basic.cc *****************************************	








