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
void FFF_PENCIL::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
{}

void FFF_PENCIL::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{
	Real Kx = kx*kfactor[1];
	Real Ky = ky*kfactor[2];
	Real Kz = kz*kfactor[3];
	
	if (kz != 0) {
		Vz = (Complex(0,Kx)*Vx + Complex(0,Ky)*Vy)/Complex(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	}
	
	else {
		if (ky != 0) { // 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = (Complex(0,Kx)*Vx)/Complex(0,-Ky); 
		}
		else {	// k = (kx,0,0); input fields are (Vy, Vz); Vx=0
			Vz = Vy;
			Vy = Vx;
			Vx = Complex(0,0);
		}
	}
}

/**********************************************************************************************
 
Dealias
 
 step 1: A(Range(Ny/3+1,2*Ny/3-1),:,:) = 0;
 
 step 2: A(:,Range(Nz/3+1,Nz/2),:) = 0;
 
 step 3: A(:,Range(0,Nz/3),Range(Nx/3+1,2*Nx/3-1)) = 0;
 
 ***********************************************************************************************/

void FFF_PENCIL::Dealias(Array<Complex,3> A)
{
	Assign_sub_array(Range(Nx/3+1,2*Nx/3-1), Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,Nz/2), A, Complex(0,0));
}


// Data resides till outer_radius in k-space
bool FFF_PENCIL::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
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


void FFF_PENCIL::Get_XY_plane(Array<Complex,3> A, Array<Complex,2> plane_xy, int kz)
{
	if (num_y_procs == 1) {
		plane_xy= A(Range::all(),Range::all(),kz);
	}
	
	else if (num_y_procs > 1) {
		int lz = Get_lz(kz);
		
		global.temp_array.plane_xy_inproc = A(Range::all(), Range::all(), lz);

		int data_size = 2*maxly*Nx;
		int full_data_size = 2*Ny*Nx;
	
		MPI_Gather(reinterpret_cast<Real*>(global.temp_array.plane_xy_inproc.data()), data_size, MPI_Real, reinterpret_cast<Real*>(plane_xy.data()), 1, MPI_Vector_resized_y_plane_block, 0, global.mpi.MPI_COMM_ROW);
		
		MPI_Bcast(reinterpret_cast<Real*>(global.temp_array.plane_xy.data()), full_data_size, MPI_Real, 0, global.mpi.MPI_COMM_ROW);
	}
}

void FFF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	int array_index_minus_kx, array_index_minus_ky;
	
	if (my_z_pcoord == 0) {
		Get_XY_plane(A, global.temp_array.plane_xy, 0);

		#pragma ivdep
		// For a given (minuskx, minusky), locate (kx,ky) and then substitute.
		// A(minuskx, minusky, 0) = conj(A(kx,ky,0))
		for (int lx=Nx/2+1; lx<Nx; lx++)  {
			array_index_minus_kx = -Get_kx(lx);         // kx<0; minuskx = -Get_kx(lx) > 0

			for (int ly=0; ly<maxly; ly++) {
				if (Get_ky(ly) != Ny/2) {  // Do not apply for ky=Ny/2
			
					array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly);
					

					A(lx,ly,0) = conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky));
				}

			}
		}
	
		// for (kx=0,kz=0) line
		// For a given (0, minusky), locate (0,ky) and then substitute.
		for (int ly=0; ly<maxly; ly++) {
			if (Get_ky(ly) < 0) {
				array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly);
				A(0,ly,0) = conj(global.temp_array.plane_xy(0,array_index_minus_ky));
			}
		}
		
	}

	// for kz=Nz/2
	if (my_z_pcoord == num_z_procs-1)
			A(Range::all(),Range::all(),maxlz-1) = 0.0;
	
}

void FFF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
	// for kz=Nz/2
	if (my_z_pcoord == num_z_procs-1)
		 A(Range::all(),Range::all(),maxlz-1) = 0.0;
}

void FFF_PENCIL::Test_reality_condition_in_Array(Array<Complex,3> A)
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
					
						
					if (abs(A(lx,ly,0)-conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky))) > MYEPS2)
							cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << Get_ky(ly) << "," << 0 << ")" << endl;
				}
				// for (ky=0,kz=0) line
				if (Get_ky(ly) < 0)
					if (abs(A(0,ly,0)-conj(global.temp_array.plane_xy(array_index_minus_kx,array_index_minus_ky))) > MYEPS2)
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


void FFF_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<Complex,3> A, Complex value)
{
	static Array<int,1> y_filter(Ny);
	static Array<int,1> z_filter(Nz/2+1);
	
	//Sanitize ranges. if last index is less than first index, modify the range (this happens for 2D)
	if (y_range.last() < y_range.first())
		y_range = Range(y_range.first(), y_range.first());

	if (z_range.last() < z_range.first())
		z_range = Range(z_range.first(), z_range.first());

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

int FFF_PENCIL::Read(Array<Complex,3> A, h5::Plan plan, string file_name, string dataset_name)
{
	if (Ny > 1) {
		if (my_z_pcoord == 0)  
			BasicIO::Read(global.temp_array.Xr_slab.data(), plan, file_name, dataset_name);

		fftk.To_pencil(global.temp_array.Xr_slab, global.temp_array.Xr);
		fftk.Transpose(global.temp_array.Xr, A);

	}
	if (Ny == 1) {
		BasicIO::Read(global.temp_array.Xr.data(), plan, file_name, dataset_name);
		fftk.Transpose(global.temp_array.Xr(Range::all(),0,Range::all()), A(Range::all(),0,Range::all()));
	}
	return 0;
}

int FFF_PENCIL::Read(Array<Real,3> Ar, h5::Plan plan, string file_name, string dataset_name)
{
	if (Ny > 1) {
		if (my_z_pcoord == 0)  
			BasicIO::Read(global.temp_array.Xr_slab.data(), plan, file_name, dataset_name);

		fftk.To_pencil(global.temp_array.Xr_slab, Ar);
	}
	if (Ny == 1) {
		BasicIO::Read(Ar.data(), plan, file_name, dataset_name);
	}
	return 0;
}


int FFF_PENCIL::Write(Array<Complex,3> A, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
{
	if (Ny > 1) {
		fftk.Transpose(A, global.temp_array.Xr);
		fftk.To_slab(global.temp_array.Xr, global.temp_array.Xr_slab);
		if (my_z_pcoord == 0)  
			BasicIO::Write(global.temp_array.Xr_slab.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}
	if (Ny == 1) {
		fftk.Transpose(A(Range::all(),0,Range::all()), global.temp_array.Xr(Range::all(),0,Range::all()));
		BasicIO::Write(global.temp_array.Xr.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}

	return 0;
}

int FFF_PENCIL::Write(Array<Real,3> Ar, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
{
	if (Ny > 1) {
		fftk.To_slab(Ar, global.temp_array.Xr_slab);
		if (my_z_pcoord == 0)  
			BasicIO::Write(global.temp_array.Xr_slab.data(), plan, access_mode, folder_name, file_name, dataset_name);

	}
	if (Ny == 1) {
		BasicIO::Write(Ar.data(), plan, access_mode, folder_name, file_name, dataset_name);
	}
	return 0;
}


//*****************************  End of four_basic.cc *****************************************	








