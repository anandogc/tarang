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

/*! \file scft_basic.cc
 *
 * @sa scft_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0 
 * @date 30 August 2008
 * @bug  No known bug
 */
 

#include "SFF_pencil.h"
#include "SFF_pencil_inline.h"


 
/**********************************************************************************************
 
 COmment
 
***********************************************************************************************/
void SFF_PENCIL::Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz)
{}

void SFF_PENCIL::Last_component(int kx, int ky, int kz, Complex &Vx, Complex &Vy, Complex &Vz)
{
	Real Kx = kx*kfactor[1];
	Real Ky = ky*kfactor[2];
	Real Kz = kz*kfactor[3];
	
	Complex dvxdx;
	
	if (global.program.sincostr_switch[0] == 'S')
		dvxdx = Complex(Kx,0)*Vx;
	
	else if (global.program.sincostr_switch[0] == 'C') 
		dvxdx = Complex(-Kx,0)*Vx;
	
	if (kz != 0) 
		Vz = (dvxdx + Complex(0,Ky)*Vy)/Complex(0,-Kz);
		// works for both 2D (Ky=0) and 3D.
	
	else {
		if (ky != 0) {// 3D: input fields are (Vx, Vz); compute Vy
			Vz = Vy;
			Vy = dvxdx/Complex(0,-Ky); 
		}
		else { // k = (kx,0,0)
			if ( (abs(imag(Vx)) > MYEPS) || (abs(imag(Vx)) > MYEPS))
				cout << "MYERROR: SCFT_SLAB::Last_component(); For (kx,0,0), the modes are purely real; Setting imag part to zero. " << endl;
			Vz = Complex(real(Vy), 0);
			Vy = Complex(real(Vx), 0);
			Vx = Complex(0,0);
		}
	}
	
}


/*******************************************************************************
 
 Dealias

	A(Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd), Range(2*Nx/3+1,toEnd)) = 0;
 
 *******************************************************************************/

void SFF_PENCIL::Dealias(Array<Complex,3> A)
{
	Assign_sub_array(Range(2*Nx/3+1,toEnd), Range(Ny/3+1,2*Ny/3-1), Range(Nz/3+1,toEnd), A, Complex(0,0));
}


// Data resides till outer_radius in k-space
bool SFF_PENCIL::Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius)
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if ((kx_max > 2*Nx/3) || (ky_max > Ny/3) || (kz_max > Nz/3))
		return true;

	return false;
}

void SFF_PENCIL::Get_XY_plane(Array<Complex,3> A, Array<Complex,2> plane_xy, int kz)
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

/**********************************************************************************************
 
 Satisfy_reality_condition_in_Array
 Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
 No need for data transfer; local..
 
f(m, -ky, 0) = conj(f(m, ky, 0))- do for kz=N[3]/2
 
 ***********************************************************************************************/

void SFF_PENCIL::Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A)
{
	int array_index_minus_ky;
	
	if (my_z_pcoord == 0) {
		Get_XY_plane(A, global.temp_array.plane_xy, 0);

		// For a given (kx, minusky), locate (kx,ky) and then subst.
		// A(kx, minusky, 0) = conj(A(kx,ky,0))
		#pragma ivdep
		for (int lx=0; lx<Nx; lx++)  {
			for (int ly=0; ly<maxly; ly++) {
				if (Get_ky(ly) <= 0) {
			
					array_index_minus_ky =  Get_iy(-Get_ky(ly));  // minusky = -Get_ky(ly);
					
					A(lx,ly,0) = conj(global.temp_array.plane_xy(lx,array_index_minus_ky));
				}
			}
		}
	}

	if (my_z_pcoord == num_z_procs-1)
         A(Range::all(),Range::all(),maxlz-1) = 0.0;
}

void SFF_PENCIL::Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A)
{
    if (my_z_pcoord == num_z_procs-1)
        A(Range::all(),Range::all(),maxlz-1) = 0.0;
}



void SFF_PENCIL::Test_reality_condition_in_Array(Array<Complex,3> A)
{
    int array_index_minus_ky;
	
	if (my_z_pcoord == 0) {
        for (int lx=0; lx<maxlx; lx++) {
           for (int ly=Ny/2+1; ly<Ny; ly++) {
                array_index_minus_ky = -Get_ky(ly);
               
                if (abs(A(lx,ly,0)-conj(A(lx,array_index_minus_ky,0))) > MYEPS2)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << ly << "," << 0 << ")" << endl; 
            }
            // for (ky=0,kz=0) line
            if (abs(imag(A(lx,0,0))) > MYEPS2)
                cout << "Reality condition voilated for (kx,ky,kz)=(" << lx <<  "," << 0 << "," << 0 << ")" << endl;
        }
	}
	
    // for kz=Nz/2
    int last_index=maxlz-1;
	if (my_z_pcoord == num_z_procs-1) {
        for (int lx=0; lx<maxlx; lx++)
			for (int ly=0; ly<Ny; ly++)
                if (abs(A(lx,ly,last_index)) > MYEPS)
                    cout << "Reality condition voilated for (kx,ky,kz)=(" << Get_kx(lx) <<  "," << Get_ky(ly) << "," << Nz/2 << ")" << endl;
    }

}


//******************************************************************************
/** @brief Set the modes to zero for (0,ky,kz) in 3D.
 *
 * @return SFF: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return CFF: \f$ V_2(0,ky,kz) = V_3(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 */

void SFF_PENCIL::Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
    // lx = 0 reside in master node
    
    Range zero(0,0);
    global.program.sincostr_switch = sincostr_switch_Vx;

	if (my_x_pcoord == 0) {
        if (global.program.sincostr_switch == "SFF")
        	Assign_sub_array(zero,Range::all(),Range::all(),Ax,Complex(0,0));
        
        else if (global.program.sincostr_switch == "CFF") {
        	Assign_sub_array(zero,Range::all(),Range::all(),Ay,Complex(0,0));
        	Assign_sub_array(zero,Range::all(),Range::all(),Az,Complex(0,0));
        }
    }
}

/** @brief Set the modes to zero for T.F(0,ky,kz) in 3D.
 *
 * Temperature same as V_1; (see the above function).
 */    
void SFF_PENCIL::Zero_modes(Array<Complex,3> F)
{
    // lx = 0 reside in master node
    Range zero(0,0);
    global.program.sincostr_switch = sincostr_switch_F;
    if ((my_x_pcoord == 0) && (global.program.sincostr_switch == "SFF"))
    	Assign_sub_array(zero,Range::all(),Range::all(),F,Complex(0,0));
}

void SFF_PENCIL::Assign_sub_array(Range x_range, Range y_range, Range z_range, Array<Complex,3> A, Complex value)
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


int SFF_PENCIL::Read(Array<Complex,3> A, h5::Plan plan, string file_name, string dataset_name)
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

int SFF_PENCIL::Read(Array<Real,3> Ar, h5::Plan plan, string file_name, string dataset_name)
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


int SFF_PENCIL::Write(Array<Complex,3> A, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
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

int SFF_PENCIL::Write(Array<Real,3> Ar, h5::Plan plan, string access_mode, string folder_name, string file_name, string dataset_name)
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

//*********************************  End of scft_basic.cc *************************************


