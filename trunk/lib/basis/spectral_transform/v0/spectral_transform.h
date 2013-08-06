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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _SPECTRAL_TRANSFORM_H
#define _SPECTRAL_TRANSFORM_H

#include "spectral_def_vars.h"


//*********************************************************************************************	

class SpectralTransform
{
public:
    int my_id;
    int numprocs;
    
    int Nx, Ny, Nz;

    
    Array<complx,3> X;
    Array<complx,3> Xr;
    

    bool fftw3D_switch;					// for fftw_original
    
    
    ptrdiff_t local_Nx_start;
    ptrdiff_t local_Ny_start;
    ptrdiff_t local_Nz_start;
    
    //For Slab
    ptrdiff_t local_Nx;
    ptrdiff_t local_Ny;
    ptrdiff_t local_Nz;

    
    // original
    FFTW_PLAN_DP r2c_xyz, c2r_xyz;
    
    
    // for slab fftw
    
    // for FT
    FFTW_PLAN_DP c2r_z_plan;
    FFTW_PLAN_DP r2c_z_plan;
    FFTW_PLAN_DP c2c_x_forward_plan, c2c_y_forward_plan;
    FFTW_PLAN_DP c2c_x_inverse_plan, c2c_y_inverse_plan;
    
    FFTW_PLAN_DP r2c_xyz_plan;
    FFTW_PLAN_DP c2r_xyz_plan;
    FFTW_PLAN_DP r2c_xz_plan;
    FFTW_PLAN_DP r2c_yz_plan;
    
    
    
    FFTW_PLAN_DP c2r_xz_plan, c2r_yz_plan;
    
    // for ST/CT
    FFTW_PLAN_DP sintr_x_plan, sintr_y_plan, sintr_z_plan;
    FFTW_PLAN_DP isintr_x_plan, isintr_y_plan, isintr_z_plan;
    FFTW_PLAN_DP costr_x_plan, costr_y_plan, costr_z_plan;
    FFTW_PLAN_DP icostr_x_plan, icostr_y_plan, icostr_z_plan;
    FFTW_PLAN_DP Chebyshevtr_x_plan;
    
    // for GP
    FFTW_PLAN_DP c2c_forward_plan, c2c_inverse_plan;
    

    
    // fft arrays
    // if 3D
    Array<complx,1> col_x;      // Nx or N1
    Array<complx,1> col_y;		// Ny or N2
    Array<complx,1> col_z_input;		// Nz/2
    Array<complx,1> col_z;		// Nz/2+1
    
    
    Array<DP,1> col_x_real;          // Nx or N1
    Array<DP,1> col_y_real;		// Ny or N2
    Array<DP,1> col_z_real;     // Nz for free-slip
    // DOES IT WORK? should it be Array<complx,1> col_zfs size Nz/2
    
    Array<DP,1> col_x_Chebyshev;          // Nx or N1
    
    Array<complx,2> plane_xz;  // N1,N3/2+1
    Array<complx,2> plane_yz;  // N2,N3/2+1
    
    Array<complx,3> Ayzx;
    Array<complx,3> Azyx;			// SSS
    
    Array<complx,3> Ayxz;
    
    Array<complx,3> Axyz_Chff;
    
    
    // for transpose
    struct slab_transpose {
        // 3D SLAB
        Array<complx,2> AxyTR;			// (Ny, local_Nx);
        Array<complx,2> AxyTR_recv;		// (Ny, local_Nx);
        
        Array<complx,2> AyxTR;			// (Nx, local_Ny);
        Array<complx,2> AyxTR_recv;		// (Nx, local_Ny);
        
        // 2D SLAB
        Array<complx,2> AxzTR;			// (Nz/2, localNx)
        Array<complx,2> AxzTR_recv;		//	(Nz/2, localNx)
        
        Array<complx,2> AzxTR;			// (Nx, localNz);
        Array<complx,2> AzxTR_recv;		// (Nx, localNz);
    } slab_transpose;
    
    
	void Transpose_array_SLAB(Array<complx,3> A, Array<complx,3> B, char X1, char X2, char X3, char fixedaxis);
	void Transpose_array_SLAB(Array<complx,2> A, Array<complx,2> B, char X1, char X2);

	
	// 3D: FFT original
	void FTr2c_xyz(Array<complx,3> A);
	void FTc2r_xyz(Array<complx, 3>);
	
	// 2D
	void FTr2c_xz(Array<complx,2> Plane);
	void FTr2c_yz(Array<complx,2> Plane);
	void FTc2r_xz(Array<complx,2> Plane);
	void FTc2r_yz(Array<complx,2> Plane);
	
	// 1D
	// c2c
	void FTc2c_x(Array<complx,1> Col);
	void FTc2c_y(Array<complx,1> Col);
	void IFTc2c_x(Array<complx,1> Col);
	void IFTc2c_y(Array<complx,1> Col);
	
	// r2c, r2c
	void FTr2c_z(Array<complx,1> Col);
	void FTc2r_z(Array<complx,1> Col);
	void Sintr_x(Array<DP,1> Col);
	void Sintr_y(Array<DP,1> Col);
	void Sintr_z(Array<DP,1> Col);
	void Costr_x(Array<DP,1> Col);
	void Costr_y(Array<DP,1> Col);
	void Costr_z(Array<DP,1> Col);
	
	void ISintr_x(Array<DP,1> Col);
	void ISintr_y(Array<DP,1> Col);
	void ISintr_z(Array<DP,1> Col);
	void ICostr_x(Array<DP,1> Col);
	void ICostr_y(Array<DP,1> Col);
	void ICostr_z(Array<DP,1> Col);
	
	void Chebyshevtr_x(Array<DP,1> Col);
	void IChebyshevtr_x(Array<DP,1> Col);
    
    // for GP
    void FTc2c_xyz(Array<complx,3> A);
    void IFTc2c_xyz(Array<complx,3> A);
	
	// Presently along x.

	void ArrayShiftLeft_SLAB(Array<complx,2> A, char shift_dirn);
	void ArrayShiftRight_SLAB(Array<complx,2> A, char shift_dirn);
	
	void ArrayShiftLeft_SLAB(Array<complx,3> A, char Xa, char Xb, char Xc, char shift_dirn);
	void ArrayShiftRight_SLAB(Array<complx,3> A, char Xa, char Xb, char Xc, char shift_dirn);
    
    
    
    void Init(string basis, string docomposition, int Nx, int Ny, int Nz, bool fftw3D_switch=false);
    //FFF_SLAB
    void Init_FFF_SLAB();
    //3D
    void Zero_pad_last_plane_FFF_SLAB(Array<complx,3> Ar);
    void Norm_FFF_SLAB(Array<complx,3> A);
    void Forward_transform_array_transpose_order_FFF_SLAB(Array<complx,3> Ar, Array<complx,3> A);
    void Forward_transform_array_FFF_SLAB(Array<complx,3> A);
    void Inverse_transform_array_transpose_order_FFF_SLAB(Array<complx,3> A, Array<complx,3> Ar);
    void Inverse_transform_array_FFF_SLAB(Array<complx,3> A);
    //2D
    void Zero_pad_last_plane_FFF_SLAB(Array<complx,2> Ar);
    void Norm_FFF_SLAB(Array<complx,2> A);
    void Forward_transform_array_transpose_order_FFF_SLAB(Array<complx,2> Ar, Array<complx,2> A);
    void Forward_transform_array_FFF_SLAB(Array<complx,2> A);
    void Inverse_transform_array_transpose_order_FFF_SLAB(Array<complx,2> A, Array<complx,2> Ar);
    void Inverse_transform_array_FFF_SLAB(Array<complx,2> A);
    
    //SFF_SLAB
    void Init_SFF_SLAB();
    //3D
    void Zero_pad_last_plane_SFF_SLAB(Array<complx,3> Ar);
    void Norm_SFF_SLAB(Array<complx,3> A);
    void Forward_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,3> Ar, Array<complx,3> A);
    void Forward_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,3> A);
    void Inverse_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,3> A, Array<complx,3> Ar);
    void Inverse_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,3> A);
    
    //2D
    void Zero_pad_last_plane_SFF_SLAB(Array<complx,2> Ar);
    void Norm_SFF_SLAB(Array<complx,2> A);
    void Forward_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,2> Ar, Array<complx,2> A);
    void Forward_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,2> A);
    void Inverse_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,2> A, Array<complx,2> Ar);
    void Inverse_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,2> A);
    
    

	/** @brief Transpose(A) -> Atr.  Only x-y components.
	 *
	 * Steps in each processor
	 * (1) Take a xy plane Axy for a given z (in each processor).
	 * (2) Transpose Axy in each processor. Axy -> Axy_tr.
	 * (3) MPI_Alltoall(Axy_tr -> temp2)
	 * (4) Transpose each transmitted piece in each processor.
	 *
	 * @param  A() complex array
	 * @param  N[] size of A: (N1, N2, N3/2+1).
	 * 
	 * @return   Atr():  Transpose of A along xy directions.
	 *
	 * @sa figure ..
	 *
	 */
		
//*********************************************************************************************
	
	/** @brief Inverse_transpose(Atr) -> A.  Only x-y components.
	 *
	 * Steps in each processor
	 * (1) Take a xy plane Atr_xy for a given z (in each processor).
	 * (2) Transpose Atr_xy in each processor. Atr_xy -> Axy.
	 * (3) MPI_Alltoall(Axy -> temp1)
	 * (4) Transpose each transmitted piece in each processor.
	 *
	 * @param  Atr() complex array
	 * @param  N[] size of A: (N1, N2, N3/2+1).
	 * 
	 * @return   A():  Transpose of Atr along xy directions.
	 *
	 * @sa figure ..
	 *
	 */
	
		// pencil
	
};

#endif


//******************************** End of field_basic.h  **************************************


