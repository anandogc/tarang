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
#include "TransposeConfig.h"


//*********************************************************************************************	

class SpectralTransform
{
public:

	int my_id;
	int numprocs;
	
	//For Pencil
	int num_p_vert;
	int num_p_hor;
	int my_vert_pcoord;
	int my_hor_pcoord;
	
	
	int Nx, Ny, Nz;

	
	Array<complx,3> X_3d;
	Array<DP,3> Xr_3d;
	Array<complx,3> A3d_interm;
	
	Array<complx,2> X_2d;
	Array<DP,2> Xr_2d;
	Array<complx,2> A2d_interm;

	//For Transpose

	string decomposition;

	enum {
		YZ_PLANE,
		ZX_PLANE,
		X,
		Y,
		Z,
	};

	//Slab
	TinyVector<int, 3> shape_vertical_array_3d;
	TinyVector<int, 3> shape_horizontal_array_3d;
	
	TinyVector<int, 2> shape_vertical_array_2d;
	TinyVector<int, 2> shape_horizontal_array_2d;
	
	//Pencil
	TinyVector<int, 3> shape_x_array_3d;
	TinyVector<int, 3> shape_y_array_3d;
	TinyVector<int, 3> shape_z_array_3d;

	TransposeConfig XY, YZ, ZX;
	TransposeConfig XZ, ZY, YX;
	TransposeConfig XY_real, YX_real; //used in 3D ChFF
	TransposeConfig XY_complex, YX_complex; //used 3D in ChFF
	TransposeConfig XZ_real, ZX_real; //used in 2D ChFF
	TransposeConfig XZ_complex, ZX_complex; //used in 2D ChFF
	// **
	
	ptrdiff_t local_Nx_start;
	ptrdiff_t local_Ny_start;
	ptrdiff_t local_Nz_start;
	
	//For Slab
	ptrdiff_t local_Nx;
	ptrdiff_t local_Ny;
	ptrdiff_t local_Nz;
    
	//For Pencil
	ptrdiff_t local_Nx_vert;
	ptrdiff_t local_Ny_vert;
	ptrdiff_t local_Nz_vert;
	
	ptrdiff_t local_Nx_hor;
	ptrdiff_t local_Ny_hor;
	ptrdiff_t local_Nz_hor;
	
	MPI_Comm MPI_COMM_CART;
	MPI_Comm MPI_COMM_HOR_SEGMENT;
	MPI_Comm MPI_COMM_VERT_SEGMENT;
	
	
	
	// for slab fftw and GP
	FFTW_PLAN c2c_xyz_forward_plan;
	FFTW_PLAN c2c_xyz_inverse_plan;
	
	// for FT
	FFTW_PLAN c2c_x_forward_plan, c2c_x_inverse_plan;
	FFTW_PLAN c2c_y_forward_plan, c2c_y_inverse_plan;	

	FFTW_PLAN r2c_xyz_plan;
	FFTW_PLAN c2r_xyz_plan;

	FFTW_PLAN r2c_xz_plan;
	FFTW_PLAN c2r_xz_plan;

	FFTW_PLAN r2c_yz_plan;
	FFTW_PLAN c2r_yz_plan;

	FFTW_PLAN r2c_z_plan;
	FFTW_PLAN c2r_z_plan;
	
	
	// for ST/CT
	FFTW_PLAN sintr_x_plan, sintr_y_plan, sintr_z_plan;
	FFTW_PLAN isintr_x_plan, isintr_y_plan, isintr_z_plan;
	FFTW_PLAN costr_x_plan, costr_y_plan, costr_z_plan;
	FFTW_PLAN icostr_x_plan, icostr_y_plan, icostr_z_plan;
	FFTW_PLAN Chebyshevtr_x_plan;

	
	template<class T1, class T2, int N_rank>
	void Transpose_array(Array<T1,N_rank> A, Array<T2,N_rank> B, TransposeConfig config);

	template <class T>
	void Set_transpose_config(int send_config_id, int recv_config_id, MPI_Comm communicator, TransposeConfig& config);

	void FTc2c_xyz(Array<complx,3> A1, Array<complx,3> A2);
	void IFTc2c_xyz(Array<complx,3> A2, Array<complx,3> A1);

	// 3D: FFT original
	void FTr2c_xyz(Array<DP,3> Ar, Array<complx,3> A);
	void FTc2r_xyz(Array<complx,3> A, Array<DP,3> Ar);
	
	// 2D
	void FTr2c_xz(Array<DP,2> Ar, Array<complx,2> A);
	void FTc2r_xz(Array<complx,2> A, Array<DP,2> Ar);
	void FTr2c_yz(Array<complx,2> Plane);
	// for GP
	void FTc2c_xyz(Array<complx,3> A);
	void IFTc2c_xyz(Array<complx,3> A);

	void FTr2c_yz(Array<DP,3> Xr);
	void FTc2r_yz(Array<DP,3> Xr);
	
	void FTr2c_yz(Array<complx,3> X);
	void FTc2r_yz(Array<complx,3> X);
	
	// 1D
	// c2c
	void FTc2c_x(Array<complx,3> A);
	void IFTc2c_x(Array<complx,3> A);
	
	void FTc2c_y(Array<complx,3> A);
	void IFTc2c_y(Array<complx,3> A);
	
	// r2c, r2c
	void FTr2c_z(Array<complx,3> A);
	void FTc2r_z(Array<complx,3> A);
	
	void FTr2c_z(Array<DP,3> A);
	void FTc2r_z(Array<DP,3> A);
	
	void FTr2c_z(Array<complx,2> A);
	void FTc2r_z(Array<complx,2> A);
	
	void FTr2c_z(Array<DP,2> A);
	void FTc2r_z(Array<DP,2> A);

	// 3D only
	void SinCostr_x(char sincostr_switch, Array<DP,3> Ar);
	void ISinCostr_x(char sincostr_switch, Array<DP,3> Ar);

	void SinCostr_y(char sincostr_switch, Array<complx,3> A);
	void ISinCostr_y(char sincostr_switch, Array<complx,3> A);
	
	void Chebyshevtr_x(Array<DP,3> A);
	void IChebyshevtr_x(Array<DP,3> A);
	
	void SinCostr_z(char sincostr_switch, Array<complx,3> A);
	void ISinCostr_z(char sincostr_switch, Array<complx,3> A);
	
	// 2D
	void SinCostr_x(char sincostr_switch, Array<DP,2> Ar);
	void ISinCostr_x(char sincostr_switch, Array<DP,2> Ar);
	void Chebyshevtr_x(Array<DP,2> A);
	void IChebyshevtr_x(Array<DP,2> A);
	void SinCostr_z(char sincostr_switch, Array<complx,2> A);
	void ISinCostr_z(char sincostr_switch, Array<complx,2> A);
	
	
	
	void ArrayShiftLeft(Array<DP,3> Ar, char shift_dirn);
	void ArrayShiftRight(Array<DP,3> A, char shift_dirn);
	void ArrayShiftLeft(Array<complx,3> A, char shift_dirn);
	void ArrayShiftRight(Array<complx,3> A, char shift_dirn);
	

	void ArrayShiftLeft(Array<DP,2> Ar, char shift_dirn);
	void ArrayShiftRight(Array<DP,2> A, char shift_dirn);
	void ArrayShiftRight(Array<complx,2> A, char shift_dirn);
	void ArrayShiftLeft(Array<complx,2> A, char shift_dirn);
	
	
	
	void Init(string basis, string docomposition, int Nx, int Ny, int Nz, int num_p_hor=0);

	//CFFF_SLAB
	void Init_CFFF_SLAB();
	//3D
	void Zero_pad_last_plane_CFFF_SLAB(Array<complx,3> Ar);
	void Norm_CFFF_SLAB(Array<complx,3> A);
	void Forward_transform_CFFF_SLAB(Array<complx,3> Ar, Array<complx,3> A);
	void Inverse_transform_CFFF_SLAB(Array<complx,3> A, Array<complx,3> Ar);

	//2D
	void Zero_pad_last_col_CFFF_SLAB(Array<complx,2> Ar);
	void Norm_CFFF_SLAB(Array<complx,2> A);
	void Forward_transform_CFFF_SLAB(Array<complx,2> Ar, Array<complx,2> A);
	void Inverse_transform_CFFF_SLAB(Array<complx,2> A, Array<complx,2> Ar);

	//FFFW_SLAB
	void Init_FFFW_SLAB();
	//3D
	void Zero_pad_last_plane_FFFW_SLAB(Array<DP,3> Ar);
	void Norm_FFFW_SLAB(Array<complx,3> A);
	void Forward_transform_FFFW_SLAB(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_FFFW_SLAB(Array<complx,3> A, Array<DP,3> Ar);

	//2D
	void Zero_pad_last_col_FFFW_SLAB(Array<DP,2> Ar);
	void Norm_FFFW_SLAB(Array<complx,2> A);
	void Forward_transform_FFFW_SLAB(Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_FFFW_SLAB(Array<complx,2> A, Array<DP,2> Ar);


	//FFF_SLAB
	void Init_FFF_SLAB();
	//3D
	void Zero_pad_last_plane_FFF_SLAB(Array<DP,3> Ar);
	void Norm_FFF_SLAB(Array<complx,3> A);
	void Forward_transform_FFF_SLAB(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_FFF_SLAB(Array<complx,3> A, Array<DP,3> Ar);

	//2D
	void Zero_pad_last_col_FFF_SLAB(Array<DP,2> Ar);
	void Norm_FFF_SLAB(Array<complx,2> A);
	void Forward_transform_FFF_SLAB(Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_FFF_SLAB(Array<complx,2> A, Array<DP,2> Ar);
	
	
	//SFF_SLAB
	void Init_SFF_SLAB();
	//3D
	void Zero_pad_last_plane_SFF_SLAB(Array<DP,3> Ar);
	void Norm_SFF_SLAB(Array<complx,3> A);
	void Forward_transform_SFF_SLAB(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SFF_SLAB(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);

	//2D
	void Zero_pad_last_col_SFF_SLAB(Array<DP,2> Ar);
	void Norm_SFF_SLAB(Array<complx,2> A);
	void Forward_transform_SFF_SLAB(string sincostr_switch, Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_SFF_SLAB(string sincostr_switch, Array<complx,2> A, Array<DP,2> Ar);
	
	//SSF_SLAB
	void Init_SSF_SLAB();
	//3D
	void Zero_pad_last_plane_SSF_SLAB(Array<DP,3> Ar);
	void Norm_SSF_SLAB(Array<complx,3> A);
	void Forward_transform_SSF_SLAB(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SSF_SLAB(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);
	
	//2D
	void Zero_pad_last_col_SSF_SLAB(Array<DP,2> Ar);
	void Norm_SSF_SLAB(Array<complx,2> A);
	void Forward_transform_SSF_SLAB(string sincostr_switch, Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_SSF_SLAB(string sincostr_switch, Array<complx,2> A, Array<DP,2> Ar);
	
	//SSS_SLAB
	void Init_SSS_SLAB();
	//3D
	void Zero_pad_last_plane_SSS_SLAB(Array<DP,3> Ar);
	void Norm_SSS_SLAB(Array<complx,3> A);
	void Forward_transform_SSS_SLAB(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SSS_SLAB(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);
	
	//2D
	void Zero_pad_last_col_SSS_SLAB(Array<DP,2> Ar);
	void Norm_SSS_SLAB(Array<complx,2> A);
	void Forward_transform_SSS_SLAB(string sincostr_switch, Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_SSS_SLAB(string sincostr_switch, Array<complx,2> A, Array<DP,2> Ar);
	
	//ChFF_SLAB
	void Init_ChFF_SLAB();
	//3D
	void Zero_pad_last_plane_ChFF_SLAB(Array<DP,3> Ar);
	void Norm_ChFF_SLAB(Array<complx,3> A);
	void Forward_transform_ChFF_SLAB(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_ChFF_SLAB(Array<complx,3> A, Array<DP,3> Ar);
	
	//2D
	void Zero_pad_last_col_ChFF_SLAB(Array<DP,2> Ar);
	void Norm_ChFF_SLAB(Array<complx,2> A);
	void Forward_transform_ChFF_SLAB(Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform_ChFF_SLAB(Array<complx,2> A, Array<DP,2> Ar);
	
	
	
	//Pencil
	
	//FFF_PENCIL
	void Init_FFF_PENCIL();
	//3D
	void Zero_pad_last_plane_FFF_PENCIL(Array<DP,3> Ar);
	void Norm_FFF_PENCIL(Array<complx,3> A);
	void Forward_transform_FFF_PENCIL(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_FFF_PENCIL(Array<complx,3> A, Array<DP,3> Ar);
	
	//SFF_PENCIL
	void Init_SFF_PENCIL();
	//3D
	void Zero_pad_last_plane_SFF_PENCIL(Array<DP,3> Ar);
	void Norm_SFF_PENCIL(Array<complx,3> A);
	void Forward_transform_SFF_PENCIL(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SFF_PENCIL(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);

	//SSF_PENCIL
	void Init_SSF_PENCIL();
	//3D
	void Zero_pad_last_plane_SSF_PENCIL(Array<DP,3> Ar);
	void Norm_SSF_PENCIL(Array<complx,3> A);
	void Forward_transform_SSF_PENCIL(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SSF_PENCIL(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);

	
	//SSS_PENCIL
	void Init_SSS_PENCIL();
	//3D
	void Zero_pad_last_plane_SSS_PENCIL(Array<DP,3> Ar);
	void Norm_SSS_PENCIL(Array<complx,3> A);
	void Forward_transform_SSS_PENCIL(string sincostr_switch, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform_SSS_PENCIL(string sincostr_switch, Array<complx,3> A, Array<DP,3> Ar);

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


