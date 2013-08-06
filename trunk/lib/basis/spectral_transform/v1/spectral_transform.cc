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
/*! \file fourier.cc 
 * 
 * @sa fourier.hFTr2c_yz(plane_yz);
 * 
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "spectral_transform.h"


//*********************************************************************************************

void SpectralTransform::Init(string basis, string decomposition, int Nx, int Ny, int Nz, int num_p_hor) {
	
    this->Nx=Nx;
    this->Ny=Ny;
    this->Nz=Nz;
	
    this->decomposition = decomposition;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (decomposition == "PENCIL" && numprocs>1)
	{
		this->num_p_hor = num_p_hor;
		if (num_p_hor>0){
			if (numprocs%num_p_hor==0){
				num_p_vert = numprocs/num_p_hor;
				
			}
			else{
				if (my_id==0)
					cerr << "SpectralTransform::Init - numprocs(=" << numprocs << ") must be divisible by num_p_hor(=" << num_p_hor << ")." << endl;
				MPI_Finalize();
				exit(1);
			}
		}
		
		int dim_size[2];
		dim_size[0] = num_p_vert;
		dim_size[1] = num_p_hor;
		
		int periods[2];
		int reorder = 0;
		int ndims = 2;
		
		periods[0] = 0;
		periods[1] = 0;
        
		MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &MPI_COMM_CART);
		
		int coord[2];		// coordinates of a proc
		MPI_Cart_coords(MPI_COMM_CART, my_id, 2, coord);
		my_hor_pcoord = coord[0];
		my_vert_pcoord = coord[1];
		
		
		MPI_Comm_split(MPI_COMM_WORLD, my_vert_pcoord, 0, &MPI_COMM_VERT_SEGMENT);
		MPI_Comm_rank(MPI_COMM_VERT_SEGMENT, &my_vert_pcoord);
		
		MPI_Comm_split(MPI_COMM_WORLD, my_hor_pcoord, 0, &MPI_COMM_HOR_SEGMENT);
		MPI_Comm_rank(MPI_COMM_HOR_SEGMENT, &my_hor_pcoord);
		
		cout << "SpectralTransform::Init - coord: " << my_id << " " << my_vert_pcoord << " " << my_hor_pcoord << endl;
	}
	

    /*if (basis=="FFFW" && docomposition=="SLAB")
        Init_FFFW_SLAB();*/
    if (basis=="FFF" && decomposition=="SLAB")
        Init_FFF_SLAB();
	else if (basis=="SFF" && decomposition=="SLAB")
        Init_SFF_SLAB();
	else if (basis=="SSF" && decomposition=="SLAB")
        Init_SSF_SLAB();
	else if (basis=="SFF" && decomposition=="SLAB")
        Init_SFF_SLAB();
	else if (basis=="SSS" && decomposition=="SLAB")
        Init_SSS_SLAB();
    else if (basis=="ChFF" && decomposition=="SLAB")
        Init_ChFF_SLAB();
	
	else if (basis=="FFF" && decomposition=="PENCIL")
        Init_FFF_PENCIL();
	else if (basis=="SFF" && decomposition=="PENCIL")
        Init_SFF_PENCIL();
	else if (basis=="SSF" && decomposition=="PENCIL")
        Init_SSF_PENCIL();
	else if (basis=="SSS" && decomposition=="PENCIL")
        Init_SSS_PENCIL();
}


//*********************************************************************************************

// 3D: FFFW original
void SpectralTransform::FTr2c_xyz(Array<complx,3> Ar, Array<complx,3> A)
{
	//FFTW_EXECUTE_DFT_R2C_DP(r2c_xyz_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2r_xyz(Array<complx,3> A, Array<complx,3> Ar)
{
	//FFTW_EXECUTE_DFT_C2R_DP(c2r_xyz_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

// 2D FFFW 
void SpectralTransform::FTr2c_xz(Array<complx,2> Ar, Array<complx,2> A)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_xz_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2r_xz(Array<complx,2> A, Array<complx,2> Ar)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_xz_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<DP*>(Ar.data()));
}

// GP
// 3D: FFT original
void SpectralTransform::FTc2c_xyz(Array<complx,3> A)
{
	FFTW_MPI_EXECUTE_DFT_DP(c2c_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::IFTc2c_xyz(Array<complx,3> A)
{
	FFTW_MPI_EXECUTE_DFT_DP(c2c_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}


//********* FFFW, CFFF closed


void SpectralTransform::FTr2c_yz(Array<complx,3> X)
{
    fftw_execute_dft_r2c(r2c_yz_plan,reinterpret_cast<DP*>(X.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()));
}



void SpectralTransform::FTc2r_yz(Array<complx,3> X)
{
    fftw_execute_dft_c2r(c2r_yz_plan,reinterpret_cast<FFTW_COMPLEX_DP*>(X.data()), reinterpret_cast<DP*>(X.data()));
}


void SpectralTransform::FTr2c_yz(Array<DP,3> Xr)
{
    fftw_execute_dft_r2c(r2c_yz_plan,reinterpret_cast<DP*>(Xr.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Xr.data()));
}



void SpectralTransform::FTc2r_yz(Array<DP,3> Xr)
{
    fftw_execute_dft_c2r(c2r_yz_plan,reinterpret_cast<FFTW_COMPLEX_DP*>(Xr.data()), reinterpret_cast<DP*>(Xr.data()));
}


// 1D
// c2c
void SpectralTransform::FTc2c_x(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_DP(c2c_x_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}


void SpectralTransform::IFTc2c_x(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_DP(c2c_x_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2c_y(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_DP(c2c_y_forward_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}


void SpectralTransform::IFTc2c_y(Array<complx,3> A)
{
	FFTW_EXECUTE_DFT_DP(c2c_y_inverse_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}



// r2c, c2c FFF_PENCIL, SFF_SLAB
void SpectralTransform::FTr2c_z(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2r_z(Array<complx,2> A)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_z_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

void SpectralTransform::FTr2c_z(Array<DP,2> A)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
}

void SpectralTransform::FTc2r_z(Array<DP,2> A)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_z_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

void SpectralTransform::FTr2c_z(Array<complx,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		FFTW_EXECUTE_DFT_R2C_DP(r2c_z_plan, reinterpret_cast<DP*>(A(ly,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A(ly,Range::all(),Range::all()).data()));
}

void SpectralTransform::FTc2r_z(Array<complx,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		FFTW_EXECUTE_DFT_C2R_DP(c2r_z_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A(ly,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(A(ly,Range::all(),Range::all()).data()));
}

void SpectralTransform::FTr2c_z(Array<DP,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		FFTW_EXECUTE_DFT_R2C_DP(r2c_z_plan, reinterpret_cast<DP*>(A(ly,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_COMPLEX_DP*>(A(ly,Range::all(),Range::all()).data()));
}

void SpectralTransform::FTc2r_z(Array<DP,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		FFTW_EXECUTE_DFT_C2R_DP(c2r_z_plan, reinterpret_cast<FFTW_COMPLEX_DP*>(A(ly,Range::all(),Range::all()).data()), reinterpret_cast<DP*>(A(ly,Range::all(),Range::all()).data()));
}


// sin/cos
// along x (3D)

void SpectralTransform::SinCostr_x(char sincostr_switch, Array<DP,3> Ar)
{
	if (sincostr_switch == 'S') {
		FFTW_EXECUTE_R2R_DP(sintr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(costr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

void SpectralTransform::ISinCostr_x(char sincostr_switch, Array<DP,3> Ar)
{
    if (sincostr_switch == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R_DP(isintr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(icostr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}



// along y (3D)
void SpectralTransform::SinCostr_y(char sincostr_switch, Array<complx,3> A)
{
	if (sincostr_switch == 'S') {
	    FFTW_EXECUTE_R2R_DP(sintr_y_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
		ArrayShiftRight(A, 'Y');
	}
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(costr_y_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}


void SpectralTransform::ISinCostr_y(char sincostr_switch, Array<complx,3> A)
{
	if (sincostr_switch == 'S') {
		ArrayShiftLeft(A, 'Y');
	    FFTW_EXECUTE_R2R_DP(isintr_y_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
	}
	else if (sincostr_switch == 'C') 
		FFTW_EXECUTE_R2R_DP(icostr_y_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}


// Chebyshev transform via fftw cosine transform
void SpectralTransform::Chebyshevtr_x(Array<DP,3> A)
{
    FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
    // Postprocess (see definition)
    A(Range::all(),Range::all(),Range(1,Nx-2)) *= 2.0;
    
}

// Inverse Chebyshev transform via fftw cosine transform
void SpectralTransform::IChebyshevtr_x(Array<DP,3> A)
{
	// Preprocess (see definition)
    A(Range::all(),Range::all(),Range(1,Nx-2)) /= 2.0;
    
	FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

// transform along z
void SpectralTransform::SinCostr_z(char sincostr_switch, Array<complx,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		SinCostr_z(sincostr_switch, A(ly,Range::all(),Range::all()));
}

void SpectralTransform::ISinCostr_z(char sincostr_switch, Array<complx,3> A)
{
	for (int ly=0; ly<A.extent(0); ly++)
		ISinCostr_z(sincostr_switch, A(ly,Range::all(),Range::all()));
}

//****************** 2D 

// Used in SFF, SSF, and SSS
void SpectralTransform::SinCostr_x(char sincostr_switch, Array<DP,2> Ar)
{
    if (sincostr_switch == 'S') {
		FFTW_EXECUTE_R2R_DP(sintr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
		ArrayShiftRight(Ar, 'X');
	}
	
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(costr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}

// SpectralTransform::inverse
void SpectralTransform::ISinCostr_x(char sincostr_switch, Array<DP,2> Ar)
{
    if (sincostr_switch == 'S') {
		ArrayShiftLeft(Ar, 'X');
		FFTW_EXECUTE_R2R_DP(isintr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
	}
	
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(icostr_x_plan, reinterpret_cast<DP*>(Ar.data()), reinterpret_cast<DP*>(Ar.data()));
}




// Chebyshev transform via fftw cosine transform
// Used in ChFF
void SpectralTransform::Chebyshevtr_x(Array<DP,2> A)
{
    FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
    // Postprocess (see definition)
    A(Range::all(),Range(1,Nx-2)) *= 2.0;
    
}

// Inverse Chebyshev transform via fftw cosine transform
void SpectralTransform::IChebyshevtr_x(Array<DP,2> A)
{
	// Preprocess (see definition)
    A(Range::all(),Range(1,Nx-2)) /= 2.0;
    
	FFTW_EXECUTE_R2R_DP(Chebyshevtr_x_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

// transform along z
void SpectralTransform::SinCostr_z(char sincostr_switch, Array<complx,2> A)
{
	if (sincostr_switch == 'S') {
		FFTW_EXECUTE_R2R_DP(sintr_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
		ArrayShiftRight(A, 'Z');
	}
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(costr_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}

void SpectralTransform::ISinCostr_z(char sincostr_switch, Array<complx,2> A)
{
	if (sincostr_switch == 'S') {
		ArrayShiftLeft(A, 'Z');
		FFTW_EXECUTE_R2R_DP(isintr_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
	}
	
	else if (sincostr_switch == 'C')
		FFTW_EXECUTE_R2R_DP(icostr_z_plan, reinterpret_cast<DP*>(A.data()), reinterpret_cast<DP*>(A.data()));
}




//****************************************************************************************
// After sin transform of A, shift right along the shift direction.
// A(Ny/Nx, Nx/Ny, Nz/2+1) except for SSS basis for which Ar(Ny, Nx, Nz/2)
// Shift direction != Xa
void SpectralTransform::ArrayShiftRight(Array<complx,3> B, char shift_dirn)
{
	Array<DP,3> A(reinterpret_cast<DP*>(B.data()), B.shape()*shape(1,2,1), neverDeleteData);
	
    if (shift_dirn == 'Y') {
		A(Range(Ny-1,1,-1),Range::all(),Range::all()) = A(Range(Ny-2,0,-1),Range::all(),Range::all());
        A(0,Range::all(),Range::all()) = 0.0;
    }
	
	else if (shift_dirn == 'Z') {   // A(:,Nz,:)=A(:,Nz+1,:)=0 by definition (zero pad)
		A(Range::all(),Range(Nz-1,1,-1),Range::all()) = A(Range::all(),Range(Nz-2,0,-1),Range::all());
        A(Range::all(),0,Range::all()) = 0.0;
    }
	
	else
		cerr<< "In void SpectralTransform::ArrayShiftRight_SLAB(Array<complx,3> B, char shift_dirn): shift_dirn must be along Y or Z " << endl;
}


void SpectralTransform::ArrayShiftLeft(Array<complx,3> B, char shift_dirn)
{
	Array<DP,3> A(reinterpret_cast<DP*>(B.data()), B.shape()*shape(1,2,1), neverDeleteData);
	
	if (shift_dirn == 'Y') {
		A(Range(0,Ny-2),Range::all(),Range::all()) = A(Range(1,Ny-1),Range::all(),Range::all());
        A(Ny-1,Range::all(),Range::all()) = 0.0;
    }
	
	else if (shift_dirn == 'Z') {   // A(:,Nz,:)=A(:,Nz+1,:)=0 by definition (zero pad)
		A(Range::all(),Range(0,Nz-2),Range::all()) = A(Range::all(),Range(1,Nz-1),Range::all());
        A(Range::all(),Nz-1,Range::all()) = 0.0;
    }
	
	else 
		cerr<< "In void SpectralTransform::ArrayShiftLeft_SLAB(Array<complx,3> B, char shift_dirn): shift_dirn must be along Y or Z" << endl;
}

//****************************************************************************************

void SpectralTransform::ArrayShiftRight(Array<DP,3> A, char shift_dirn)
{	
	if (shift_dirn == 'X') {
		A(Range::all(),Range::all(),Range(Nx-1,1,-1)) = A(Range::all(),Range::all(),Range(Nx-2,0,-1));
        A(Range::all(),Range::all(),0) = 0.0;
    }
	
    else 
		cerr<< "In void SpectralTransform::ArrayShiftRight_SLAB(Array<DP,3> A, char shift_dirn): shift_dirn must be along X" << endl;
}


// Before sin transform, shift left
void SpectralTransform::ArrayShiftLeft(Array<DP,3> A, char shift_dirn)
{	
	if (shift_dirn == 'X') {
		A(Range::all(),Range::all(),Range(0,Nx-2)) = A(Range::all(),Range::all(),Range(1,Nx-1));
        A(Range::all(),Range::all(),Nx-1) = 0.0;
    }
	
    else 
		cerr<< "In void SpectralTransform::ArrayShiftLeft_SLAB(Array<DP,3> A, char shift_dirn): shift_dirn must be along X" << endl;
}



//****************************************************************************************

// 2D
// Shift direction = second direction
void SpectralTransform::ArrayShiftRight(Array<DP,2> A, char shift_dirn)
{
	
	if (shift_dirn == 'X') {
		A(Range::all(),Range(Nx-1,1,-1)) = A(Range::all(),Range(Nx-2,0,-1));
        A(Range::all(),0) = 0.0;
    }
	
	else
		cerr<< "In void SpectralTransform::ArrayShiftRight_SLAB(Array<DP,2> A, char shift_dirn): shift_dirn must be along X" << endl;
}

// 2D: shift_dirn = second dirn
void SpectralTransform::ArrayShiftLeft(Array<DP,2> A, char shift_dirn)
{
	
	if (shift_dirn == 'X') {
		A(Range::all(),Range(0,Nx-2)) = A(Range::all(),Range(1,Nx-1));
        A(Range::all(),Nx-1) = 0.0;
    }
	
	else
		cerr<< "In void SpectralTransform::ArrayShiftLeft_SLAB(Array<DP,2> A, char shift_dirn): shift_dirn must be along X" << endl;
}


// 2D
// Shift direction = second direction
void SpectralTransform::ArrayShiftRight(Array<complx,2> B, char shift_dirn)
{
	
	Array<DP,2> A(reinterpret_cast<DP*>(B.data()), B.shape()*shape(2,1), neverDeleteData);
	
	if (shift_dirn == 'Z') {   // A(:,Nz,:)=A(:,Nz+1,:)=0 by definition (zero pad)
		A(Range(Nz-1,1,-1),Range::all()) = A(Range(Nz-2,0,-1),Range::all());
        A(0,Range::all()) = 0.0;
	}
	
	else
		cerr << "In void SpectralTransform::ArrayShiftRight_SLAB(Array<complx,2> B, char shift_dirn): shift_dirn must be along Z. " << endl;
}


// Shift direction = second direction
void SpectralTransform::ArrayShiftLeft(Array<complx,2> B, char shift_dirn)
{
	
	Array<DP,2> A(reinterpret_cast<DP*>(B.data()), B.shape()*shape(2,1), neverDeleteData);
	
	if (shift_dirn == 'Z') {   // A(:,Nz,:)=A(:,Nz+1,:)=0 by definition (zero pad)
		A(Range(0,Nz-2),Range::all()) = A(Range(1,Nz-1),Range::all());
        A(Nz-1,Range::all()) = 0.0;
	}
	
	else
		cerr << "In void SpectralTransform::ArrayShiftLeft_SLAB(Array<complx,2> B, char shift_dirn): shift_dirn must be along Z. " << endl;
}




//********************************	End of four_tr.cc *****************************************


