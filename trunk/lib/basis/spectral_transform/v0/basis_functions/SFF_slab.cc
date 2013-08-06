/* Tarang-2
 *
 * Copyright_PENCIL( (C) 2008, 2009  Mahendra K. Verma
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

/*! \file scft_tr.cc 
 * 
 * @sa scft_tr.h
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "spectral_transform.h"



//*********************************************************************************************
void SpectralTransform::Init_SFF_SLAB()
{
	if (Ny > 1) {
		local_Nx = Nx/numprocs;
		local_Ny = Ny/numprocs;
		local_Nz = Nz/2+1;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = my_id * local_Ny;
		local_Nz_start = 0;
	}
	
	else if (Ny == 1) {
		local_Nx = Nx/numprocs;
		local_Ny = 1;
		local_Nz = (Nz/2+1)/numprocs;
		
		local_Nx_start = my_id * local_Nx;
		local_Ny_start = 0;
		local_Nz_start = my_id * local_Nz;
	}
    
    X.resize(local_Nx, Ny, Nz/2+1);
    Xr.resize(local_Ny,Nx,Nz/2+1);
    
    col_x_real.resize(Nx);
	plane_yz.resize(Ny,Nz/2+1);
	
	slab_transpose.AxyTR.resize(Ny,local_Nx);
	slab_transpose.AxyTR_recv.resize(Ny,local_Nx);
	
	slab_transpose.AyxTR.resize(Nx,local_Ny);
	slab_transpose.AyxTR_recv.resize(Nx,local_Ny);
	

    //Initializing plans
	
	sintr_x_plan = FFTW_PLAN_R2R_1D_DP(Nx,(DP*)(col_x_real.data()), (DP*)(col_x_real.data()), FFTW_RODFT10, FFTW_PLAN_FLAG);
	
	costr_x_plan = FFTW_PLAN_R2R_1D_DP(Nx,(DP*)(col_x_real.data()), (DP*)(col_x_real.data()), FFTW_REDFT10, FFTW_PLAN_FLAG);
	
	isintr_x_plan = FFTW_PLAN_R2R_1D_DP(Nx,(DP*)(col_x_real.data()), (DP*)(col_x_real.data()), FFTW_RODFT01, FFTW_PLAN_FLAG);
	
	icostr_x_plan = FFTW_PLAN_R2R_1D_DP(Nx,(DP*)(col_x_real.data()), (DP*)(col_x_real.data()), FFTW_REDFT01, FFTW_PLAN_FLAG);
	
	r2c_yz_plan = FFTW_PLAN_DFT_R2C_2D_DP(Ny, Nz, reinterpret_cast<DP*>(plane_yz.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(plane_yz.data()), FFTW_PLAN_FLAG);
	
	c2r_yz_plan = FFTW_PLAN_DFT_C2R_2D_DP(Ny, Nz, reinterpret_cast<FFTW_COMPLEX_DP*>(plane_yz.data()), reinterpret_cast<DP*>(plane_yz.data()),FFTW_PLAN_FLAG);
    
    
    // For 2D
    if (Ny==1) {
        col_z.resize(Nz/2+1);
        
        r2c_z_plan = FFTW_PLAN_DFT_R2C_1D_DP(Nz,reinterpret_cast<DP*>(col_z.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(col_z.data()),FFTW_PLAN_FLAG);
        
        c2r_z_plan = FFTW_PLAN_DFT_C2R_1D_DP(Nz, reinterpret_cast<FFTW_COMPLEX_DP*>(col_z.data()), reinterpret_cast<DP*>(col_z.data()),FFTW_PLAN_FLAG);
        
		slab_transpose.AxzTR.resize(Nz/2+1,local_Nx);
		slab_transpose.AxzTR_recv.resize(Nz/2+1,local_Nx);
        
        slab_transpose.AzxTR.resize(Nx,local_Nz);
        slab_transpose.AzxTR_recv.resize(Nx,local_Nz);
	}
	


}
//*********************************************************************************************


void SpectralTransform::Zero_pad_last_plane_SFF_SLAB(Array<complx,3> Ar)
{
    Ar(Range::all(),Range::all(),Nz/2) = 0.0;
}


void SpectralTransform::Zero_pad_last_plane_SFF_SLAB(Array<complx,2> Ar)
{
    if (my_id == numprocs-1)
        Ar(local_Nz-1, 0, Range::all()) = 0.0;
    // Zero pad Ar(Nz/2+1,0,:)
}



//*********************************************************************************************

void SpectralTransform::Norm_SFF_SLAB(Array<complx,3> A) 
{  
	A = A/(2*DP(Nx)*DP(Ny)*DP(Nz));
}


void SpectralTransform::Norm_SFF_SLAB(Array<complx,2> A)
{
	A = A/(2*DP(Nx)*DP(Nz));
}

/**********************************************************************************************

	SFT(Ar) = A; Ar(Ny, Nx, Nz/2+1) is in transposed order

	Note: z= Nz/2 plane does contain any real data (zero everywhere).
	
***********************************************************************************************/

void SpectralTransform::Forward_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,3> Ar, Array<complx,3> A)
{
        
        Zero_pad_last_plane_SFF_SLAB(Ar);							 // Zero_pad the row=Nz/2
  			
		// Array Ar is transposed N2, N1, N3
		// Sin transform along x for all y,z of Ar
		for (int ly=0; ly<local_Ny; ly++)
			for (int lz = 0; lz < Nz/2; lz++)   {
				col_x_real = real(Ar(ly, Range::all(), lz)); 
                
                if (sincostr_switch[0] == 'S')
                    Sintr_x(col_x_real);  
                else if (sincostr_switch[0] == 'C')
                    Costr_x(col_x_real); 
                
				real(Ar(ly, Range::all(), lz)) = col_x_real;
		
				col_x_real = imag(Ar(ly, Range::all(), lz)); 
				
                if (sincostr_switch[0] == 'S')
                    Sintr_x(col_x_real);  
                else if (sincostr_switch[0] == 'C')
                    Costr_x(col_x_real);
                
				imag(Ar(ly, Range::all(), lz)) = col_x_real;
			}
        
        if (sincostr_switch[0] == 'S')
            ArrayShiftRight_SLAB(Ar, 'Y', 'X', 'Z', 'X');	
	  
		Transpose_array_SLAB(Ar, A, 'Y', 'X', 'Z', 'Z');
	  
        // FT;  Array A is normal ordered
        // FFT along y-z planes of  A
        for (int lx = 0; lx < local_Nx; lx++)  {					
            plane_yz = A(lx, Range::all(), Range::all());
            FTr2c_yz(plane_yz);
            A(lx, Range::all(), Range::all()) = plane_yz;
        }
		
		Norm_SFF_SLAB(A);
}

void SpectralTransform::Forward_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,2> Ar, Array<complx,2> A)
{
		
		Zero_pad_last_plane_SFF_SLAB(Ar);
		
		// Sin transform along x for all z of Ar
		for (int lz=0; lz<local_Nz; lz++) {
			col_x_real = real(Ar(lz,Range::all())); 
            
            if (sincostr_switch[0] == 'S')
                Sintr_x(col_x_real);  
            else if (sincostr_switch[0] == 'C')
                Costr_x(col_x_real);
            
			real(Ar(lz,Range::all())) = col_x_real;
			
			col_x_real = imag(Ar(lz,Range::all())); 
            
            if (sincostr_switch[0] == 'S')
                Sintr_x(col_x_real);
            else if (sincostr_switch[0] == 'C')
                Costr_x(col_x_real);

			imag(Ar(lz,Range::all())) = col_x_real;
		}
        
        if (sincostr_switch[0] == 'S')
            ArrayShiftRight_SLAB(Ar, 'X');
		
		Transpose_array_SLAB(Ar, A, 'Z', 'X'); // zx->xz
		
		// Array A is normal ordered
		// FFT along z col of  A
		for (int lx=0; lx<local_Nx; lx++)  {
			col_z = A(lx,Range::all());
            FTr2c_z(col_z);
			A(lx,Range::all()) = col_z;
		}
		
		Norm_SFF_SLAB(A);  
}


//
//
void SpectralTransform::Forward_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,3> A)
{
    Transpose_array_SLAB(A, Xr, 'X', 'Y', 'Z', 'Z');
        
	Forward_transform_array_transpose_order_SFF_SLAB(sincostr_switch, Xr, A);
}

void SpectralTransform::Forward_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,2> A)
{
    Transpose_array_SLAB(A, Xr(Range::all(),0,Range::all()), 'X', 'Z'); // xz -> zx
    
	Forward_transform_array_transpose_order_SFF_SLAB(sincostr_switch, Xr(Range::all(),0,Range::all()), A);
}


/**********************************************************************************************

	Inverse SFT (A)  = Ar 
	IFT along perp dirn and SIN transform along x dirn of A 

***********************************************************************************************/

void SpectralTransform::Inverse_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,3> A, Array<complx,3> Ar)
{
    // Array A is normal  ordered
    // IFFT along y-z planes of  A
    for (int lx=0; lx < local_Nx; lx++) {
        plane_yz = A(lx, Range::all(), Range::all());
        FTc2r_yz(plane_yz);
        X(lx, Range::all(), Range::all()) = plane_yz;
    }
    

    Transpose_array_SLAB(X, Ar, 'X', 'Y', 'Z', 'Z');
    
   if (sincostr_switch[0] =='S')
        ArrayShiftLeft_SLAB(Ar, 'Y', 'X', 'Z', 'X');
    
    // Array Ar is transposed N2, N1, N3
    // ISin transform along x for all y,z of Ar
    for (int ly=0; ly < local_Ny; ly++)											
         for (int lz=0; lz<Nz/2; lz++) {        
            col_x_real = real(Ar(ly, Range::all(), lz)); 
             
            if (sincostr_switch[0] == 'S') 
                ISintr_x(col_x_real);
            else if (sincostr_switch[0] == 'C')
                ICostr_x(col_x_real);
                 
            real(Ar(ly, Range::all(), lz)) = col_x_real ;
    
            col_x_real = imag(Ar(ly, Range::all(), lz));
             
            if (sincostr_switch[0] == 'S')
                ISintr_x(col_x_real);
            else if (sincostr_switch[0] == 'C')
                ICostr_x(col_x_real);

            imag(Ar(ly, Range::all(), lz)) = col_x_real ;
        }
    Zero_pad_last_plane_SFF_SLAB(Ar);	
}


void SpectralTransform::Inverse_transform_array_transpose_order_SFF_SLAB(string sincostr_switch, Array<complx,2> A, Array<complx,2> Ar)
{
		
    // Array A is normal  ordered
    // IFFT along z col of  A
    for (int lx=0; lx < local_Nx; lx++) {
        col_z = A(lx,Range::all());
        FTc2r_z(col_z);
        X(lx,0,Range::all()) = col_z;
    }

    Transpose_array_SLAB(X(Range::all(),0,Range::all()), Ar, 'X', 'Z'); // xz -> zx

    if (sincostr_switch[0] =='S')
        ArrayShiftLeft_SLAB(Ar, 'X');
    // Array Ar is transposed N3/2+1, 0, local_Nx
    // ISin transform along x for all z of Ar
    for (int lz=0; lz < local_Nz; lz++)	{        
        col_x_real = real(Ar(lz,Range::all()));
        
        if (sincostr_switch[0] == 'S')
            ISintr_x(col_x_real);
        else if (sincostr_switch[0] == 'C')
            ICostr_x(col_x_real);
        
        real(Ar(lz,Range::all())) = col_x_real;
        
        col_x_real = imag(Ar(lz,Range::all()));
        
        if (sincostr_switch[0] == 'S')
            ISintr_x(col_x_real);
        else if (sincostr_switch[0] == 'C')
            ICostr_x(col_x_real);
        
        imag(Ar(lz,Range::all())) = col_x_real;
    }

    Zero_pad_last_plane_SFF_SLAB(Ar); 
}

//
//
void SpectralTransform::Inverse_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,3> A)
{
	Inverse_transform_array_transpose_order_SFF_SLAB(sincostr_switch, A, Xr);
    
    Transpose_array_SLAB(Xr, A, 'Y', 'X', 'Z', 'Z');
}


void SpectralTransform::Inverse_transform_array_SFF_SLAB(string sincostr_switch, Array<complx,2> A)
{
	Inverse_transform_array_transpose_order_SFF_SLAB(sincostr_switch, A, Xr(Range::all(),0,Range::all()));
    
    Transpose_array_SLAB(Xr(Range::all(),0,Range::all()), A, 'Z', 'X');  //zx->xz
}

//******************************** End of SFF_slab_tr.cc  *****************************************
