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


/*! \file  RSprod.cc
 * FluidVF& U
 * @brief  Compute the diagonal and nondiagonal terms of the real space products 
 *			(NORMAL ORDER).
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug	No known bugs
 */


#include "Nlin.h"


//******************************************************************************


void Nlin_incompress::Compute_nlin_diag(FluidVF& U)
{
		
    global.program.sincostr_switch = global.program.sincostr_switch_Visqr;
    ArrayOps::Real_space_multiply(U.rvf.V1r, U.rvf.V1r, global.temp_array.Xr);
	
	universal->Forward_transform(global.temp_array.Xr, U.nlin1);
    universal->Xderiv(U.nlin1, U.nlin1);

	if (!global.program.two_dimension) {  // for 3d and 2.5d
        ArrayOps::Real_space_multiply(U.rvf.V2r, U.rvf.V2r, global.temp_array.Xr);
        universal->Forward_transform(global.temp_array.Xr, U.nlin2);
		universal->Yderiv(U.nlin2, U.nlin2);
    } 
    
	ArrayOps::Real_space_multiply(U.rvf.V3r, U.rvf.V3r, global.temp_array.Xr);
    universal->Forward_transform(global.temp_array.Xr, U.nlin3);
    universal->Zderiv(U.nlin3, U.nlin3);
}


void Nlin_incompress::Compute_nlin_offdiag(FluidVF& U)
{
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		
		
		global.temp_array.Xr = U.rvf.V1r;
        ArrayOps::Real_space_multiply(U.rvf.V1r, U.rvf.V2r, U.rvf.V1r);			// V1r= V1rV2r
        ArrayOps::Real_space_multiply(U.rvf.V2r, U.rvf.V3r, U.rvf.V2r);			// V2r= V2rV3r
		ArrayOps::Real_space_multiply(U.rvf.V3r, global.temp_array.Xr, U.rvf.V3r);			// V3r= V3rV1r
		
        // Forward transform(V1rV2r) and derivatives.
        global.program.sincostr_switch = global.program.sincostr_switch_VxVy;
        universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
		
		universal->Add_Yderiv(global.temp_array.X, U.nlin1);
		universal->Add_Xderiv(global.temp_array.X, U.nlin2);
        
        // Forward transform(V2rV3r) and derivatives.
        global.program.sincostr_switch = global.program.sincostr_switch_VyVz;
        universal->Forward_transform(U.rvf.V2r, global.temp_array.X);	
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin2);
		universal->Add_Yderiv(global.temp_array.X, U.nlin3);
        
        // Forward transform(V1rV3r) and derivatives.
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);
	}
	
	// 2D...
	// *****
	else  {
        ArrayOps::Real_space_multiply(U.rvf.V3r, U.rvf.V1r, U.rvf.V3r);			// V3r= V3rV1r 
        
        // Forward transform(V1rV3r) and derivatives.
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);	
        
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);
	}
}


	// For T.nlin only
void Nlin_incompress::Compute_nlin_VxT(FluidVF& U, FluidSF& T)
{
    global.program.sincostr_switch = global.program.sincostr_switch_FVx;
    ArrayOps::Real_space_multiply(U.rvf.V1r, T.rsf.Fr, global.temp_array.Xr);
    universal->Forward_transform(global.temp_array.Xr, T.nlin);
    universal->Xderiv(T.nlin, T.nlin);
    
    if (!global.program.two_dimension) {  // for 3d and 2.5d
        global.program.sincostr_switch = global.program.sincostr_switch_FVy;
        ArrayOps::Real_space_multiply(U.rvf.V2r, T.rsf.Fr, global.temp_array.Xr);           
        universal->Forward_transform(global.temp_array.Xr, global.temp_array.X);
        
        universal->Add_Yderiv(global.temp_array.X, T.nlin);
    }
    
    global.program.sincostr_switch = global.program.sincostr_switch_FVz;
    ArrayOps::Real_space_multiply(U.rvf.V3r, T.rsf.Fr, global.temp_array.Xr);
    universal->Forward_transform(global.temp_array.Xr, global.temp_array.X);
    
    universal->Add_Zderiv(global.temp_array.X, T.nlin);
}


void Nlin_incompress::Compute_nlin_diag(FluidVF& U, FluidVF& W)           
{
	if (global.program.kind != "KEPLERIAN") {
        global.program.sincostr_switch = global.program.sincostr_switch_Visqr;
        
        ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V1r, global.temp_array.Xr);
        universal->Forward_transform(global.temp_array.Xr, U.nlin1);
        universal->Xderiv(U.nlin1, U.nlin1);  // parity unimportant here.
        
        if (!global.program.two_dimension) {  // for 3d and 2.5d
            ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V2r, global.temp_array.Xr);
            universal->Forward_transform(global.temp_array.Xr, U.nlin2);
            universal->Yderiv(U.nlin2, U.nlin2);  // parity unimportant here.
        }
        
        ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V3r, global.temp_array.Xr); 
        universal->Forward_transform(global.temp_array.Xr, U.nlin3);
        universal->Zderiv(U.nlin3, U.nlin3);	// parity unimportant here.
        
        W.nlin1 = U.nlin1;  
        
        if (!global.program.two_dimension)   // for 3d and 2.5d
            W.nlin2 = U.nlin2;	
        
        W.nlin3 = U.nlin3;
    }
    
    else {
     /*   // for Keplerian: (Basis=FFF with no-transpose option)
        // d/dx -> d/dx + q omega t d/dy
        Real omega_keplerian = global.force.double_para(0);
        Real q_keplerian = global.force.double_para(1);
        
        Real q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
        
        global.program.sincostr_switch = global.program.sincostr_switch_Visqr;
        
        ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V1r, U.nlin1);
        universal->Forward_transform_array(U.nlin1);
        
        global.temp_array.X = U.nlin1;
        universal->Yderiv(global.temp_array.X,  global.temp_array.X);
        
        universal->Xderiv(U.nlin1, U.nlin1);
        U.nlin1 = U.nlin1 + complex<Real>(q_omega_t,0)*global.temp_array.X;
        
        if (!global.program.two_dimension) {  // for 3d and 2.5d
            ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V2r, U.nlin2);
            universal->Forward_transform_array(U.nlin2);
            universal->Yderiv(U.nlin2, U.nlin2);
        }
        
        ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V3r, U.nlin3);
        universal->Forward_transform_array(U.nlin3);
        universal->Zderiv(U.nlin3, U.nlin3);
        
        W.nlin1 = U.nlin1;
        W.nlin2 = U.nlin2;  // MRI is 2.5d
        W.nlin3 = U.nlin3; */
    }
}


void Nlin_incompress::Compute_nlin_offdiag(FluidVF& U, FluidVF& W)
{
    if (!global.program.two_dimension) {  // for 3d and 2.5d
        
        global.temp_array.Xr = U.rvf.V1r;
        global.temp_array.Xr2 = W.rvf.V1r;							// store V1r and W1r in X,X2
        
        ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V2r, U.rvf.V1r);	// V1r = V1r*W2r
        ArrayOps::Real_space_multiply(W.rvf.V1r, U.rvf.V2r, W.rvf.V1r);	// W1r = W1r*V2r
        
        ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V3r, U.rvf.V2r);	// V2r = V2r*W3r
        ArrayOps::Real_space_multiply(W.rvf.V2r, U.rvf.V3r, W.rvf.V2r);	// W2r = W2r*V3r
        
        ArrayOps::Real_space_multiply(U.rvf.V3r, global.temp_array.Xr2, U.rvf.V3r);	// V3r = V3r*W1r
        ArrayOps::Real_space_multiply(W.rvf.V3r, global.temp_array.Xr, W.rvf.V3r);	// W3r = W3r*V1r
		
        // Transform and derivative
        global.program.sincostr_switch = global.program.sincostr_switch_VxVy;
        universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
		
		universal->Add_Yderiv(global.temp_array.X, U.nlin1);	// N1 += D2(W2r* V1r)
		universal->Add_Xderiv(global.temp_array.X, W.nlin2);	// W.N2 += D1(V1r*W2r)
	
        
        universal->Forward_transform(W.rvf.V1r, global.temp_array.X);
		
		universal->Add_Yderiv(global.temp_array.X, W.nlin1);	// W.N1 += D2(V2r* W1r)
		universal->Add_Xderiv(global.temp_array.X, U.nlin2);	// N2 += D1(W1r*V2r)
		
        
        
        global.program.sincostr_switch = global.program.sincostr_switch_VyVz;
        universal->Forward_transform(U.rvf.V2r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin2);	// N2 += D3(W3r*V2r)
		universal->Add_Yderiv(global.temp_array.X, W.nlin3);	// W.N3 += D2(W3r*V2r)
		
		
        universal->Forward_transform(W.rvf.V2r, global.temp_array.X);
		
		
		universal->Add_Zderiv(global.temp_array.X, W.nlin2);	//W.N2 += D3(V3r*W2r)
		universal->Add_Yderiv(global.temp_array.X, U.nlin3);	//N3 += D2(V3r*W2r)
        
        
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, W.nlin1);	// W.N1 += D3(V3r*W1r)
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);	// N3 += D1(V3r*W1r)
        
        universal->Forward_transform(W.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);	// N1 += D3(W3r*V1r)
		universal->Add_Xderiv(global.temp_array.X, W.nlin3);	// W.N3 += D1(W3r*V1r)
	
        
        // For Keplerian: Notranspose order, FFF basis only
        if (global.program.kind == "KEPLERIAN") {
       /*     Real omega_keplerian = global.force.double_para(0);
            Real q_keplerian = global.force.double_para(1);
            
            Real q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
            
            universal->Yderiv(W.rvf.V1r, global.temp_array.X);
            U.nlin2 = U.nlin2 + Complex(q_omega_t,0)*global.temp_array.X;
            
            universal->Yderiv(U.rvf.V3r, global.temp_array.X);
            U.nlin3 = U.nlin3 + Complex(q_omega_t,0)*global.temp_array.X;
            
            universal->Yderiv(U.rvf.V1r, global.temp_array.X);
            W.nlin2 = W.nlin2 + Complex(q_omega_t,0)*global.temp_array.X;
            
            universal->Yderiv(W.rvf.V3r, global.temp_array.X);
            W.nlin3 = W.nlin3 + Complex(q_omega_t,0)*global.temp_array.X; */
        }
    }
    
    // 2D...
    else {
        ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V1r, U.rvf.V3r);	// V3r = V3r*W1r
        ArrayOps::Real_space_multiply(W.rvf.V3r, U.rvf.V1r, W.rvf.V3r);	// W3r = W3r*V1r
        
        // Transform and derivative
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, W.nlin1);	// W.N1 += D3(V3r*W1r)
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);	// N3 += D1(V3r*W1r)

        
        universal->Forward_transform(W.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);	// N1 += D3(W3r*V1r)
		universal->Add_Xderiv(global.temp_array.X, W.nlin3);	// W.N3 += D1(W3r*V1r)
    }
}

//****************************  FOR ENERGY TRANSPOSE.cc **********************************************

// Treat W.V1 as T.F (scalar)
// U.nlin1 = FT(Dj (Uj W1))
void Nlin_incompress::Compute_nlin_first_component_product(FluidVF& U, PlainFluidVF& W)
{
    global.program.sincostr_switch = global.program.sincostr_switch_FVx;
    ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V1r, global.temp_array.Xr);
    universal->Forward_transform(global.temp_array.Xr, U.nlin1);
    universal->Xderiv(U.nlin1, U.nlin1);	
    
    if (!global.program.two_dimension) {  // for 3d and 2.5d
        global.program.sincostr_switch = global.program.sincostr_switch_FVy;
        ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V1r, global.temp_array.Xr);
        universal->Forward_transform(global.temp_array.Xr, global.temp_array.X);
        universal->Yderiv(global.temp_array.X, global.temp_array.X);	
        U.nlin1 = U.nlin1 + global.temp_array.X;
    }
    
    global.program.sincostr_switch = global.program.sincostr_switch_FVz;
    ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V1r, global.temp_array.Xr);			
    universal->Forward_transform(global.temp_array.Xr, global.temp_array.X);
    universal->Zderiv(global.temp_array.X, global.temp_array.X);	
    U.nlin1 = U.nlin1 + global.temp_array.X;
	
}


// U.nlin_i = FT(Di (U_i W_i))
void Nlin_incompress::Compute_nlin_diag(FluidVF& U, PlainFluidVF& W)           
{
   
    global.program.sincostr_switch = global.program.sincostr_switch_Visqr;
    
    ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V1r, global.temp_array.Xr);
	
    universal->Forward_transform(global.temp_array.Xr, U.nlin1);
    universal->Xderiv(U.nlin1, U.nlin1);  // parity unimportant here.
    
    if (!global.program.two_dimension) {  // for 3d and 2.5d
        ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V2r, global.temp_array.Xr);
        universal->Forward_transform(global.temp_array.Xr, U.nlin2);
        universal->Yderiv(U.nlin2, U.nlin2);  // parity unimportant here.
    }	
    
    ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V3r, global.temp_array.Xr); 
    universal->Forward_transform(global.temp_array.Xr, U.nlin3);
	universal->Zderiv(U.nlin3, U.nlin3);	// parity unimportant here.
}


// U.nlin_i = FT(Dj (U_j W_i)) 
void Nlin_incompress::Compute_nlin_offdiag(FluidVF& U, PlainFluidVF& W)
{
	if (!global.program.two_dimension) {  // for 3d and 2.5d
			
        global.temp_array.Xr = U.rvf.V1r;
        global.temp_array.Xr2 = W.rvf.V1r;							// store V1r and W1r in X,X2
        
        ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V2r, U.rvf.V1r);	// V1r = V1r*W2r
        ArrayOps::Real_space_multiply(W.rvf.V1r, U.rvf.V2r, W.rvf.V1r);	// W1r = W1r*V2r
        
        ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V3r, U.rvf.V2r);	// V2r = V2r*W3r
        ArrayOps::Real_space_multiply(W.rvf.V2r, U.rvf.V3r, W.rvf.V2r);	// W2r = W2r*V3r
        
        ArrayOps::Real_space_multiply(U.rvf.V3r, global.temp_array.Xr2, U.rvf.V3r);	// V3r = V3r*W1r
        ArrayOps::Real_space_multiply(W.rvf.V3r, global.temp_array.Xr, W.rvf.V3r);	// W3r = W3r*V1r	
        
            // Transform and derivative
        global.program.sincostr_switch = global.program.sincostr_switch_VxVy;
        universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
		universal->Add_Xderiv(global.temp_array.X, U.nlin2);	// N2 += D1(V1r*W2r)
	
        
        universal->Forward_transform(W.rvf.V1r, global.temp_array.X);		
		universal->Add_Yderiv(global.temp_array.X, U.nlin1);	// N1 += D2(V2r*W1r)
		
        
        
        global.program.sincostr_switch = global.program.sincostr_switch_VyVz;
        universal->Forward_transform(U.rvf.V2r, global.temp_array.X);
		universal->Add_Yderiv(global.temp_array.X, U.nlin3);	// N3 += D2(V2r*W3r)
	
        
        universal->Forward_transform(W.rvf.V2r, global.temp_array.X);
		universal->Add_Zderiv(global.temp_array.X, U.nlin2);	//N2 += D3(V3r*W2r)
		
        
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);	// N1 += D3(V3r*W1r)
		
        
        universal->Forward_transform(W.rvf.V3r, global.temp_array.X);
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);	// N3 += D1(V1r*W3r)
	}
	
	// 2D...
	// *****
	else  {
        ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V1r, U.rvf.V3r);	// V3r = V3r*W1r
        ArrayOps::Real_space_multiply(W.rvf.V3r, U.rvf.V1r, W.rvf.V3r);	// W3r = W3r*V1r	
        
        // Transform and derivative
        global.program.sincostr_switch = global.program.sincostr_switch_VxVz;
        universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
		
		universal->Add_Zderiv(global.temp_array.X, U.nlin1);	// N1 += D3(V3r*W1r)
        
        universal->Forward_transform(W.rvf.V3r, global.temp_array.X);
		
		universal->Add_Xderiv(global.temp_array.X, U.nlin3);	//N3 += D1(V1r*W3r)
	}
}


// GP
void Nlin_incompress::Compute_nlin_Tsqr(FluidSF& T)
{
    
/*    Real C=1.0;
    Real Potential_amp = 0.0;
    
    //	global.program.sincostr_switch = global.program.sincostr_switch_Visqr;
    
    global.temp_array.X = T.rsf.Fr;
    
    real(T.rsf.Fr) = sqr(abs(T.rsf.Fr));

    // Compute Vtrap
    int rx,ry,rz;
    Real x,y,z;

    
      for (int lx=0; lx<local_Nx; lx++) {
        rx = universal->Get_rx_real_space(lx);
        x = rx*global.field.L[1]/Nx;
        
        for (int ly=0; ly<Ny; ly++) {
            ry = universal->Get_ry_real_space(ly);
            y = ry*global.field.L[2]/Ny;
            
            for (int lz=0; lz<Nz; lz++) {
                universal->Get_rz_real_space(lz, rz, 0);
                z = rz*global.field.L[3]/Nz;
                
                real(T.nlin(lx,ly,lz)) += Potential_amp*((x-M_PI)*(x-M_PI)+ (y-M_PI)*(y-M_PI) + (z-M_PI)*(z-M_PI));
            }
        }
    } */
    
    // (|psi|^2+V)*(-i C psi)
 /*   real(T.rsf.Fr) = real(T.rsf.Fr) * (C*imag(global.temp_array.X));
    imag(T.rsf.Fr) = real(T.rsf.Fr) * (-C*real(global.temp_array.X));
    
    universal->Forward_transform(T.rsf.Fr, T.nlin); */
}



void Nlin_incompress::Compute_nlin_helical(FluidVF& U, PlainFluidVF& W)
{
  if (!global.program.two_dimension) {  // for 3d and 2.5d
    
    U.Inverse_transform();
    W.Inverse_transform();

    
    
    ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V3r, global.temp_array.Xr);	// V1r = V1r*W2r
    ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V2r, global.temp_array.Xr2);
    
    global.temp_array.Xr = global.temp_array.Xr-global.temp_array.Xr2;
    universal->Forward_transform(global.temp_array.Xr, U.nlin1);
    
    
    ArrayOps::Real_space_multiply(U.rvf.V3r, W.rvf.V1r, global.temp_array.Xr);
    ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V3r, global.temp_array.Xr2);
    
    global.temp_array.Xr = global.temp_array.Xr-global.temp_array.Xr2;
    universal->Forward_transform(global.temp_array.Xr, U.nlin2);
    
    ArrayOps::Real_space_multiply(U.rvf.V1r, W.rvf.V2r, global.temp_array.Xr);
    ArrayOps::Real_space_multiply(U.rvf.V2r, W.rvf.V1r, global.temp_array.Xr2);
    
    global.temp_array.Xr = global.temp_array.Xr-global.temp_array.Xr2;
    universal->Forward_transform(global.temp_array.Xr, U.nlin3);
    
  }
  
  // 2D...
  else {
    cout<<"2D will be implemented latter. "<<endl;

  }
}

//****************************  End of RSprod.cc **********************************************



