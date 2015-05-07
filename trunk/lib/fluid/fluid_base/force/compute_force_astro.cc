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
/*! \file  compute_force_RB.cc
 * 
 * @brief Force for RB convection
 *
 *		For Pr_large, u_small:		 F1 = Ra*Pr*T(k);  
 *		For Pr_large, u_large:		 F1 = Ra*T(k);    	
 *		For Pr_small or 0, u_small:  F1 = Ra*T(k);	   
 *		For Pr_small or 0, u_large:  F1 = Pr*T(k)
 *
 *		For Pr_large				Ftheta = V1;    	
 *		For Pr_small				Ftheta = V1/Pr;	   
 *		For Pr = 0					Ftheta = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 
#include "FORCE.h"


//*********************************************************************************************


void FORCE::Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	Real Schwarzschild = global.force.double_para(0);
	Real Peclet_astro = global.force.double_para(1);
	
	U.Force1 =  (global.PHYSICS.Rayleigh/global.PHYSICS.Chandrasekhar)*(T.csf.F);
	
	T.Force = (-Schwarzschild)*(U.cvf.V1);
	
	temp1.resize(shape_real_array);
	temp2.resize(shape_complex_array);
	
	// Compute unit vector(B)
	W.Inverse_transform();
	
	Real abs_Bxyz;
	for (int ly=0; ly<=(W.rvf.V1r).extent(0); ly++)
        for (int lz=0; lz<=(W.rvf.V1r).extent(1); lz++)
            for (int lx=0; lx<=(W.rvf.V1r).extent(2); lx++) {
				
				abs_Bxyz = sqrt(my_pow(W.rvf.V1r(ly,lz,lx),2)+my_pow(W.rvf.V2r(ly,lz,lx),2)+my_pow(W.rvf.V3r(ly,lz,lx),2));
				if (abs_Bxyz> MYEPS2) {
					W.rvf.V1r(ly,lz,lx) /= abs_Bxyz;
					W.rvf.V2r(ly,lz,lx) /= abs_Bxyz;
					W.rvf.V3r(ly,lz,lx) /= abs_Bxyz;
				}				
			}
	
	// COmpute Qtheta
	// Compute \hat B_j \delta_j \theta
	universal->Xderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	global.temp_array.Xr -= 1.0;  // To care of Tc
	ArrayOps::Real_space_multiply(W.rvf.V1r, global.temp_array.Xr, temp1);
	
	universal->Yderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V2r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	universal->Zderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V3r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	ArrayOps::Real_space_multiply(W.rvf.V1r, temp1, U.rvf.V1r);
	ArrayOps::Real_space_multiply(W.rvf.V2r, temp1, U.rvf.V2r);
	ArrayOps::Real_space_multiply(W.rvf.V3r, temp1, U.rvf.V3r);
	// U.rvf.Vir = Qtheta
	
	universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
	universal->Xderiv(global.temp_array.X, temp2);
	
	universal->Forward_transform(U.rvf.V2r, global.temp_array.X);
	universal->Yderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	
	universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
	universal->Zderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	// temp2 = div(Qtheta)
	
	T.Force = T.Force + (temp2)/(Peclet_astro);
}


void FORCE::Compute_force_astro(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	
	Real Schwarzschild = global.force.double_para(0);
	Real Peclet_astro = global.force.double_para(1);
	Real Peclet_c_astro = global.force.double_para(1);
	
	U.Force1 =  (global.PHYSICS.Rayleigh/global.PHYSICS.Chandrasekhar)*(T.csf.F);
	
	T.Force = (-Schwarzschild)*(U.cvf.V1);
	
	C.Force = (U.cvf.V1);
	
	temp1.resize(shape_real_array);
	temp2.resize(shape_complex_array);
	
	// Compute unit vector(B)
	W.Inverse_transform();
	Real abs_Bxyz;
	for (int ly=0; ly<=(W.rvf.V1r).extent(0); ly++)
        for (int lz=0; lz<=(W.rvf.V1r).extent(1); lz++)
            for (int lx=0; lx<=(W.rvf.V1r).extent(2); lx++) {
				
				abs_Bxyz = sqrt(my_pow(W.rvf.V1r(ly,lz,lx),2)+my_pow(W.rvf.V2r(ly,lz,lx),2)+my_pow(W.rvf.V3r(ly,lz,lx),2));
				if (abs_Bxyz> MYEPS2) {
					W.rvf.V1r(ly,lz,lx) /= abs_Bxyz;
					W.rvf.V2r(ly,lz,lx) /= abs_Bxyz;
					W.rvf.V3r(ly,lz,lx) /= abs_Bxyz;
				}
			}
	
	// Compute Qtheta
	// Compute \hat B_j \delta_j \theta
	universal->Xderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	global.temp_array.Xr -= 1.0;  // To care of Tc
	ArrayOps::Real_space_multiply(W.rvf.V1r, global.temp_array.Xr, temp1);
	
	universal->Yderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V2r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	universal->Zderiv(T.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V3r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	ArrayOps::Real_space_multiply(W.rvf.V1r, temp1, U.rvf.V1r);
	ArrayOps::Real_space_multiply(W.rvf.V2r, temp1, U.rvf.V2r);
	ArrayOps::Real_space_multiply(W.rvf.V3r, temp1, U.rvf.V3r);
	// U.rvf.Vir = Qtheta
	
	// Div calc
	universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
	universal->Xderiv(global.temp_array.X, temp2);
	
	universal->Forward_transform(U.rvf.V2r, global.temp_array.X);
	universal->Yderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	
	universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
	universal->Zderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	// temp2 = div(Qtheta)
	
	T.Force += temp2/(Peclet_astro);
	
	
	// Compute Qc
	// Compute \hat B_j \delta_j \theta
	universal->Xderiv(C.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	global.temp_array.Xr -= 1.0;  // To care of concentration gradient
	ArrayOps::Real_space_multiply(W.rvf.V1r, global.temp_array.Xr, temp1);
	
	universal->Yderiv(C.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V2r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	universal->Zderiv(C.csf.F, global.temp_array.X);
	universal->Inverse_transform(global.temp_array.X, global.temp_array.Xr);
	ArrayOps::Real_space_multiply(W.rvf.V3r, global.temp_array.Xr, global.temp_array.Xr2);
	temp1 += global.temp_array.Xr2;
	
	ArrayOps::Real_space_multiply(W.rvf.V1r, temp1, U.rvf.V1r);
	ArrayOps::Real_space_multiply(W.rvf.V2r, temp1, U.rvf.V2r);
	ArrayOps::Real_space_multiply(W.rvf.V3r, temp1, U.rvf.V3r);
	// U.rvf.Vir = Qtheta
	
	// Div calc
	universal->Forward_transform(U.rvf.V1r, global.temp_array.X);
	universal->Xderiv(global.temp_array.X, temp2);
	
	universal->Forward_transform(U.rvf.V2r, global.temp_array.X);
	universal->Yderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	
	universal->Forward_transform(U.rvf.V3r, global.temp_array.X);
	universal->Zderiv(global.temp_array.X, global.temp_array.X2);
	temp2 += global.temp_array.X2;
	// temp2 = div(Qc)
	
	C.Force += temp2/Peclet_c_astro;
}

//************************ End of compute_force_RB.cc *****************************************

