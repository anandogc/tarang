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
/*! \file  compute_force_modes.cc
 * 
 * @brief Force given set of modes only
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  INITIALIZW Fx, Fy, Fz, F etc.
 */

#include "FORCE.h"

//// NOTE:: If that particular force does not exist.. make the components zero.
//// Read at t=0.


//****************************************************************************************

/** @brief Read F(k) for some modes as initial condition.  
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ F_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Fy computed using \f$ F_y = -D_x F_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ F_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * Complex conjugate modes are automatically added.
 */	
void FORCE::Compute_force_given_modes(FluidVF& U)
{
	if (!global.force.modes.read_done) {

		int kx, ky, kz;	
		Complex Fx_complex, Fy_complex, Fz_complex;
		Real Fx_real, Fy_real, Fz_real;
		
		(U.Force1) = 0.0; 
		(U.Force2) = 0.0; 
		(U.Force3) = 0.0;
		
		for (int mode = 0; mode < global.force.modes.number; mode++) {
			kx = global.force.modes.coords(mode,1); 
			ky = global.force.modes.coords(mode,2); 
			kz = global.force.modes.coords(mode,3);
			
			if (basis_type == "SSS") {
				Fx_real = global.force.modes.field_array_real(mode,0); 
				Fy_real = global.force.modes.field_array_real(mode,1); 
				
				if (global.field.incompressible) 
					universal->Last_component(kx, ky, kz, Fx_real, Fy_real, Fz_real);
				
				else
					Fz_real = global.force.modes.field_array_real(mode,2); 
				
				if (U.force_switch)
					Assign_force(U, kx, ky, kz, Fx_real, Fy_real, Fz_real);
				
				if (my_id == master_id) 
					cout	<< "k, V : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_real	<<	" "		<< Fy_real	<<	" "	<<	Fz_real	<<  endl; 
			}
			
			else {
				Fx_complex = global.force.modes.field_array_complex(mode,0); 
				Fy_complex = global.force.modes.field_array_complex(mode,1); 
				
				if (global.field.incompressible) 
					universal->Last_component(kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
				
				else
					Fz_complex = global.force.modes.field_array_complex(mode,2); 
				
				if (U.force_switch)
					Assign_force_and_comp_conj(U, kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);				
				
				if (my_id == master_id) 
					cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_complex	<<	" "		<< Fy_complex	<<	" "	<<	Fz_complex	<<  endl; 
			}
		}
	}
	
	global.force.modes.read_done = true;
}

//*********************************************************************************************


/** @brief Read V(k) and T.F(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Scalar field is added simply.
 */	 
void FORCE::Compute_force_given_modes(FluidVF& U, FluidSF& T)
{	
	if (!global.force.modes.read_done) {

		int kx, ky, kz;	
		Complex Fx_complex, Fy_complex, Fz_complex, G_complex;
		Real Fx_real, Fy_real, Fz_real, G_real;
		
		(U.Force1) = 0.0; 
		(U.Force2) = 0.0; 
		(U.Force3) = 0.0;
		(T.Force) = 0.0;
		
		for (int mode = 0; mode < global.force.modes.number; mode++) {
			kx = global.force.modes.coords(mode,1); 
			ky = global.force.modes.coords(mode,2); 
			kz = global.force.modes.coords(mode,3);
			
			if (basis_type == "SSS") {
				Fx_real = global.force.modes.field_array_real(mode,0); 
				Fy_real = global.force.modes.field_array_real(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_real, Fy_real, Fz_real);
					G_real = global.force.modes.field_array_real(mode,2);
				}
				
				else {
					Fz_real = global.force.modes.field_array_real(mode,2); 
					G_real = global.force.modes.field_array_real(mode,3);
				}
				
				if (U.force_switch)
					Assign_force(U, kx, ky, kz, Fx_real, Fy_real, Fz_real);
				
				if (T.force_switch)
					Assign_force(T, kx, ky, kz, G_real);
				
				if (my_id == master_id) 
					cout	<< "k, Force, T.Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_real	<<	" "		<< Fy_real	<<	" "	<<	Fz_real	<<  "   " << G_real << endl; 
			}
			
			else {
				Fx_complex = global.force.modes.field_array_complex(mode,0); 
				Fy_complex = global.force.modes.field_array_complex(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
					G_complex = global.force.modes.field_array_complex(mode,2); 
				}
				
				else {
					Fz_complex = global.force.modes.field_array_complex(mode,2);
					G_complex = global.force.modes.field_array_complex(mode,3);
				}
				
				if (U.force_switch)
					Assign_force_and_comp_conj(U, kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
				
				if (T.force_switch)
					Assign_force_and_comp_conj(T, kx, ky, kz, G_complex);
				
				if (my_id == master_id) 
					cout	<< "k, Force, T.Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_complex	<<	" "		<< Fy_complex	<<	" "		<<	Fz_complex	<< " " << G_complex << endl; 
			}
		}
	}
	
	global.force.modes.read_done = true;
}


//*********************************************************************************************

/** @brief Read V(k) and W(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Vector field W is added in a similar manner.
 */		
void FORCE::Compute_force_given_modes(FluidVF& U, FluidVF& W)
{	
	if (!global.force.modes.read_done) {
		
		int kx, ky, kz;	
		Complex Fx_complex, Fy_complex, Fz_complex;
		Complex FxW_complex, FyW_complex, FzW_complex;
		Real Fx_real, Fy_real, Fz_real;
		Real FxW_real, FyW_real, FzW_real;
		
		(U.Force1) = 0.0; 
		(U.Force2) = 0.0; 
		(U.Force3) = 0.0;
		(W.Force1) = 0.0; 
		(W.Force2) = 0.0; 
		(W.Force3) = 0.0;
		
		for (int mode = 0; mode < global.force.modes.number; mode++) {
			kx = global.force.modes.coords(mode,1); 
			ky = global.force.modes.coords(mode,2); 
			kz = global.force.modes.coords(mode,3);
			
			if (basis_type == "SSS") {
				Fx_real = global.force.modes.field_array_real(mode,0); 
				Fy_real = global.force.modes.field_array_real(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_real, Fy_real, Fz_real);
					
					FxW_real = global.force.modes.field_array_real(mode,2); 
					FyW_real = global.force.modes.field_array_real(mode,3);
					
					universal->Last_component(kx, ky, kz, FxW_real, FyW_real, FzW_real);
				}
				
				else {
					Fz_real = global.force.modes.field_array_real(mode,2); 
					
					FxW_real = global.force.modes.field_array_real(mode,3); 
					FyW_real = global.force.modes.field_array_real(mode,4); 
					FzW_real = global.force.modes.field_array_real(mode,5);
				}
				
				if (U.force_switch)
					Assign_force(U, kx, ky, kz, Fx_real, Fy_real, Fz_real);
				
				if (W.force_switch)
					Assign_force(W, kx, ky, kz, FxW_real, FyW_real, FzW_real);
				
				if (my_id == master_id) 
					cout	<< "k, Force, W.Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_real	<<	" " << Fy_real	<<	" "	<<	Fz_real	<<  "   " << FxW_real << " " << FyW_real <<	" "	<<	FzW_real <<	endl; 
			}
			
			else {
				Fx_complex = global.force.modes.field_array_complex(mode,0); 
				Fy_complex = global.force.modes.field_array_complex(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
					
					FxW_complex = global.force.modes.field_array_complex(mode,2); 
					FyW_complex = global.force.modes.field_array_complex(mode,3);
					universal->Last_component(kx, ky, kz, FxW_complex, FyW_complex, FzW_complex);
				}
				
				else {
					Fz_complex = global.force.modes.field_array_complex(mode,2);
					
					FxW_complex = global.force.modes.field_array_complex(mode,3); 
					FyW_complex = global.force.modes.field_array_complex(mode,4);
					FzW_complex = global.force.modes.field_array_complex(mode,5);
				}
				
				if (U.force_switch)
					Assign_force_and_comp_conj(U, kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
				
				if (W.force_switch)
					Assign_force_and_comp_conj(W, kx, ky, kz, FxW_complex, FyW_complex, FzW_complex);
				
				if (my_id == master_id) 
					cout	<< "k, Force, W.Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_complex	<<	" "		<< Fy_complex	<<	" "		<<	Fz_complex	<< " " << FxW_complex	<<	" "		<< FyW_complex	<<	" "		<<	FzW_complex	 << endl; 
			}
		}
	}
	
	global.force.modes.read_done = true;
}


//*********************************************************************************************

/** @brief Read V(k), W(k), and T.F(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Vector field W and scalar field T are added in a similar manner.
 */
void FORCE::Compute_force_given_modes(FluidVF& U, FluidVF& W, FluidSF& T)
{	
	if  (!global.force.modes.read_done) {
		
		int kx, ky, kz;	
		Complex Fx_complex, Fy_complex, Fz_complex;
		Complex FxW_complex, FyW_complex, FzW_complex, G_complex;
		Real Fx_real, Fy_real, Fz_real;
		Real FxW_real, FyW_real, FzW_real, G_real;
		
		(U.Force1) = 0.0; 
		(U.Force2) = 0.0; 
		(U.Force3) = 0.0;
		(W.Force1) = 0.0; 
		(W.Force2) = 0.0; 
		(W.Force3) = 0.0;
		(T.Force) = 0.0;
		
		for (int mode = 0; mode < global.force.modes.number; mode++) {
			kx = global.force.modes.coords(mode,1); 
			ky = global.force.modes.coords(mode,2); 
			kz = global.force.modes.coords(mode,3);
			
			if (basis_type == "SSS") {
				Fx_real = global.force.modes.field_array_real(mode,0); 
				Fy_real = global.force.modes.field_array_real(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_real, Fy_real, Fz_real);
					
					FxW_real = global.force.modes.field_array_real(mode,2); 
					FyW_real = global.force.modes.field_array_real(mode,3);
					universal->Last_component(kx, ky, kz, FxW_real, FyW_real, FzW_real);
					
					G_real = global.force.modes.field_array_real(mode,4);
				}
				
				else {
					Fz_real = global.force.modes.field_array_real(mode,2); 
					
					FxW_real = global.force.modes.field_array_real(mode,3); 
					FyW_real = global.force.modes.field_array_real(mode,4); 
					FzW_real = global.force.modes.field_array_real(mode,5);
					
					G_real = global.force.modes.field_array_real(mode,6);
				}
				
				if (U.force_switch)
					Assign_force(U, kx, ky, kz, Fx_real, Fy_real, Fz_real);
				
				if (W.force_switch)
					Assign_force(W, kx, ky, kz, FxW_real, FyW_real, FzW_real);
				
				if (T.force_switch)
					Assign_force(T, kx, ky, kz, G_real);
				
				if (my_id == master_id) 
					cout	<< "k, Force, ForceW, T.Force : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_real	<<	" " << Fy_real	<<	" "	<<	Fz_real	<<  "   " << FxW_real << " " << FyW_real <<	" "	<<	FzW_real <<  " " << G_real <<	endl; 
			}
			
			else {
				Fx_complex = global.force.modes.field_array_complex(mode,0); 
				Fy_complex = global.force.modes.field_array_complex(mode,1); 
				
				if (global.field.incompressible) {
					universal->Last_component(kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
					
					FxW_complex = global.force.modes.field_array_complex(mode,2); 
					FyW_complex = global.force.modes.field_array_complex(mode,3);
					universal->Last_component(kx, ky, kz, FxW_complex, FyW_complex, FzW_complex);
					
					G_complex = global.force.modes.field_array_complex(mode,4);
				}
				
				else {
					Fz_complex = global.force.modes.field_array_complex(mode,2);
					
					FxW_complex = global.force.modes.field_array_complex(mode,3); 
					FyW_complex = global.force.modes.field_array_complex(mode,4);
					FzW_complex = global.force.modes.field_array_complex(mode,5);
					
					G_complex = global.force.modes.field_array_complex(mode,6);
				}
				
				if (U.force_switch)
					Assign_force_and_comp_conj(U, kx, ky, kz, Fx_complex, Fy_complex, Fz_complex);
				
				if (W.force_switch)
					Assign_force_and_comp_conj(W, kx, ky, kz, FxW_complex, FyW_complex, FzW_complex);
				
				if (T.force_switch)
					Assign_force_and_comp_conj(T, kx, ky, kz, G_complex);
				
				if (my_id == master_id) 
					cout	<< "k, ForceV, ForceW, ForceT : (" << kx << " " << ky  << " " <<kz << ") " <<	Fx_complex	<<	" "		<< Fy_complex	<<	" "	<<	Fz_complex	<< " " << FxW_complex	<<	" "	<< FyW_complex	<<	" "		<<	FzW_complex	 << " " << G_complex << endl; 
			}
		}
	}
	
	global.force.modes.read_done = true;
}


//*******************************  End of compute_force_modes.cc  *****************************


