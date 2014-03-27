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

/*! \file  init_cond_modes.cc
 * 
 * @brief Input several modes as initial condition
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */


#include "FluidIO.h"

//****************************************************************************************


/** @brief Read V(k) for some modes as initial condition.  
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
 * Complex conjugate modes are automatically added.
 */	
void  FluidIO::Init_cond_modes(FluidVF &U)
{
	if (basis_type.find("Ch") != string::npos) {
        Init_cond_modes_Chebyshev(U);
		return;
	}
	
	int kx, ky, kz;	
	complx Vx_complex, Vy_complex, Vz_complex;
	DP Vx_real, Vy_real, Vz_real;
	
	(U.cvf.V1) = 0.0; 
	(U.cvf.V2) = 0.0; 
	(U.cvf.V3) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1); 
		ky = global.io.init_cond_modes.coords(mode,2); 
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0); 
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1); 
			
			if (global.field.incompressible) 
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			
			else
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2); 
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			
			if (my_id == master_id) 
				cout << "k, V : (" << kx << " " << ky  << " " <<kz << ") " << Vx_real	<<	" "	<< Vy_real	<<	" "		<<	Vz_real	<<  endl; 
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0); 
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1); 
			
			if (global.field.incompressible) 
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			
			else
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2); 
			
			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);				
			
			if (my_id == master_id) 
				cout << "k, V : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_complex	<<	" "	<< Vy_complex	<<	" "	<<	Vz_complex	<<  endl; 
		}
	}

}

//*********************************************************************************************


/** @brief Read V(k) and T.csf.F(k) for some modes as initial condition.
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

void  FluidIO::Init_cond_modes(FluidVF&U, FluidSF& T)
{
	
	if (global.program.kind == "INC_SCALAR")
		Init_cond_modes_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Init_cond_modes_RBC(U, T);

}



void  FluidIO::Init_cond_modes_scalar(FluidVF&U, FluidSF& T)
{

	if (basis_type.find("Ch") != string::npos) {
        Init_cond_modes_scalar_Chebyshev(U, T);
		return;
	}
	
	int kx, ky, kz;	
	complx Vx_complex, Vy_complex, Vz_complex, G_complex;
	DP Vx_real, Vy_real, Vz_real, G_real;
	
	(U.cvf.V1) = 0.0; 
	(U.cvf.V2) = 0.0; 
	(U.cvf.V3) = 0.0;
	(T.csf.F) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1); 
		ky = global.io.init_cond_modes.coords(mode,2); 
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0); 
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
				G_real = global.io.init_cond_modes.field_array_real(mode,2);
			}
			
			else {
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2); 
				G_real = global.io.init_cond_modes.field_array_real(mode,3); 
			}
			
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			T.Assign_field(kx, ky, kz, G_real);
			
			if (my_id == master_id)
				cout << "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " << Vx_real  << " "  << Vy_real << " " << Vz_real << " " << G_real <<  endl;
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0); 
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
				G_complex = global.io.init_cond_modes.field_array_complex(mode,2);
			}
			
			else {
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2); 
				G_complex = global.io.init_cond_modes.field_array_complex(mode,3);
			}

			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			T.Assign_field_and_comp_conj(kx, ky, kz, G_complex);
			
			if (my_id == master_id)
				cout << "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " << Vx_complex  << " "  << Vy_complex << " " << Vz_complex << " " << G_complex <<  endl;
		}
	}	
	
}

// RBC
void  FluidIO::Init_cond_modes_RBC(FluidVF&U, FluidSF& T)
{
	if (basis_type.find("Ch") != string::npos) {
        Init_cond_modes_RBC_Chebyshev(U, T);
		return;
	}
	
	
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_modes(U);		
		
		U.Zero_Prandtl_number_compute_temperature(T);
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		
		int kx, ky, kz;
		complx G_complex;
		DP G_real;
		
		(U.cvf.V1) = 0.0; 
		(U.cvf.V2) = 0.0; 
		(U.cvf.V3) = 0.0;
		(T.csf.F) = 0.0;
		
		for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
			kx = global.io.init_cond_modes.coords(mode,1); 
			ky = global.io.init_cond_modes.coords(mode,2); 
			kz = global.io.init_cond_modes.coords(mode,3);
			
			if (basis_type == "SSS") {
				G_real = global.io.init_cond_modes.field_array_real(mode,0);
				T.Assign_field(kx, ky, kz, G_real);
				
				if (my_id == master_id)
					cout << "k, G : (" << kx << " " << ky  << " " <<kz << ") "  << G_real << endl;
			}
				
			else {
				G_complex = global.io.init_cond_modes.field_array_complex(mode,0);
				T.Assign_field_and_comp_conj(kx, ky, kz, G_complex);
				
				if (my_id == master_id)
					cout << "k, G : (" << kx << " " << ky  << " " <<kz << ") "  << G_complex << endl;
			}
		}
		
		U.Infinite_Prandtl_number_compute_velocity(T);
		
	}
	
	else
		Init_cond_modes_scalar(U, T);
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
void  FluidIO::Init_cond_modes(FluidVF& U, FluidVF& W)
{
	
	int kx, ky, kz;	
	complx Vx_complex, Vy_complex, Vz_complex;
	DP Vx_real, Vy_real, Vz_real;
	
	complx Wx_complex, Wy_complex, Wz_complex;
	DP Wx_real, Wy_real, Wz_real;
	
	(U.cvf.V1) = 0.0; 
	(U.cvf.V2) = 0.0; 
	(U.cvf.V3) = 0.0;
	
	(W.cvf.V1) = 0.0; 
	(W.cvf.V2) = 0.0; 
	(W.cvf.V3) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1); 
		ky = global.io.init_cond_modes.coords(mode,2); 
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0); 
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,2); 
				Wy_real = global.io.init_cond_modes.field_array_real(mode,3);
				universal->Last_component(kx, ky, kz, Wx_real, Wy_real, Wz_real);
			}
			
			else {
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,3); 
				Wy_real = global.io.init_cond_modes.field_array_real(mode,4);
				Wz_real = global.io.init_cond_modes.field_array_real(mode,5);
			}
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			W.Assign_field(kx, ky, kz, Wx_real, Wy_real, Wz_real);
			
			if (my_id == master_id)
				cout	<< "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_real	<<	" "		<< Vy_real	<<	" "	<<	Vz_real <<	" || " <<	Wx_real	<<	" "		<< Wy_real	<<	" "	<<	Wz_real	 <<  endl; 
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0); 
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,2); 
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,3);
				universal->Last_component(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
			}
			
			else {
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,3); 
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,4);
				Wz_complex = global.io.init_cond_modes.field_array_complex(mode,5);
			}
			
			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			W.Assign_field_and_comp_conj(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
			
			if (my_id == master_id)
				cout	<< "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_complex	<<	" "		<< Vy_complex	<<	" "	<<	Vz_complex <<	" || " <<	Wx_complex	<<	" "		<< Wy_complex	<<	" "	<<	Wz_complex	 <<  endl; 
		}
	}

}


//*********************************************************************************************

/** @brief Read V(k), W(k), and T.csf.F(k) for some modes as initial condition.
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
void  FluidIO::Init_cond_modes(FluidVF& U, FluidVF& W, FluidSF& T)
{
	int kx, ky, kz;	
	complx Vx_complex, Vy_complex, Vz_complex;
	DP Vx_real, Vy_real, Vz_real;
	
	complx Wx_complex, Wy_complex, Wz_complex;
	DP Wx_real, Wy_real, Wz_real;
	
	complx G_complex;
	DP G_real;
	
	(U.cvf.V1) = 0.0; 
	(U.cvf.V2) = 0.0; 
	(U.cvf.V3) = 0.0;
	
	(W.cvf.V1) = 0.0; 
	(W.cvf.V2) = 0.0; 
	(W.cvf.V3) = 0.0;
	
	(T.csf.F) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1); 
		ky = global.io.init_cond_modes.coords(mode,2); 
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0); 
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,2); 
				Wy_real = global.io.init_cond_modes.field_array_real(mode,3);
				universal->Last_component(kx, ky, kz, Wx_real, Wy_real, Wz_real);
				
				G_real = global.io.init_cond_modes.field_array_real(mode,4);
			}
			
			else {
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,3); 
				Wy_real = global.io.init_cond_modes.field_array_real(mode,4);
				Wz_real = global.io.init_cond_modes.field_array_real(mode,5);
				
				G_real = global.io.init_cond_modes.field_array_real(mode,6);
			}
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			W.Assign_field(kx, ky, kz, Wx_real, Wy_real, Wz_real);
			T.Assign_field(kx, ky, kz, G_real);
			
			if (my_id == master_id)
				cout << "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_real	<<	" "		<< Vy_real	<<	" "	<<	Vz_real <<	" || " <<	Wx_real	<<	" "		<< Wy_real	<<	" "	<<	Wz_real << " " << G_real << endl; 
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0); 
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,2); 
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,3);
				universal->Last_component(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
				
				G_complex = global.io.init_cond_modes.field_array_complex(mode,4);
			}
			
			else {
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,3); 
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,4);
				Wz_complex = global.io.init_cond_modes.field_array_complex(mode,5);
				
				G_complex = global.io.init_cond_modes.field_array_complex(mode,6);
			}
			
			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			W.Assign_field_and_comp_conj(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
			T.Assign_field_and_comp_conj(kx, ky, kz, G_complex);
			
			if (my_id == master_id)
				cout << "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_complex	<<	" "		<< Vy_complex	<<	" "	<<	Vz_complex <<	" || " <<	Wx_complex	<<	" "		<< Wy_complex	<<	" "	<<	Wz_complex << " " << G_complex << endl; 
		}
	}
	
}

//**************************************************************************************

void  FluidIO::Init_cond_modes(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	int kx, ky, kz;
	complx Vx_complex, Vy_complex, Vz_complex;
	DP Vx_real, Vy_real, Vz_real;
	
	complx Wx_complex, Wy_complex, Wz_complex;
	DP Wx_real, Wy_real, Wz_real;
	
	complx G_complex, H_complex;
	DP G_real, H_real;
	
	(U.cvf.V1) = 0.0;
	(U.cvf.V2) = 0.0;
	(U.cvf.V3) = 0.0;
	
	(W.cvf.V1) = 0.0;
	(W.cvf.V2) = 0.0;
	(W.cvf.V3) = 0.0;
	
	(T.csf.F) = 0.0;
	(C.csf.F) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1);
		ky = global.io.init_cond_modes.coords(mode,2);
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0);
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1);
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,2);
				Wy_real = global.io.init_cond_modes.field_array_real(mode,3);
				universal->Last_component(kx, ky, kz, Wx_real, Wy_real, Wz_real);
				
				G_real = global.io.init_cond_modes.field_array_real(mode,4);
				H_real = global.io.init_cond_modes.field_array_real(mode,5);
			}
			
			else {
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2);
				
				Wx_real = global.io.init_cond_modes.field_array_real(mode,3);
				Wy_real = global.io.init_cond_modes.field_array_real(mode,4);
				Wz_real = global.io.init_cond_modes.field_array_real(mode,5);
				
				G_real = global.io.init_cond_modes.field_array_real(mode,6);
				H_real = global.io.init_cond_modes.field_array_real(mode,7);
			}
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			W.Assign_field(kx, ky, kz, Wx_real, Wy_real, Wz_real);
			T.Assign_field(kx, ky, kz, G_real);
			C.Assign_field(kx, ky, kz, H_real);
			
			if (my_id == master_id)
				cout << "k, V, W, G : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_real	<<	" "		<< Vy_real	<<	" "	<<	Vz_real <<	" || " <<	Wx_real	<<	" "		<< Wy_real	<<	" "	<<	Wz_real << " " << G_real << " " << H_real << endl;
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0);
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1);
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,3);
				universal->Last_component(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
				
				G_complex = global.io.init_cond_modes.field_array_complex(mode,4);
				H_complex = global.io.init_cond_modes.field_array_complex(mode,5);
			}
			
			else {
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				
				Wx_complex = global.io.init_cond_modes.field_array_complex(mode,3);
				Wy_complex = global.io.init_cond_modes.field_array_complex(mode,4);
				Wz_complex = global.io.init_cond_modes.field_array_complex(mode,5);
				
				G_complex = global.io.init_cond_modes.field_array_complex(mode,6);
				H_complex = global.io.init_cond_modes.field_array_complex(mode,7);
			}
			
			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			W.Assign_field_and_comp_conj(kx, ky, kz, Wx_complex, Wy_complex, Wz_complex);
			T.Assign_field_and_comp_conj(kx, ky, kz, G_complex);
			C.Assign_field_and_comp_conj(kx, ky, kz, H_complex);
			
			if (my_id == master_id)
				cout << "k, V, W, C : (" << kx << " " << ky  << " " <<kz << ") " <<	Vx_complex	<<	" "		<< Vy_complex	<<	" "	<<	Vz_complex <<	" || " <<	Wx_complex	<<	" "		<< Wy_complex	<<	" "	<<	Wz_complex << " " << G_complex << " " << H_complex << endl;
		}
	}
	
}


//******************************************************************************

void  FluidIO::Init_cond_modes(FluidVF&U, FluidSF& T1, FluidSF& T2)
{
	
	int kx, ky, kz;	
	complx Vx_complex, Vy_complex, Vz_complex, G1_complex, G2_complex;
	DP Vx_real, Vy_real, Vz_real, G1_real, G2_real;
	
	(U.cvf.V1) = 0.0; 
	(U.cvf.V2) = 0.0; 
	(U.cvf.V3) = 0.0;
	(T1.csf.F) = 0.0;
	(T2.csf.F) = 0.0;
	
	for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1); 
		ky = global.io.init_cond_modes.coords(mode,2); 
		kz = global.io.init_cond_modes.coords(mode,3);
		
		if (basis_type == "SSS") {
			Vx_real = global.io.init_cond_modes.field_array_real(mode,0); 
			Vy_real = global.io.init_cond_modes.field_array_real(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_real, Vy_real, Vz_real);
				G1_real = global.io.init_cond_modes.field_array_real(mode,2);
				G2_real = global.io.init_cond_modes.field_array_real(mode,3);
			}
			
			else {
				Vz_real = global.io.init_cond_modes.field_array_real(mode,2); 
				G1_real = global.io.init_cond_modes.field_array_real(mode,3); 
				G2_real = global.io.init_cond_modes.field_array_real(mode,4); 
			}
			
			U.Assign_field(kx, ky, kz, Vx_real, Vy_real, Vz_real);
			T1.Assign_field(kx, ky, kz, G1_real);
			T2.Assign_field(kx, ky, kz, G2_real);
			
			if (my_id == master_id)
				cout	<< "k, V, G1, G2 : (" << kx << " " << ky  << " " <<kz << ") " << Vx_real  << " "  << Vy_real << " " << Vz_real << " " << G1_real << " " << G2_real << endl;
		}
		
		else {
			Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0); 
			Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1); 
			
			if (global.field.incompressible) {
				universal->Last_component(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
				G1_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				G2_complex = global.io.init_cond_modes.field_array_complex(mode,3);
			}
			
			else {
				Vz_complex = global.io.init_cond_modes.field_array_complex(mode,2); 
				G1_complex = global.io.init_cond_modes.field_array_complex(mode,3);
				G2_complex = global.io.init_cond_modes.field_array_complex(mode,4);
			}
			
			U.Assign_field_and_comp_conj(kx, ky, kz, Vx_complex, Vy_complex, Vz_complex);
			T1.Assign_field_and_comp_conj(kx, ky, kz, G1_complex);
			T2.Assign_field_and_comp_conj(kx, ky, kz, G2_complex);
			
			if (my_id == master_id)
				cout	<< "k, V, G1, G2 : (" << kx << " " << ky  << " " <<kz << ") " << Vx_complex  << " "  << Vy_complex << " " << Vz_complex << " " << G1_complex <<  " " << G2_complex <<  endl;
		}
	}		
}


//******************************************************************************
//******************************************************************************

void  FluidIO::Init_cond_modes_Chebyshev(FluidVF&U)
{
	int kx, ky, kz;
    DP Kx,Ky,Kz;
    complx Vx_complex, Vy_complex, Vz_complex;
    DP Vx_real, Vy_real, Vz_real, dvxdx;
    int mode;
	
	(U.cvf.V1) = 0.0;
	(U.cvf.V2) = 0.0;
	(U.cvf.V3) = 0.0;
	
    if (basis_type == "ChSS") {
        for (mode = 0; mode < global.io.init_cond_modes.number; mode++) {
            kx = global.io.init_cond_modes.coords(mode,1);
            ky = global.io.init_cond_modes.coords(mode,2);
            kz = global.io.init_cond_modes.coords(mode,3);
            
            Vx_real = global.io.init_cond_modes.field_array_real(mode,0);
            
            if (ky>=0)
                universal->Assign_spectral_field(kx,ky,kz, U.cvf.V1, Vx_real);
            else
                cerr << "Init_cond_modes_Chebyshev(U) demands that ky>=0) " << endl;
        }
        
        universal->Xderiv(U.cvf.V1, global.temp_array.X);
        
        for (int ly=0; ly<(U.cvf.V1).extent(0); ly++)
            for (int lz=0; lz<2*(U.cvf.V1).extent(1); lz++)
                for (int lx=0; lx<(U.cvf.V1).extent(2); lx++) {
                    kx = universal->Get_kx(lx);
                    ky = universal->Get_ky(ly);
                    kz = universal->Get_kz(lz);
                    
                    Ky = ky*kfactor[2];
                    Kz = kz*kfactor[3];
                    
                    mode = global.io.init_cond_modes.In_table(kx,ky,kz);
                    // if (kx,ky,kz) is input mode, then mode yields the table index.
                    
                    if (mode > -1)
                        Vy_real = global.io.init_cond_modes.field_array_real(mode,1);
                    else // no input mode
                        Vy_real = 0.0;
                    
                    if (lz%2 == 0)
                        dvxdx = real(global.temp_array.X(ly,lz,lx));
                    else
                        dvxdx = imag(global.temp_array.X(ly,lz,lx));
                    
                    universal->Last_component(kx, ky, kz, dvxdx, Vy_real, Vz_real);
                    
                    universal->Assign_local_spectral_field(lx,ly,lz,U.cvf.V2, Vy_real);
                    universal->Assign_local_spectral_field(lx,ly,lz,U.cvf.V3, Vz_real);
                }
    }
    
    else if ((basis_type == "ChFF") || (basis_type == "ChSF"))  {
        for (mode = 0; mode < global.io.init_cond_modes.number; mode++) {
            kx = global.io.init_cond_modes.coords(mode,1);
            ky = global.io.init_cond_modes.coords(mode,2);
            kz = global.io.init_cond_modes.coords(mode,3);
            
            Vx_complex = global.io.init_cond_modes.field_array_complex(mode,0);
            
            if (ky>=0)
                universal->Assign_spectral_field(kx,ky,kz, U.cvf.V1, Vx_complex);
            else
                cerr << "Init_cond_modes_Chebyshev(U) demands that ky>=0) " << endl;
        }
        
        universal->Xderiv(U.cvf.V1, global.temp_array.X);
        
        for (int ly=0; ly<(U.cvf.V1).extent(0); ly++)
            for (int lz=0; lz<(U.cvf.V1).extent(1); lz++)
                for (int lx=0; lx<(U.cvf.V1).extent(2); lx++) {
                    kx = universal->Get_kx(lx);
                    ky = universal->Get_ky(ly);
                    kz = universal->Get_kz(lz);
                    
                    Ky = ky*kfactor[2];
                    Kz = kz*kfactor[3];
                    
                    mode = global.io.init_cond_modes.In_table(kx,ky,kz);
                    // if (kx,ky,kz) is input mode, then mode yields the table index.
                    
                    if (mode > -1)
                        Vy_complex = global.io.init_cond_modes.field_array_complex(mode,1);
                    else // no input mode
                        Vy_complex = complx(0.0,0.0);
                        
                    universal->Last_component(kx, ky, kz, global.temp_array.X(ly,lz,lx), Vy_complex, Vz_complex);
                    
                    universal->Assign_local_spectral_field(lx,ly,lz,U.cvf.V2, Vy_complex);
                    universal->Assign_local_spectral_field(lx,ly,lz,U.cvf.V3, Vz_complex);
                }
    }
	
	// V(kx,-ky,0) = V^*(kx,ky,0)
	U.Satisfy_strong_reality_condition_field();
}

//******************************************************************************
void  FluidIO::Init_cond_modes_scalar_Chebyshev(FluidVF&U, FluidSF& T)
{
	Init_cond_modes_Chebyshev(U);
    
   	(T.csf.F) = 0.0;
    
    int kx, ky, kz;
    complx G_complex;
	DP G_real;
	
    for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
		kx = global.io.init_cond_modes.coords(mode,1);
		ky = global.io.init_cond_modes.coords(mode,2);
		kz = global.io.init_cond_modes.coords(mode,3);
        
		if ((basis_type == "ChFF") || (basis_type == "ChSF")) {
			G_complex = global.io.init_cond_modes.field_array_complex(mode,2);
			
			universal->Assign_spectral_field(kx,ky,kz, T.csf.F, G_complex);
			
			if (master)
				cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " <<  G_complex <<  endl;
		}
		else if (basis_type == "ChSS") {
			G_real = global.io.init_cond_modes.field_array_real(mode,2);
			
			universal->Assign_spectral_field(kx,ky,kz, T.csf.F, G_real);
			
			if (master)
				cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " <<  G_real <<  endl;
		}
    }
}

//******************************************************************************
void  FluidIO::Init_cond_modes_RBC_Chebyshev(FluidVF&U, FluidSF& T)
{
	
	if (global.PHYSICS.Pr_option == "PRZERO") {
		Init_cond_modes_Chebyshev(U);
		
		T.csf.F = 0.0;
	}
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {
		(U.cvf.V1) = 0.0;
		(U.cvf.V2) = 0.0;
		(U.cvf.V3) = 0.0;
		
		int kx, ky, kz;
		complx G_complex;
		DP G_real;
		
		for (int mode = 0; mode < global.io.init_cond_modes.number; mode++) {
			kx = global.io.init_cond_modes.coords(mode,1);
			ky = global.io.init_cond_modes.coords(mode,2);
			kz = global.io.init_cond_modes.coords(mode,3);
			
			if ((basis_type == "ChFF") || (basis_type == "ChSF")) {
				G_complex = global.io.init_cond_modes.field_array_complex(mode,2);
				
				universal->Assign_spectral_field(kx,ky,kz, T.csf.F, G_complex);
				
				if (master)
					cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " <<  G_complex <<  endl;
			}
			else if (basis_type == "ChSS") {
				G_real = global.io.init_cond_modes.field_array_real(mode,2);
				
				universal->Assign_spectral_field(kx,ky,kz, T.csf.F, G_real);
				
				if (master)
					cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " <<  G_real <<  endl;
			}
		}
	}
	
	else
		Init_cond_modes_scalar_Chebyshev(U, T);
	
}

//******************************  End of Init_cond_modes.cc  **********************************



