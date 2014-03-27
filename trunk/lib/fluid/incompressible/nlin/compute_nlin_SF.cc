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

/*! \file  compute_nlin_SF.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j V_j V_i) \f$ and 
 *		\f$ N_i^\psi = \mathcal{F} (D_j V_j \psi) \f$.
 *
 *	Steps <BR>
 * # Inverse transform of V & T: \f$ \vec{Vr} = \mathcal{F}(\vec{V}) \f$  & 
 *			\f$ Fr = \mathcal{F}(F) \f$  <BR>
 * # T.nlin \f$ N^\psi = \mathcal{F} (D_j V_j \psi) \f$. <BR>
 * # nlin \f$ N_i  = (Vr_i)^2 \f$  <BR>
 * # nlin \f$ N_i  = \mathcal{F}(D_i (Vr_i)^2) \f$ <BR>
 * # compute off-diagonal terms and put them in \f$ \vec{Vr} \f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j V_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *  
 * # RB Convection: If Pr==0: nlin(), else nlin(T).
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   Transpose order needs to tested
 */

  
#include "Nlin.h"
#include "IncIO.h"
#include "extern_var_incompress.h"

//**************************************************************************************************

void Nlin_incompress::Compute_nlin(FluidVF& U, FluidSF& T) 
{
	if (global.program.kind == "INC_SCALAR")
		Compute_nlin_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind ==  "STRATIFIED")
		Compute_nlin_RBC(U, T);
}

//*********************************************************************************************									

void Nlin_incompress::Compute_nlin_scalar(FluidVF& U, FluidSF& T) 
{
	
	// Dealiasing reqd before the real-space product.
    if (global.program.alias_option == "DEALIAS") {
        U.cvf.Dealias();
        T.csf.Dealias();
    }
    
	
	U.Inverse_transform();	// Vir = Inv_transform(Vi)
	T.Inverse_transform();	// Fr = Inv_transform(F)

	global.io.real_space_field_available = true;
    
    if (!global.time.dt_computation_done) {
		global.time.dt = U.Get_dt(T);
		global.time.now = global.time.now + global.time.dt;
		global.time.dt_computation_done = true;
	}
    
    // Output real field
    if (!global.io.output_real_field_done) {
		fluidIO_incompress.Output_real_field(U, T);
		fluidIO_incompress.Output_field_r(U, T);
		fluidIO_incompress.Output_global_real_space(U, T);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U, T);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
	
	global.io.real_space_field_available = false;
	
	Compute_nlin_VxT(U, T);  // Computes only T.nlin = FT(Dj(ujT))
    
	Compute_nlin_diag(U);
	
	Compute_nlin_offdiag(U);
	
	if (!global.io.output_field_k_done) {
		fluidIO_incompress.Output_field_k(U, T);
        fluidIO_incompress.Output_Tk_shell_spectrum(U, T);
        fluidIO_incompress.Output_Tk_ring_spectrum(U, T);
        fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, T);
		global.io.output_field_k_done = true;
	} 
    
	

/*	global.io.real_space_field_available = false;
	
    // Output real field
    if (!global.io.output_real_field_done) {
	//	fluidIO_incompress.Output_real_field(U, T);
	//	fluidIO_incompress.Output_field_r(U, T);
		fluidIO_incompress.Output_global_real_space(U, T);  // Only for Cheby basis
		fluidIO_incompress.Output_cout_real_space(U, T);	 // Only for Cheby basis
		global.io.output_real_field_done = true;
	}
	
	global.time.now = global.time.now + global.time.dt;
	
    U.nlin1 = 0.0;
    U.nlin2 = 0.0;
    U.nlin3 = 0.0;
    T.nlin = 0.0; 
 */
}


//*********************************************************************************************
// For RB convection
//

void Nlin_incompress::Compute_nlin_RBC(FluidVF& U, FluidSF& T) 
{
	// For chebyshev
	if (global.program.basis_type.find("Ch") != string::npos)
	{
		if (global.PHYSICS.Pr_option == "PRZERO") {
			Compute_nlin(U);
			T.nlin = 0;
		}
		// if real-field is to be written, only V(r)
		
		else if (global.PHYSICS.Pr_option == "PRINFTY")	{
			// nlin1 = 0.0; nlin2 = 0.0; nlin3 = 0.0;
			
			// Dealiasing reqd before the real-space product.
			if (global.program.alias_option == "DEALIAS")
				T.csf.Dealias();
			
			T.Inverse_transform();	// Fr = Inv_transform(F)
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			
			if (!global.time.dt_computation_done) {
				global.time.dt = U.Get_dt();
				global.time.now = global.time.now + global.time.dt;
				global.time.dt_computation_done = true;
			}
			
			if (!global.io.output_real_field_done) {
				fluidIO_incompress.Output_real_field(U, T);
				fluidIO_incompress.Output_field_r(U, T);
				fluidIO_incompress.Output_global_real_space(U);  // Only for Cheby basis
				fluidIO_incompress.Output_cout_real_space(U);	 // Only for Cheby basis
				global.io.output_real_field_done = true;
			}
			
			Compute_nlin_VxT(U, T);		// For T.nlin only
			
			U.nlin1 = 0.0;
			U.nlin2 = 0.0;
			U.nlin3 = 0.0;
		}
		
		else
			Compute_nlin_scalar(U, T);
        
        if (!global.io.output_field_k_done) {
            fluidIO_incompress.Output_field_k(U, T);
            fluidIO_incompress.Output_Tk_shell_spectrum(U, T);
            fluidIO_incompress.Output_Tk_ring_spectrum(U, T);
            fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, T);
            global.io.output_field_k_done = true;
        }
		
#ifdef GROSSMANN_LOHSE	
		if (!global.io.output_nlin_magnitude_done && global.time.now >= global.io.time.global_save_next) {

            DP visc1, visc2, visc3, Tvisc;
            universal->Laplacian(1.0, U.cvf.V1, global.temp_array.X);
            visc1 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            
            if (!global.program.two_dimension) {
                universal->Laplacian(1.0, U.cvf.V2, global.temp_array.X);
                visc2 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            }
            
            universal->Laplacian(1.0, U.cvf.V3, global.temp_array.X);
            visc3 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            
            universal->Laplacian(1.0, T.csf.F, global.temp_array.X);
            Tvisc = (-T.diffusion_coefficient)* sum(abs(global.temp_array.X));
            
            
            
            DP U_nlin1_local = sum(sqr(abs(U.nlin1)));
            DP U_nlin2_local = sum(sqr(abs(U.nlin2)));
            DP U_nlin3_local = sum(sqr(abs(U.nlin3)));
       
            DP U_nlin1_total;
            DP U_nlin2_total;
            DP U_nlin3_total;
            
            DP T_nlin_local = sum(sqr(abs(T.nlin)));
            DP T_nlin_total;
            
            DP visc1_total;
            DP visc2_total;
            DP visc3_total;
            
            DP Tvisc_total;
              
			MPI_Reduce(&U_nlin1_local, &U_nlin1_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&U_nlin2_local, &U_nlin2_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&U_nlin3_local, &U_nlin3_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&T_nlin_local, &T_nlin_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&visc1, &visc1_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&visc2, &visc2_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&visc3, &visc3_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&Tvisc, &Tvisc_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			U_nlin1_total = sqrt(U_nlin1_total);
			U_nlin2_total = sqrt(U_nlin2_total);
			U_nlin3_total = sqrt(U_nlin3_total);
			
            T_nlin_total = sqrt(T_nlin_local);
			
			fluidIO_incompress.misc_file << "rms(nlin): " << global.time.now << " " << U_nlin1_total << " " << U_nlin2_total << " " << U_nlin3_total << " " << T_nlin_total << " " << visc1_total << " " << visc2_total << " " << visc3_total << " " << Tvisc_total << endl;
			
			global.io.output_nlin_magnitude_done = true;
		}
#endif
	}
	
	//**********
	// Other than ChFF basis
	else {
		if (global.PHYSICS.Pr_option == "PRZERO")
			Compute_nlin(U);   
			// if real-field is to be written, only V(r)
		
		else if (global.PHYSICS.Pr_option == "PRINFTY")	{
			// nlin1 = 0.0; nlin2 = 0.0; nlin3 = 0.0;
			
			// Dealiasing reqd before the real-space product.
			if (global.program.alias_option == "DEALIAS") 
				T.csf.Dealias();
		
			T.Inverse_transform();	// Fr = Inv_transform(F)
			U.Infinite_Prandtl_number_compute_velocity(T); // find U
			U.Inverse_transform();	// Vir = Inv_transform(Vi)
			
			if (!global.time.dt_computation_done) {
				global.time.dt = U.Get_dt();
				global.time.now = global.time.now + global.time.dt;
				global.time.dt_computation_done = true;
			}
			
			if (!global.io.output_real_field_done) {
				fluidIO_incompress.Output_real_field(U, T);
				fluidIO_incompress.Output_field_r(U, T);
				global.io.output_real_field_done = true;
			}
			
			Compute_nlin_VxT(U, T);		// For T.nlin only
		}
		
		else
			Compute_nlin_scalar(U, T);
        
        if (!global.io.output_field_k_done) {
            fluidIO_incompress.Output_field_k(U, T);
            fluidIO_incompress.Output_Tk_shell_spectrum(U, T);
            fluidIO_incompress.Output_Tk_ring_spectrum(U, T);
            fluidIO_incompress.Output_Tk_cylindrical_ring_spectrum(U, T);
            global.io.output_field_k_done = true;
        } 

#ifdef GROSSMANN_LOHSE         
		if (!global.io.output_nlin_magnitude_done && global.time.now >= global.io.time.global_save_next) {
			         
			DP visc1, visc2, visc3, Tvisc;
            universal->Laplacian(1.0, U.cvf.V1, global.temp_array.X);
            visc1 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            
            if (!global.program.two_dimension) {
                universal->Laplacian(1.0, U.cvf.V2, global.temp_array.X);
                visc2 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            }
            
            universal->Laplacian(1.0, U.cvf.V3, global.temp_array.X);
            visc3 = (-U.dissipation_coefficient)* sum(abs(global.temp_array.X));
            
            universal->Laplacian(1.0, T.csf.F, global.temp_array.X);
            
            Tvisc = (-T.diffusion_coefficient)* sum(abs(global.temp_array.X));
            
            //cout << my_id << " " << global.time.now << " " << sum(abs(T.csf.F)) << endl;
            
			DP U_nlin1_local = sum(sqr(abs(U.nlin1)));
            DP U_nlin2_local = sum(sqr(abs(U.nlin2)));
            DP U_nlin3_local = sum(sqr(abs(U.nlin3)));
       
            DP U_nlin1_total;
            DP U_nlin2_total;
            DP U_nlin3_total;
            
            DP T_nlin_local = sum(sqr(abs(T.nlin)));
            DP T_nlin_total;
            
            DP visc1_total;
            DP visc2_total;
            DP visc3_total;
            
            DP Tvisc_total;
            
			MPI_Reduce(&U_nlin1_local, &U_nlin1_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&U_nlin2_local, &U_nlin2_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&U_nlin3_local, &U_nlin3_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&T_nlin_local, &T_nlin_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			
			
			MPI_Reduce(&visc1, &visc1_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&visc2, &visc2_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&visc3, &visc3_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&Tvisc, &Tvisc_total, 1, MPI_DP, MPI_SUM, 0, MPI_COMM_WORLD);
			
			U_nlin1_total = sqrt(U_nlin1_total);
			U_nlin2_total = sqrt(U_nlin2_total);
			U_nlin3_total = sqrt(U_nlin3_total);
			
            T_nlin_total = sqrt(T_nlin_local);

			
			// For RBC GL scaling analysis only.  Comment else
			fluidIO_incompress.misc_file << "rms(nlin): " << global.time.now << " " << U_nlin1_total << " " << U_nlin2_total << " " << U_nlin3_total << " " << T_nlin_total << " " << visc1_total << " " << visc2_total << " " << visc3_total << " " << Tvisc_total << endl;
			
			global.io.output_nlin_magnitude_done = true;
		}
#endif
    }
	
/*	// Linear theory
	U.nlin1 = 0.0;
	U.nlin2 = 0.0;
	U.nlin3 = 0.0;
	T.nlin = 0.0; */
	
}		

//***************************** fn compute_nlin_SF ends ***************************************


