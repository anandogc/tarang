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

/*! \file  Cvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "FluidSF.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
   
**********************************************************************************************/

FluidSF::FluidSF
(
	Real diffusion_coefficient, 
	Real hyper_diffusion_coefficient, 
	int hyper_diffusion_exponent,
	bool force_switch,
	string field_name
):csf(field_name), rsf(field_name)
{			

	this->diffusion_coefficient		= diffusion_coefficient;
	this->hyper_diffusion_coefficient = hyper_diffusion_coefficient;
	this->hyper_diffusion_exponent = hyper_diffusion_exponent;
	
	if (hyper_diffusion_exponent == 0)
        hyper_diffusion_switch = false;

	this->force_switch = force_switch;
    
    nlin.resize(shape_complex_array);
	
	if (force_switch)
        Force.resize(shape_complex_array);
}


// rsf.F = inv_transform(csf.F)	
void FluidSF::Inverse_transform()
{
    rsf.Inverse_transform(csf); 

}


// csf.F = Forward_transform(rsf.F)	
void FluidSF::Forward_transform()
{
	csf.Forward_transform(rsf); 

}	

	
//******************************************************************************

void FluidSF::Satisfy_strong_reality_condition_field()
{
	universal->Satisfy_strong_reality_condition_in_Array(csf.F);
}

void FluidSF::Satisfy_strong_reality_condition_force_field()
{
	universal->Satisfy_strong_reality_condition_in_Array(Force);
}

void FluidSF::Satisfy_weak_reality_condition_field()
{
	universal->Satisfy_weak_reality_condition_in_Array(csf.F);
}

void FluidSF::Satisfy_weak_reality_condition_force_field()
{
	universal->Satisfy_weak_reality_condition_in_Array(Force);
}


void FluidSF::Test_reality_condition_field()
{
	universal->Test_reality_condition_in_Array(csf.F);
}

void FluidSF::Test_reality_condition_force_field()
{
	universal->Test_reality_condition_in_Array(Force);
}

//******************************************************************************

void FluidSF::Dealias_force_field()
{
	if (global.program.alias_option == "DEALIAS") 
		universal->Dealias(Force);
}

/**********************************************************************************************
 Copy_field_to(T):	  T.F <- F
 Copy_field_from(T):   csf.F <- T.F
***********************************************************************************************/

void FluidSF::Copy_field_to(PlainCSF& T)
{
	T.F = csf.F;
} 

void FluidSF::Copy_field_to(CSF& T)
{
	T.F = csf.F;
}

void FluidSF::Copy_field_from(PlainCSF& T)
{
	csf.F = T.F;
}  

void FluidSF::Copy_field_from(CSF& T)
{
	csf.F = T.F;
} 

/**********************************************************************************************
 
 Function 1:  csf.F = F + nlin*global.time.dt
 Function 2:  new CSF:  S.F = S.F + T.nlin*factor*global.time.dt
 
***********************************************************************************************/

void FluidSF::Add_nlin_factor_dt(Real factor)
{
	if (abs(factor) > MYEPS2)
		csf.F = csf.F + complex<Real>(factor*global.time.dt,0)*(nlin);
}

void FluidSF::Add_field_nlin_factor_dt(PlainCSF& S, Real factor)
{  
	if (abs(factor) > MYEPS2)						  
		S.F = S.F + complex<Real>(factor*global.time.dt,0) * (nlin);
}

//*********************************************************************************************
/** @brief Multiply scalar field by \f$ \exp(-\kappa K^2 global.time.dt) \f$ or 
 *			\f$ \exp(-\kappa K^2 global.time.dt) + \exp(-\kappa_h)  K^4 dt) \f$
 * 
 *  @param global.time.dt
 *
 * @return \f$ csf.F(k) = F(k) * \exp(-\kappa K^2 global.time.dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ csf.F(k) = F(k) * \exp(-\kappa K^2 global.time.dt) + \exp(-\kappa_h)  K^4 dt) \f$
 */
void FluidSF::Mult_field_exp_ksqr_dt(Real a)
{
    if (global.program.kind != "GP") {
        if (abs(a) > MYEPS) {
            if (!hyper_diffusion_switch){
                universal->Array_exp_ksqr(csf.F, -diffusion_coefficient*a*global.time.dt); 
        	}
            else {
                universal->Array_exp_ksqr(csf.F, -diffusion_coefficient*a*global.time.dt, -hyper_diffusion_coefficient*a*global.time.dt, hyper_diffusion_exponent);
            }
        }
    }
	
    else {
       universal->Array_exp_ksqr(csf.F, a);
    }
}

//*********************************************************************************************
/** @brief Multiply nlin of scalar field by \f$ \exp(-\kappa K^2 global.time.dt) \f$ 
 *			\f$ \exp(-\kappa K^2 global.time.dt) + \exp(-\kappa_h)  K^4 dt) \f$
 * 
 *  @param global.time.dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\kappa K^2 global.time.dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ csf.F(k) = F(k) * \exp(-\kappa K^2 global.time.dt) + \exp(-\kappa_h)  K^4 dt) \f$
 */
void FluidSF::Mult_nlin_exp_ksqr_dt(Real a)
{
    if (global.program.kind != "GP") {
        if (abs(a) > MYEPS) {
            if (!hyper_diffusion_switch) {					  
                universal->Array_exp_ksqr(nlin, -diffusion_coefficient*a*global.time.dt); 
        	}
            else {
                universal->Array_exp_ksqr(nlin, -diffusion_coefficient*a*global.time.dt, -hyper_diffusion_coefficient*a*global.time.dt, hyper_diffusion_exponent);
            }
        }
    }
    
    else {
      universal->Array_exp_ksqr(nlin, a);
    }
}


//******************************************************************************

/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once \f$ csf.F(kx, ky, kz) \f$ is given.
 *  
 * @bug  Add_complex_conj does not work across the processors.
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FFF: \f$ csf.F(-kx, -ky, kz) = F^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ csf.F(kx, -ky, kz) = F^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(csf.F(kx, 0, kz))  = 0\f$,
 *				
 */
void FluidSF::Add_complex_conj(int kx, int ky, int kz, Complex localG)					
{	
	// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FFF" || basis_type == "FFFW") && ((kz == 0) || (kz == N[3]/2)) )	{          
		universal->Assign_spectral_field(-kx, -ky, kz, csf.F, localG);
		if (master)
			cout << "Complex-conj(F) added for k = (" << -kx << "," << -ky << "," << kz << ")" << endl;	
	}
	
	else if ((basis_type == "SFF") && ((kz == 0) || (kz == N[3]/2)) ) 	{         
		universal->Assign_spectral_field(kx, -ky, kz, csf.F, localG);
		if (master)
			cout << "Complex-conj(F) added for k = (" << kx << "," << -ky << "," << kz << ")" << endl;	
	}
}


/** @brief D=3:  \f$ csf.F(kx, ky, kz) = G \f$ and add complex conjugate.
 *  
 * @sa Add_complex_conj(int kx, int ky,int kz,  Complex G)
 */
void FluidSF::Assign_field_and_comp_conj(int kx, int ky, int kz, Complex localG)
{
	universal->Assign_spectral_field(kx, ky, kz, csf.F, localG);
	
	Add_complex_conj(kx, ky, kz, localG);	
}


void FluidSF::Assign_field(int kx, int ky, int kz, Real localG)
{
	universal->Assign_spectral_field(kx, ky, kz, csf.F, localG);
}

/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range 
void FluidSF::Assign_random_complex_scalar(int kx, int ky, int kz, Real rand_range)
{
	Complex localG = Complex( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
	
	Assign_field_and_comp_conj(kx, ky, kz, localG);
}

/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range 
void FluidSF::Assign_random_real_scalar(int kx, int ky, int kz, Real rand_range)
{
	Real localG = 2*rand_range*(SPECrand.random()-0.5);
	
	Assign_field(kx, ky, kz, localG);
}	

/// 3D:  Return Tk = Real(-nlin(k). conj(csf.F(k)) for scalar T
Real FluidSF::Get_Tk(int kx, int ky, int kz)
{
	Complex localG_complex, localnlin_complex;
	Real localG_real, localnlin_real;
	
	if (universal->Probe_in_me(kx,ky,kz)) {
		if (basis_type == "SSS") {
			localG_real = real(universal->Get_spectral_field(kx, ky, kz, csf.F));
			localnlin_real = real(universal->Get_spectral_field(kx, ky, kz, nlin));
			return (localG_real*localnlin_real);
		}
		
		else {
			localG_complex = universal->Get_spectral_field(kx, ky, kz, csf.F);
			localnlin_complex = universal->Get_spectral_field(kx, ky, kz, nlin);
			return real_cprod(localG_complex, localnlin_complex);
		}
	}

	return 0;
}


Real FluidSF::Get_dt()
{
	
	if (global.program.dt_option == 0)
		return global.time.dt_fixed;
	
}


//************************ END of CVF class Definitions ***************************************
