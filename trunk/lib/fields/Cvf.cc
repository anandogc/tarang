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

#include "Cvf.h"
#include "Rvf.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
 
**********************************************************************************************/

CVF::CVF(string field_name)
{
	this->field_name=field_name;
	
	total_energy = total_k2energy = entropy = 0.0;
}

/**********************************************************************************************

		   Copies Arrays Vi's to to.Vi

**********************************************************************************************/

void CVF::Copy_to(CVF& to)
{
	to.V1 = V1;   
	to.V2 = V2;   
	to.V3 = V3;
}

void CVF::Copy_to(PlainCVF& to)
{
	to.V1 = V1;   
	to.V2 = V2;   
	to.V3 = V3;
}

void CVF::Copy_from(CVF& from)
{
	V1 = from.V1;   
	V2 = from.V2;   
	V3 = from.V3;
}

void CVF::Copy_from(PlainCVF& from)
{
	V1 = from.V1;   
	V2 = from.V2;   
	V3 = from.V3;
}

		
/**********************************************************************************************

      Inplace Forward transform of all the components   
	 -- temp_r is temporary location 

**********************************************************************************************/


/* void CVF::Forward_transform()
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Forward_transform_array(V1);			
	
	if (!global.program.two_dimension) { // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Forward_transform_array(V2);		
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Forward_transform_array(V3);			
}
 */


/**********************************************************************************************

      Forward_transform_transpose_order(Vitr) = Vi 
		Vir unchanged
		temp_r = N2  x N1 x N3

**********************************************************************************************/


void CVF::Forward_transform(RVF rvf)
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Forward_transform(rvf.V1r, V1);
	
	if (!global.program.two_dimension) { // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Forward_transform(rvf.V2r, V2);
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Forward_transform(rvf.V3r, V3);
}



/**********************************************************************************************

       Inplace Inverse transform of all the components   
	 -- temp_r is temporary space

**********************************************************************************************/


/*void CVF::Inverse_transform()
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Inverse_transform_array(V1);			
	
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Inverse_transform_array(V2);	
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Inverse_transform_array(V3);			
}*/



/**********************************************************************************************

       Inverse_transform(*Vi) = Vitr 
		*Vi unchanged
		temp = N1 x N2 x N3

**********************************************************************************************/


void CVF::Inverse_transform(RVF rvf)
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Inverse_transform(V1, rvf.V1r);
	
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Inverse_transform(V2, rvf.V2r);
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Inverse_transform(V3, rvf.V3r);
}


/**********************************************************************************************
 
 F(k) = F(k)/k^2
 
 **********************************************************************************************/

void CVF::Divide_ksqr()
{
	
	universal->Array_divide_ksqr(V1); 
	universal->Array_divide_ksqr(V2); 
	universal->Array_divide_ksqr(V3);
	
}


/**********************************************************************************************

       Computes the total energy sum[|V(k)|^2/2]
			 total dissipation sum[k^2|V(k)|^2]
			 Entropy

**********************************************************************************************/

void CVF::Compute_total_energy()
{

	total_E1 = universal->Get_total_energy(V1);
	total_E2 = universal->Get_total_energy(V2);
	total_E3 = universal->Get_total_energy(V3);
	
	total_energy = total_E1 + total_E2 + total_E3;
}


void CVF::Compute_total_k2energy()
{
	total_k2energy = TWO*(universal->Get_total_Sn(V1, 2)+ universal->Get_total_Sn(V2, 2) + universal->Get_total_Sn(V3, 2));
}

void CVF::Compute_total_kn_energy(int n, Real &result)
{
	result = 2*(universal->Get_total_Sn(V1, n)+ universal->Get_total_Sn(V2, n) + universal->Get_total_Sn(V3, n));		
}


void CVF::Compute_total_k2Hc(Array<Complex,3> W1, Array<Complex,3> W2, Array<Complex,3> W3)
{
	total_k2Hc = TWO*(universal->Get_total_Sn(V1, W1, 2) + universal->Get_total_Sn(V2, W2, 2) + universal->Get_total_Sn(V3, W3, 2));
}


//
//

void CVF::Compute_entropy()
{
	 entropy = universal->Get_entropy(V1, V2, V3);
			
}


/**********************************************************************************************

       Computes moaal energy, modal helicity, modal vorticity.


**********************************************************************************************/



Real CVF::Modal_energy(int l1, int l2, int l3)
{
	return (universal->Modal_energy(l1, l2, l3, V1) + universal->Modal_energy(l1, l2, l3, V2) + universal->Modal_energy(l1, l2, l3, V3));
}				

//
//

void CVF::Compute_Modal_vorticity(int l1, int l2, int l3, TinyVector<Complex,3> &vorticity)
{
	universal->Compute_Modal_vorticity(l1, l2, l3, V1, V2, V3, vorticity);
}

void CVF::Compute_Modal_vorticity_y_component(int l1, int l2, int l3, Complex &vort_y)
{
	universal->Compute_Modal_vorticity_y_component(l1, l2, l3, V1, V2, V3, vort_y);
}

//
//

Real CVF::Modal_helicity(int l1, int l2, int l3)
{
	return universal->Get_Modal_helicity(l1, l2, l3, V1, V2, V3);
}

/**********************************************************************************************

       Computes total helicity, and helicity spectra


**********************************************************************************************/


void CVF::Compute_total_helicity()
{

	universal->Compute_total_helicity(V1, V2, V3,total_helicity1, total_helicity2,total_k2H1, total_k2H2);
										
}


//*********************************************************************************************

/*
void CVF::Compute_structure_function()
{
	
	if (structure_fn_switch == 1)	
		Compute_structure_function(basis_type, Ncv, *V1, *V2, *V3, structurefn_q_min, 
								   structurefn_q_max, *St, xfactor);		
}


void CVF::Compute_planar_structure_function()
{
	if (planar_structure_fn_switch == 1)
		Compute_planar_structure_function(basis_type, Ncv, *V1, *V2, *V3, 
									structurefn_q_min, structurefn_q_max, 
									*st_planar, xfactor);	
																			
}
 */

//******************************************************************************
/** @brief Put a vector at (lx,0,lz) given the amp and phase V_c1 (Craya-Herring)
 *  For 2D
 *
 *  @param lx, lz  (ly=0): location where the V is to be assigned.
 *	@param amp		Amplitude of the vector V_c1 (Craya-Herring).
 *	@param phase Phase of the vector V_c1 (Craya-Herring)
 *  @param add_flag: if yes,it adds the vector to existing V, otherwise assigns V at (lx,0,lz)
 *
 *  @detail If k=(0,0,0), set Force=(0,0,0)
 *  @detail For 2D, (xz) plane of the code is to be treaed as (x,y) of Craya-Herring decomposition.  For 2D, the anisotropy direction is ignored, and Vx = V_c1 sin(phi), Vz= -V_c1 cos(phi).
 *  @detail: (0,0,0) assigned before...while calling Init_cond_energy_helicity_spectrum()
 */
void CVF::Put_or_add_vector(int lx, int lz,  Real amp, Real phase, bool add_flag)
{
    Complex V_ch1;  // Amplitude along Craya-herring basis vector e1.
    Real kx, kz, Kmag, phi;
    
    TinyVector<Complex,3> VFour, Vlocal_complex;
    TinyVector<Real,3> Vlocal_real;
    
    // k != (0,0,0) : PS:(0,0,0) assigned before...while calling Init_cond_energy_helicity_spectrum()
    Kmag  = universal->Kmagnitude(lx, 0, lz);
    if (Kmag > MYEPS) {
        V1 = amp * exp(I*phase);
        kx = universal->Get_kx(lx) * kfactor[1];
        kz = universal->Get_kz(lz) * kfactor[3];
        phi = Get_azimuthal_angle(kx, kz);
        
        VFour = V_ch1*sin(phi), 0, -V_ch1*cos(phi);
        
        if (basis_type == "FFF" || basis_type == "FFFW") {
            Vlocal_complex = VFour;
            
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_complex);
            else
                universal->Add_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_complex);
        }
        
        else if (basis_type == "SSS") {
            Convert_from_Fourier_space(VFour, Vlocal_real);
            
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_real);
            else
                universal->Add_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_real);
        }
        
        else {
            Convert_from_Fourier_space(VFour, Vlocal_complex);
            
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_complex);
            else
                universal->Add_local_spectral_field(lx, 0, lz, V1, V2, V3, Vlocal_complex);
        }
    }
}




void CVF::Put_or_add_vector(int lx, int ly, int lz,  Real amp_plus, Real amp_minus, Real phase_plus, Real phase_minus, bool add_flag)
{
    Complex u_plus, u_minus;   // Amplitudes along helical basis
    Complex u_ch1, u_ch2;     // Amplitudes along Craya-Herring basis
    Complex vpll, vh1, vh2;
    Real kx,kz,phi,Kmag;
    
    TinyVector<Complex,3> VFour, Vlocal_complex;
    TinyVector<Real,3> Vlocal_real;
 
    // k != (0,0,0) : PS:(0,0,0) assigned before...while calling Init_cond_energy_helicity_spectrum()
    
    Kmag  = universal->Kmagnitude(lx, ly, lz);
    if (Kmag > MYEPS) {
        u_plus = amp_plus * exp(I*phase_plus);
        u_minus = amp_minus * exp(I*phase_minus);
        universal->Helical_to_Craya(u_plus,u_minus,u_ch1,u_ch2);
        
        if (Ny == 1) {  // 2D3C
            kx = universal->Get_kx(lx) * kfactor[1];
            kz = universal->Get_kz(lz) * kfactor[3];
            phi = Get_azimuthal_angle(kx, kz);  // tan^{-1}(kz/kx)
            
            VFour = u_ch1*sin(phi), u_ch2, -u_ch1*cos(phi);
        }
        else { // 3D
            universal->Craya_to_cartesian(lx, ly, lz, u_ch1, u_ch2, vpll, vh1, vh2);
            // Gets components pll, fh1, fh2 in Cartesian basis.
            
            if (global.field.anisotropy_dirn == 1)
                VFour = vpll, vh1, vh2;
            
            else if (global.field.anisotropy_dirn == 2)
                VFour = vh2, vpll, vh1;
            
            else if (global.field.anisotropy_dirn == 3)
                VFour = vh1, vh2, vpll;
        }
    
        
        if (basis_type == "FFF" || basis_type == "FFFW") {
            Vlocal_complex = VFour;
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_complex);
            else
                universal->Add_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_complex);
        }
        
        else if (basis_type == "SSS") {
            Convert_from_Fourier_space(VFour, Vlocal_real);
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_real);
            else
                universal->Add_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_real);
        }
    
        
        else {
            Convert_from_Fourier_space(VFour, Vlocal_complex);
            if (!add_flag)
                universal->Assign_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_complex);
            else
                universal->Add_local_spectral_field(lx, ly, lz, V1, V2, V3, Vlocal_complex);
        }
    }
}

/**********************************************************************************************

       Dealias CVF


**********************************************************************************************/


void CVF::Dealias()
{
	if (global.program.alias_option == "DEALIAS") {
		universal->Dealias(V1);
		universal->Dealias(V2);
		universal->Dealias(V3);
	}
	
}


/**********************************************************************************************

      Input and Output of all the components (Arrays)
 
 File name:  
 dir=time 
 
 Filenames: V1, V2, V3, same as fields of CVF, RVF, CSF etc.
 Reduced: V1reduced, V2reduced
 If plane V3kz0
 
 Do not write field-frequent and regular file at the same time..
 Overwrite field frequent, but not the regular..
 May be just add frequent label with time in the dirname, i.e., frquent.time
 
use IOdir as global variable and control...


**********************************************************************************************/

// field_name = U or W
void CVF::Write_complex_field()
{
	string folder_name = "time_" + To_string(global.time.now);

    universal->Write(V1, universal->H5_full, "w", folder_name, field_name+".V1");
	
    if (!global.program.two_dimension)
        universal->Write(V2, universal->H5_full, "w", folder_name, field_name+".V2");
	
    if (global.field.incompressible && global.io.output_vx_vy_switch)
        universal->Write(V3, universal->H5_kz0_full, "w", folder_name, field_name+".V3kz0");
    else
        universal->Write(V3, universal->H5_full, "w", folder_name, field_name+".V3");
}



void CVF::Write_reduced_complex_field()
{
	string folder_name="reduced_" + To_string(global.time.now);

    universal->Write(V1, universal->H5_out_reduced, "w", folder_name, field_name+".V1");
    
    if (!global.program.two_dimension)
        universal->Write(V2, universal->H5_out_reduced, "w", folder_name, field_name+".V2");
	
	if (global.field.incompressible && global.io.output_vx_vy_switch)
        universal->Write(V3, universal->H5_out_kz0_reduced, "w", folder_name, field_name+".V3kz0");
    else
        universal->Write(V3, universal->H5_out_reduced, "w", folder_name, field_name+".V3");
}


//
void CVF::Read_complex_field()
{
	universal->Read(V1, universal->H5_full, field_name+".V1");

    if (!global.program.two_dimension)
        universal->Read(V2, universal->H5_full, field_name+".V2");

	if (global.field.incompressible && global.io.input_vx_vy_switch)
        universal->Read(V3, universal->H5_kz0_full, field_name+".V3kz0");
    else
        universal->Read(V3, universal->H5_full, field_name+".V3");
}

void CVF::Read_reduced_complex_field()
{
    V1 = 0.0;
    V2 = 0.0;
    V3 = 0.0;
    
	universal->Read(V1, universal->H5_in_reduced, field_name+".V1");
    
    if (!global.program.two_dimension)
        universal->Read(V2, universal->H5_in_reduced, field_name+".V2");
	
	if (global.field.incompressible && global.io.input_vx_vy_switch)
        universal->Read(V3, universal->H5_in_kz0_reduced, field_name+".V3kz0");
    else
        universal->Read(V3, universal->H5_in_reduced, field_name+".V3");
}


//************************ END of CVF class Definitions ***************************************
