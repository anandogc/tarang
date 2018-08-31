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

/*! \file  compute_flux.cc
 * 
 * @brief  Computes energy flux from inside/outside a sphere of field V/W to inside/outsude
 *			of field V/W.
 *
 *	The Giver field is filled inside/outside the sphere (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "EnergyTr.h"

//*********************************************************************************************

// Fluid
void EnergyTr::Compute_flux(FluidVF &U)
{	
	
	flux_self = 0.0;

	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		Fill_in_sphere(sphere_index, U);	
		
		Nlin_incompress::Compute_nlin(U, Giver);
			//	U.nlin = U.grad Giver<	
		
		flux_self(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, U);	
		// flux_self = -(U.grad U<). U> = U.nlin< . U>
	} 
}


//*********************************************************************************************

// Scalar
//

void EnergyTr::Compute_flux(FluidVF& U, FluidSF& T)
{
	if (global.program.kind == "INC_SCALAR")
		Compute_flux_scalar(U, T);
	
	else if (global.program.kind == "RBC" || global.program.kind == "STRATIFIED")
		Compute_flux_RBC(U, T);
}


void EnergyTr::Compute_flux_scalar(FluidVF& U, FluidSF& T)
{

	flux_SF = 0.0;
	
	Compute_flux(U);
	// flux_self = (U.grad U<). U>	
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_in_sphere(sphere_index, T);	
					
		Nlin_incompress::Compute_nlin_first_component(U, Giver);												
		// U.nlin1 = U.grad T<	
		
		flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
		// T.flux = (U.grad T<). T>		
							
	}
}

//
// RB convection
void EnergyTr::Compute_flux_RBC(FluidVF& U, FluidSF& T)
{
	
	flux_self = 0.0;
	flux_SF = 0.0;
	
	if (global.PHYSICS.Pr_option == "PRZERO") 
		Compute_flux(U);
	
	else if (global.PHYSICS.Pr_option == "PRINFTY") {		// fill only Temperature flux
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
			
			Fill_in_sphere(sphere_index, T);	
			
			Nlin_incompress::Compute_nlin_first_component(U, Giver);						
			// U.nlin1 = U.grad T<	
			
			flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
			// T.flux = -(U.grad T<). T>		
		}
	}
	
	else
		Compute_flux_scalar(U, T);	
}


//*********************************************************************************************

//
// Vector
void EnergyTr::Compute_flux(FluidVF& U, FluidVF& W)
{
	Compute_flux(U);									
	// flux_self = (U.grad U<). U>
	
	flux_VF_Uin_Wout = 0.0;
	flux_VF_Uin_Win = 0.0;
	
	// U< to W>;  U< to W<
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_in_sphere(sphere_index, U);				
		
		Nlin_incompress::Compute_nlin(W, Giver);												
		// W.nlin = W.grad U<
		
		flux_VF_Uin_Wout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W, W);	
		// (W.graad U<). W>
		
		flux_VF_Uin_Win(sphere_index) = Prod_in_sphere_nlinV(sphere_index, W, W);		
		// (W.graad U<). W<
	}
	
	// U> to W>
	flux_VF_Uout_Wout = 0.0;
    flux_VF_Uout_Win = 0.0;//Abhishek 6th June 2017
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {		
		Fill_out_sphere(sphere_index, U);										
			// G = U>
		
		Nlin_incompress::Compute_nlin(W, Giver);											
			// W.nlin = W.grad U>
		
		flux_VF_Uout_Wout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W, W);	
			// (W.graad U>). W>
        flux_VF_Uout_Win(sphere_index) = Prod_in_sphere_nlinV(sphere_index, W, W);
        //  (W.graad U>). W<… This is  Pi U> to  W<
	}
	
	
	// W< to W>;  W< to U>
	flux_VF_Win_Wout = 0.0;
	//flux_VF_Win_Uout = 0.0;
  
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {	
		Fill_in_sphere(sphere_index, W);										
		// G = W<
		
		Nlin_incompress::Compute_nlin(U, Giver);												
		// U.nlin = U.grad W<
		
		flux_VF_Win_Wout(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, W);		
		// -(U.graad W<). W>
		
		//flux_VF_Win_Uout(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, U);
		//  (U.graad W<). U>
      //(U is nlin argument, U>)
	}
	
	//
	// Flux for Elsasser vars
	//
	
	if (global.energy_transfer.Elsasser) {
		flux_Elsasser_plus = 0.0;
		flux_Elsasser_minus = 0.0;
		
		
		MHD::UB_to_Elsasser_field(U, W);												
		// U=Zp=(U+B); B=Zm=(U-B);
		
		// Flux: Zp to Zp
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {		
			Fill_in_sphere(sphere_index, U);											
			// G = Zp<
			
			Nlin_incompress::Compute_nlin(W, Giver);												
			// W.nlin = Zm.grad Zp<
			
			flux_Elsasser_plus(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, W, U);	
			// (Zm.graad Zp<). Zp>
		
		}	
		
		
		// Flux: Zm to Zm
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {
			Fill_in_sphere(sphere_index, W);										
			// G = Zm<
			
			Nlin_incompress::Compute_nlin(U, Giver);												
			// U.nlin = Zp.grad Zm<
			
			flux_Elsasser_minus(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, W);	
			// (Zp.graad Zm<). Zm>
		
		}																																			// Flux: Zp to Zp
		MHD::Elsasser_to_UB_field(U, W);													
		// Back to U, B vars
	} // end of if(Elsasser)
}


//*********************************************************************************************
//
// Magnetoconvection

void EnergyTr::Compute_flux(FluidVF& U, FluidVF& W, FluidSF& T)
{
	
	flux_SF = 0.0;
	Compute_flux(U, W);
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++)  {		
		Fill_in_sphere(sphere_index, T);	
					
		Nlin_incompress::Compute_nlin_first_component(U, Giver);											
		// U.nlin = U.grad T<	
		
		flux_SF(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, T);		
		// T.flux = (U.grad T<). T>		
	}
}
/*********************************************************************************************

 HELICITY FLUX etc..
***********************************************************************************************/

// Fluid

void EnergyTr::Compute_kinetic_helicity_flux(FluidVF& U, FluidVF& helicalU)
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
  
  
  
  flux_VF = 0.0;
  
  
  for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
    
    Fill_in_sphere(sphere_index, helicalU);
    // nlin = U.grad B<
    Nlin_incompress::Compute_nlin_helical(U, Giver);
    flux_VF(sphere_index) = Prod_out_sphere_nlinV(sphere_index, U, helicalU);
  }
  
}


/*void EnergyTr::Compute_kinetic_helicity_flux_old(FluidVF& U, FluidVF& helicalU)
{
	universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());

	flux_VF_Uin_Uout = 0.0;
	flux_VF_Uin_Wout = 0.0;
	flux_VF_Win_Uout = 0.0;

	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
		// (helicalU.graad U<). U>
		Fill_in_sphere(sphere_index, U);
		// helicalU.nlin = helicalU.grad U<
		Nlin_incompress::Compute_nlin(helicalU, Giver);
		flux_VF_Uin_Uout(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, helicalU, U);


		// (U.graad U<). helicalU>
		Fill_in_sphere(sphere_index, U);
		// U.nlin = U.grad U<
		Nlin_incompress::Compute_nlin(U, Giver);
		flux_VF_Uin_Wout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, helicalU);


		// (U.graad helicalU<). U>
		Fill_in_sphere(sphere_index, helicalU);
		// U.nlin = U.grad helicalU<
		Nlin_incompress::Compute_nlin(U, Giver);
		flux_VF_Win_Uout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, U);
	}
}*/

//**************************** Kinetic Helicity in MHD ****************************************
// Satyajit Franck Abhishek ** 21st May. 2017.

void EnergyTr::Compute_kinetic_helicity_flux(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW)//Rename this to old
{
  universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
  universal->Compute_vorticity(W.cvf.V1, W.cvf.V2, W.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
  
  flux_VF_U = 0.0;
  flux_VF_B = 0.0;
  
  
  for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {
    
    Fill_in_sphere(sphere_index, helicalU);
    // nlin = U.grad B<
    Nlin_incompress::Compute_nlin_helical(U, Giver);
    flux_VF_U(sphere_index) = Prod_out_sphere_nlinV(sphere_index, U, helicalU);
    
    Fill_in_sphere(sphere_index, helicalW);
    // nlin = U.grad B<
    Nlin_incompress::Compute_nlin_helical(W, Giver);
    flux_VF_B(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, W, helicalU);
    

  }  
}


/*void EnergyTr::Compute_kinetic_helicity_flux(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW)//Rename this to old
{
	universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
	universal->Compute_vorticity(W.cvf.V1, W.cvf.V2, W.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());


	flux_VF_Uin_Uout = 0.0;
	flux_VF_Uin_Wout = 0.0;
	flux_VF_Bin_Uout = 0.0;
	flux_VF_Bin_Wout = 0.0;
	flux_VF_Win_Uout = 0.0;
	flux_VF_Jin_Uout = 0.0;

	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {

		// (U.graad U<). helicalU>
		Fill_in_sphere(sphere_index, U);
		// U.nlin = U.grad U<
		Nlin_incompress::Compute_nlin(U, Giver);
		flux_VF_Uin_Wout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, helicalU);



		// (B.graad B<). helicalU>
		Fill_in_sphere(sphere_index, W);
		// nlin = B.grad B<
		Nlin_incompress::Compute_nlin(W, Giver);
		flux_VF_Bin_Wout(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, W, helicalU);



		// (U.graad helicalU<). U>
		Fill_in_sphere(sphere_index, helicalU);
		// U.nlin = U.grad helicalU<
		Nlin_incompress::Compute_nlin(U, Giver);
		flux_VF_Win_Uout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, U);


		// (helicalU.graad U<). U>
		Fill_in_sphere(sphere_index, U);
		// helicalU.nlin = helicalU.grad U<
		Nlin_incompress::Compute_nlin(helicalU, Giver);
		flux_VF_Uin_Uout(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, helicalU, U);

		// (J.graad B<). U>
		Fill_in_sphere(sphere_index, W);
		// nlin = J.grad B<
		Nlin_incompress::Compute_nlin(helicalW, Giver);
		flux_VF_Bin_Uout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, helicalW, U);



		// (B.graad J<). U>
		Fill_in_sphere(sphere_index, helicalW);
		// nlin = W.grad J<
		Nlin_incompress::Compute_nlin(W, Giver);
		flux_VF_Jin_Uout(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, W, U);
	}  
}*/

//*********************************************************************************************

//**************************** Magnetic Helicity in MHD ****************************************
// Satyajit and Abhishek ** 14th April 2017.

void EnergyTr::Compute_magnetic_helicity_flux(FluidVF& U, FluidVF& W)
{

  flux_HM = 0.0;
  
  
  for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {

    Fill_in_sphere(sphere_index, W);
    // nlin = U.grad B<
    Nlin_incompress::Compute_nlin_helical(U, Giver);
    flux_HM(sphere_index) = Prod_out_sphere_nlinV(sphere_index, U, W);
  }  
  
  
}

void EnergyTr::Compute_magnetic_helicity_flux(FluidVF& U, FluidVF& W, FluidVF& helicalW)//Old
{
	// New part added on 13th Feb. 2017
	// Here helicalW --> Magnetic vector field(A)
	universal->Compute_vorticity(W.cvf.V1, W.cvf.V2, W.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
	
	helicalW.cvf.Divide_ksqr();
	
	// First we need to compute "A" from an universal function defined properly in specific place. Do it here!!!!!!!
	
	// Initializations of variables, defined and resized properly
	flux_VF_Bin_Aout = 0.0;
	flux_VF_Uin_Aout_1 = 0.0;
	flux_VF_Ain_Bout = 0.0;
	flux_VF_Uin_Aout_2 = 0.0;
	
	for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {

		// (U.grad B<). helicalW>
		Fill_in_sphere(sphere_index, W);
		// nlin = U.grad B<
		Nlin_incompress::Compute_nlin(U, Giver);
		flux_VF_Bin_Aout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, helicalW);



		// (B.graad U<). helicalW>
		Fill_in_sphere(sphere_index, U);
		// nlin = B.grad U<
		Nlin_incompress::Compute_nlin(W, Giver);
		flux_VF_Uin_Aout_1(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, W, helicalW);
		
		// (U.graad A<). W>
		Fill_in_sphere(sphere_index, helicalW);
		// nlin = U.grad A<
		Nlin_incompress::Compute_nlin(U, Giver);
		

				
		
		flux_VF_Ain_Bout(sphere_index) = -0.5 * Prod_out_sphere_nlinV(sphere_index, U, W);
		
		
       // A little different nonlinear term //
       // Appropriate nonlinear calculator should be called //
       
		// (U_j.del_iA_j<). W>
		Fill_in_sphere(sphere_index, helicalW);
		// nlin = U_j.del_iA_j<
		Nlin_incompress::Compute_nlin_vector_potential(U, Giver); // Allready exist or not

				
		flux_VF_Uin_Aout_2(sphere_index) = 0.5 * Prod_out_sphere_nlinV(sphere_index, U, W); // Question regarding the arguments
		
	}  
	

}



//*********************************************************************************************

// Fluid 2D
/*void EnergyTr::Compute_enstrophy_flux(FluidVF& U)
{	

	if (N[2] == 1)
	{
		flux_hk = 0.0;
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere_vorticity(sphere_index, U);
			// G = i k x u(k): 2D
			
			Compute_nlin(U);										
			// U.nlin = U.grad omega	
					
			flux_hk(sphere_index) =   - Prod_out_sphere_nlin_vorticity(sphere_index, U, U);	
			// -(U.grad omega). omega> : 2D
			// Can make this efficient by mult omega_y nlin_y
		}
	}
}*/

void EnergyTr::Compute_enstrophy_flux(FluidVF& U)
{
	
	//Compute_flux(U);									
	// flux_self = (U.grad U<). U>
	if (Ny>1) {
		FluidVF helicalU("helicalU");

		universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());

		cout << "max abs U = " << max(abs(U.cvf.V1)) << " " << max(abs(U.cvf.V2)) << " " << max(abs(U.cvf.V3)) << endl;
		cout << "max abs W = " << max(abs(helicalU.cvf.V1)) << " " << max(abs(helicalU.cvf.V2)) << " " << max(abs(helicalU.cvf.V3)) << endl;
		
		flux_VF_Win_Wout = 0.0;
		flux_VF_Uin_Wout = 0.0;
		// U< to helicalU>;  U< to helicalU<
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {	

			Fill_in_sphere(sphere_index, helicalU);
			
			Nlin_incompress::Compute_nlin(U, Giver);
			// U.nlin = U.grad helicalU<


			flux_VF_Win_Wout(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, U, helicalU);	
			// (U.graad helicalU<). helicalU>

			 Fill_in_sphere(sphere_index, U);

			 Nlin_incompress::Compute_nlin(helicalU, Giver);												
			//  helicalU.nlin = helicalU.grad U<
		
			 flux_VF_Uin_Wout(sphere_index) = Prod_out_sphere_nlinV(sphere_index, helicalU, helicalU);
			// (helicalU.graad U<). helicalU>
		  
		
		}
	}
	else {
	
		flux_VF_Win_Wout = 0.0;
		flux_VF_Uin_Wout = 0.0;
		
		FluidVF helicalU("helicalU");
		FluidSF helicalUT("helicalUT");

		universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
		
		helicalUT.csf.F =  helicalU.cvf.V2;
		
		// flux_self = (U.grad U<). U>	
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) {	

			Fill_in_sphere(sphere_index, helicalUT);	
						
			Nlin_incompress::Compute_nlin_first_component(U, Giver);												
			// U.nlin1 = U.grad helicalU<


			flux_VF_Win_Wout(sphere_index) = -Prod_out_sphere_nlinT(sphere_index, U, helicalUT);	
			// (U.graad helicalU<). helicalU>
			
			
		}
	}
}


//*********************************************************************************************

// MHD 2D
void EnergyTr::Compute_magnetic_enstrophy_flux(FluidVF& U, FluidVF& W)
{	
	/*
	if (N[2] == 1)
	{
		flux_Whk = 0.0;
		
		for (int sphere_index = 1; sphere_index <= global.energy_transfer.flux.no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere(sphere_index, W);
			// G = W<
			
			Compute_nlin_UcrossB(U);
			// U.nlin = u x W< : 2D
			
			flux_Whk(sphere_index) = Prod_out_sphere_nlin_vector_potential(sphere_index, U, W);	
			// add (U x W<). B>
		}
	}
	 */
}

		

//****************************  End of compute_flux.cc ****************************************







