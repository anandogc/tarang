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

/*! \file  ic_model.cc
 * 
 * @brief Initial condition where low k modes are forced: TG, ABC flows etc.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug	 No known bugs
 */


#include "FluidIO.h"



//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x) sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void FluidIO::Setup_Taylor_Green_field(FluidVF& U, int k0, Real amp)
{
	U.cvf.V1 = 0.0; 
	U.cvf.V2 = 0.0;  
	U.cvf.V3 = 0.0;		
	
	// CHOOSE SCC type
	if (basis_type == "SSS") {
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V1, amp/8);
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V2, -amp/8);
	}
	
	// CHOOSE SCF type
	else if (basis_type == "SSF") {
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V1, amp/8);
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V2, -amp/8);
	}
	
	// CHOOSE SFF type
	else if (basis_type == "SFF") {

		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V1, Complex(amp/8, 0.0));
		universal->Assign_spectral_field(k0, -k0, k0, U.cvf.V1, Complex(amp/8, 0.0));
		
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V2, Complex(0, amp/8));
		universal->Assign_spectral_field(k0, -k0, k0, U.cvf.V2, Complex(0, -amp/8));
	}
	
	else if (basis_type == "FFF" || basis_type == "FFFW") {
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V1, Complex(0, -amp/8));
		universal->Assign_spectral_field(k0, -k0, k0, U.cvf.V1, Complex(0, -amp/8));
		universal->Assign_spectral_field(-k0, k0, k0, U.cvf.V1, Complex(0, amp/8));
		universal->Assign_spectral_field(-k0, -k0, k0, U.cvf.V1, Complex(0, amp/8));
		
		universal->Assign_spectral_field(k0, k0, k0, U.cvf.V2, Complex(0, amp/8));
		universal->Assign_spectral_field(-k0, k0, k0, U.cvf.V2, Complex(0, amp/8));
		universal->Assign_spectral_field(k0, -k0, k0, U.cvf.V2, Complex(0, -amp/8));
		universal->Assign_spectral_field(-k0, -k0, k0, U.cvf.V2, Complex(0, -amp/8));
	}
}

//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @param k0
 *	@param amp  
 *
 */


void FluidIO::Setup_ABC_field(FluidVF& U, int k0, Real amp, Real A, Real B, Real C)
{

	U.cvf.V1 = 0.0; 
	U.cvf.V2 = 0.0;  
	U.cvf.V3 = 0.0;			
	
	if (basis_type == "FFF" || basis_type == "FFFW") {
		universal->Assign_spectral_field(0, k0, 0, U.cvf.V1, Complex(amp*B/2, 0));
		universal->Assign_spectral_field(0, -k0, 0, U.cvf.V1, Complex(amp*B/2, 0));
		universal->Assign_spectral_field(0, 0, k0, U.cvf.V1, Complex(0, -amp*C));
		
		universal->Assign_spectral_field(k0, 0, 0, U.cvf.V2, Complex(0, -amp*A/2));
		universal->Assign_spectral_field(-k0, 0, 0, U.cvf.V2, Complex(0, amp*A/2));
		universal->Assign_spectral_field(0, 0, k0, U.cvf.V2, Complex(amp*C, 0));
		
		universal->Assign_spectral_field(k0, 0, 0, U.cvf.V3, Complex(amp*A/2, 0));
		universal->Assign_spectral_field(-k0, 0, 0, U.cvf.V3, Complex(amp*A/2, 0));
		universal->Assign_spectral_field(0, k0, 0, U.cvf.V3, Complex(0, -amp*B));
	}
	
	else {
		cout << "ABC field not allowed in the used basis " << endl;
		exit(1);
	}
}



//*********************************************************************************************

	// WORK ON IT...
void FluidIO::Setup_SIX_MODE_field(FluidVF& U, int k0, Real amp101, Real amp011, Real amp112, Real h)
{
/*
	V1 = 0.0; 
	V2 = 0.0;  
	V3 = 0.0; 
	// initialize

	if (basis_type == "FFF" || basis_type == "FFFW") {
		int lx_k0 = Get_lx(k0);
		int lx_minus_k0 = Get_lx(-k0);
		int ly_k0 = Get_ly3D(k0);
		int ly_minus_k0 = Get_ly3D(-k0);
		
		Real factor1 = 2.0/sqrt(2+h*h) * amp101;
		Real factor2 = 2.0/sqrt(2+h*h) * amp011;
		Real factor3 = 4.0/sqrt(6+10*h*h) * amp112;
		
		if (my_id == master_id)  
		{
			(V1)(0, ly_k0, k0) = Complex(-h*(factor2/4), 0.0);
			(V1)(0, ly_minus_k0, k0) = Complex(h*(factor2/4), 0.0);
			
			(V2)(0, ly_k0, k0)  = (I) * (factor2/4);
			(V2)(0, ly_minus_k0, k0) = (-I) * (factor2/4);
			
			(V3)(0, ly_k0, k0)  = (-I) * (factor2/4);
			(V3)(0, ly_minus_k0, k0) = (-I) * (factor2/4); 
		}
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) )
		{
			(V1)(lx_k0, 0, k0)  = I * (factor1/4);
			(V2)(lx_k0, 0, k0)  = Complex(-h*(factor1/4), 0.0);
			(V3)(lx_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(V1)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8);
			(V1)(lx_k0, ly_minus_k0, 2*k0) = I * (factor3/8);
			
			(V2)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8) + Complex(h*factor3/4, 0.0);
			(V2)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + Complex(h*factor3/4, 0.0);
			
			(V3)(lx_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) - Complex(h*factor3/8, 0.0);
			(V3)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + Complex(h*factor3/8, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(V1)(lx_minus_k0, 0, k0)  = -I * (factor1/4);
			(V2)(lx_minus_k0, 0, k0)  = Complex(h*(factor1/4), 0.0);
			(V3)(lx_minus_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(V1)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8);
			(V1)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8);	
					
			(V2)(lx_minus_k0, ly_k0, 2*k0)  = I * (factor3/8) - Complex(h*factor3/4, 0.0);
			(V2)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- Complex(h*factor3/4, 0.0);
			
			(V3)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) + Complex(h*factor3/8, 0.0);
			(V3)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- Complex(h*factor3/8, 0.0);
		}	
		
	}
	
	else if ((basis_type == "SCFT") && (my_id == master_id))
	{
		cout << "Setup_SIX_MODE_field Not allowed in SCFT basis_type " << endl;
		exit(0);
	}
*/
}


/**********************************************************************************************
 
 Based on Pope's book Turbulent flow, page 232.
 
 Model energy spectrum for using Pope's formula
 Sk(k) = C epsilon^2/3 k^-5/3 f_L(kL) f_eta(k eta)
 f_L = [(kL)/sqrt{(kL)^2 + c_L}]^{5/3+p_0}  
 f_eta = exp[ -beta {[(k eta)^4+c_eta^4]^1/4 - c_eta} ]
 c_L = 6.78, c_eta = 0.40, beta = 5.2, p0=2
 C = 1.5; Kolmogorov's constant
 
 Implemented on 9 Oct 2011.
 
 ***********************************************************************************************/


void FluidIO::Model_initial_using_shell_spectrum_Pope(Real dissipation_coefficient, Real epsilon, Array<Real,1> Sk)
{
	
	Real k;	// radius
	Real L;	// size of the box
	Real fL, feta, fL_arg, feta_arg;
	Real p0 = 2;
	Real cL = 6.78;
	Real ceta = 0.40;
	Real beta = 5.2;
	Real Kolmogrov_const = 1.5;
	
	Sk = 0.0;
	
	if (N[2] > 1) 
		L =  pow(global.field.L[1]*global.field.L[2]*global.field.L[3], (Real) 1./3);  
	else if (N[2] == 1)
		L =  sqrt(global.field.L[1]*global.field.L[3]);
		
	Real eta = 0.42*pow(my_pow(dissipation_coefficient,3)/epsilon, (Real) 1./4);
	// Kolmogrov's length.. The factor 0.42 to make Pi(k)-> 0 as k-> infty.
	// Verma's on LM
	
	Sk(0) = 0.0;
	for (int i=1; i < (Sk.length())(0); i++) {
		k = ((Real) i);
		fL_arg = k*L/sqrt(my_pow(k*L,2)+cL*cL);
		fL = pow(fL_arg, (Real) 5.0/3+p0);
		feta_arg = pow( my_pow(k*eta,4) + my_pow(ceta,4), (Real) 1.0/4) - ceta;
		feta = exp(-beta *feta_arg);
		
		Sk(i) = Kolmogrov_const*pow(epsilon, (Real) 1.0/3)*pow(k, (Real) -5.0/3)*fL*feta;
	}
}


void FluidIO::Model_initial_using_shell_spectrum_Corrsin(Real a, Array<Real,1> Sk)
{
	Real q = 1.5;
	Real b = 0.02;
	Real c = 1+2.8/12;
	Real k;	// radius

	Sk = 0.0;

	for (int i=0; i < (Sk.length())(0); i++) {
		k = 1.0*i;
		Sk(i) = a * my_pow(k,4) * exp(- b * pow(k, (Real) 1.1)) /pow( (my_pow(k,4)+my_pow(q,4)), (Real) c);
    }	
}



//*******************************  End of ic_model.cc *****************************************




	
