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

/*! \file  force_model.cc
 * 
 * @brief Forcing when low k modes are forced: TG, ABC flows etc.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "FORCE.h"


//*********************************************************************************************

/** @brief Taylor Green flow  Forcing
 * 
 * @note Force_x = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Force_y = -amp*cos(k0 x)sin(k0 y) cos(k0 z)
 *		Force_z = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void FORCE::Setup_Taylor_Green_force_field(FluidVF& U, int k0, Real amp)
{
	if (U.force_switch == true) {
		U.Force1 = 0.0; 
		U.Force2 = 0.0;  
		U.Force3 = 0.0;			
		
		// CHOOSE SCC type
		if (basis_type == "SSS") {
			universal->Assign_spectral_field(k0, k0, k0, U.Force1, amp/8);
			universal->Assign_spectral_field(k0, k0, k0, U.Force2, -amp/8);
		}
		
		// CHOOSE SCF type
		else if (basis_type == "SSF") {
			universal->Assign_spectral_field(k0, k0, k0, U.Force1, amp/8);
			universal->Assign_spectral_field(k0, k0, k0, U.Force2, -amp/8);
		}
		
		// CHOOSE SFF type
		else if (basis_type == "SFF") {
			
			universal->Assign_spectral_field(k0, k0, k0, U.Force1, Complex(amp/8, 0.0));
			universal->Assign_spectral_field(k0, -k0, k0, U.Force1, Complex(amp/8, 0.0));
			
			universal->Assign_spectral_field(k0, k0, k0, U.Force2, Complex(0, amp/8));
			universal->Assign_spectral_field(k0, -k0, k0, U.Force2, Complex(0, -amp/8));
		}
		
		else if (basis_type == "FFF" || basis_type == "FFFW") {
			universal->Assign_spectral_field(k0, k0, k0, U.Force1, Complex(0, -amp/8));
			universal->Assign_spectral_field(k0, -k0, k0, U.Force1, Complex(0, -amp/8));
			universal->Assign_spectral_field(-k0, k0, k0, U.Force1, Complex(0, amp/8));
			universal->Assign_spectral_field(-k0, -k0, k0, U.Force1, Complex(0, amp/8));
			
			universal->Assign_spectral_field(k0, k0, k0, U.Force2, Complex(0, amp/8));
			universal->Assign_spectral_field(-k0, k0, k0, U.Force2, Complex(0, amp/8));
			universal->Assign_spectral_field(k0, -k0, k0, U.Force2, Complex(0, -amp/8));
			universal->Assign_spectral_field(-k0, -k0, k0, U.Force2, Complex(0, -amp/8));
		}
	}
}



//*********************************************************************************************

/** @brief ABC forcing
 * 
 * @note \f$ Force_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ Force_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ Force_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @param k0
 *	@param amp  
 *
 */
void FORCE::Setup_ABC_force_field(FluidVF& U, int k0, Real amp, Real A, Real B, Real C)
{
	if (U.force_switch == true) {
		U.Force1 = 0.0; 
		U.Force2 = 0.0;  
		U.Force3 = 0.0; 
		// initialize

		if (basis_type == "FFF" || basis_type == "FFFW") {
			universal->Assign_spectral_field(0, k0, 0, U.Force1, Complex(amp*B/2, 0));
			universal->Assign_spectral_field(0, -k0, 0, U.Force1, Complex(amp*B/2, 0));
			universal->Assign_spectral_field(0, 0, k0, U.Force1, Complex(0, -amp*C));
			
			universal->Assign_spectral_field(k0, 0, 0, U.Force2, Complex(0, -amp*A/2));
			universal->Assign_spectral_field(-k0, 0, 0, U.Force2, Complex(0, amp*A/2));
			universal->Assign_spectral_field(0, 0, k0, U.Force2, Complex(amp*C, 0));
			
			universal->Assign_spectral_field(k0, 0, 0, U.Force3, Complex(amp*A/2, 0));
			universal->Assign_spectral_field(-k0, 0, 0, U.Force3, Complex(amp*A/2, 0));
			universal->Assign_spectral_field(0, k0, 0, U.Force3, Complex(0, -amp*B));
		}
		
		else {
			cout << "ABC field not allowed in the used basis " << endl;
			exit(1);
		}
	}	
}


//*********************************************************************************************

	// WORK ON IT>>
void FORCE::Setup_SIX_MODE_force_field(FluidVF& U, int k0, Real amp101, Real amp011, Real amp112, Real h)
{
/*
	if (force_switch == true) {
		U.Force1 = 0.0;
		U.Force2 = 0.0;  
		U.Force3 = 0.0; // initialize

		if (basis_type == "FFF" || basis_type == "FFFW")
		{
			int lx_k0 = Get_lx("FFF", k0, N);
			int lx_minus_k0 = Get_lx("FFF", -k0, N);
			int ly_k0 = Get_ly3D("FFF", k0, N);
			int ly_minus_k0 = Get_ly3D("FFF", -k0, N);
			
			Real factor1 = 2.0/sqrt(2+h*h) * amp101;
			Real factor2 = 2.0/sqrt(2+h*h) * amp011;
			Real factor3 = 4.0/sqrt(6+10*h*h) * amp112;
		
			if (my_id == master_id)  
			{
				(U.Force1)(0, ly_k0, k0)  = Complex(-h*(factor2/4), 0.0);
				(U.Force1)(0, ly_minus_k0, k0) = Complex(h*(factor2/4), 0.0);
				
				(U.Force2)(0, ly_k0, k0)  = (I) * (factor2/4);
				(U.Force2)(0, ly_minus_k0, k0) = (-I) * (factor2/4);
				
				(U.Force3)(0, ly_k0, k0)  = (-I) * (factor2/4);
				(U.Force3)(0, ly_minus_k0, k0) = (-I) * (factor2/4);
			}
			
			if ( (lx_k0 >= 0) && (lx_k0 < local_N1) )
			{
				(U.Force1)(lx_k0, 0, k0)  = I * (factor1/4);
				(U.Force2)(lx_k0, 0, k0)  = Complex(-h*(factor1/4), 0.0);
				(U.Force3)(lx_k0, 0, k0)  =  (-I) * (factor1/4);
				
				(U.Force1)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8);
				(U.Force1)(lx_k0, ly_minus_k0, 2*k0) = I * (factor3/8);
				
				(U.Force2)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8) + Complex(h*factor3/4, 0.0);
				(U.Force2)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + Complex(h*factor3/4, 0.0);
				
				(U.Force3)(lx_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) - Complex(h*factor3/8, 0.0);
				(U.Force3)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + Complex(h*factor3/8, 0.0);	
			}
		
			if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
			{
				(U.Force1)(lx_minus_k0, 0, k0)  = -I * (factor1/4);		
				(U.Force2)(lx_minus_k0, 0, k0)  = Complex(h*(factor1/4), 0.0);
				(U.Force3)(lx_minus_k0, 0, k0)  =  (-I) * (factor1/4);
				
				(U.Force1)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8);
				(U.Force1)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8);
			
				(U.Force2)(lx_minus_k0, ly_k0, 2*k0)  = I * (factor3/8) - Complex(h*factor3/4, 0.0);
				(U.Force2)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
																- Complex(h*factor3/4, 0.0);
			
				(U.Force3)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) + Complex(h*factor3/8, 0.0);
				(U.Force3)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
																- Complex(h*factor3/8, 0.0);
			}	
			
		}
		
		else if ((basis_type == "SCFT") && (my_id == master_id))
		{
			cout << "Setup_SIX_MODE_field Not allowed in SCFT basis_type " << endl;
			exit(0);
		}
	}
 */
}

//*********************************************************************************************
/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once force field \f$ \vec{F}(kx, ky, kz) \f$ is given.
 *  
 * @bug Does not work across the processors.
 *
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FFF: \f$ \vec{F}(-kx, -ky, kz) = \vec{F}^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ \vec{F}(kx, -ky, kz) = \vec{F}^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(F_y(kx, 0, kz)) = imag(F_z(kx, 0, kz)) = 0\f$,
 *					and \f$ F_x(kx, 0, kz) = 0 \f$.
 *				
 */
void FORCE::Add_complex_conj_force(FluidVF& U, int kx, int ky, int kz, Complex Fx, Complex Fy, Complex Fz)	
{	
	TinyVector<Complex,3> localF;
	
	// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FFF" || basis_type == "FFFW") && ((kz == 0) || (kz == N[3]/2)) ) {
		localF = conj(Fx), conj(Fy), conj(Fz);
		universal->Assign_spectral_field(-kx, -ky, kz, U.Force1, U.Force2, U.Force3, localF);
		if (master)
			cout << "Complex-conj(Force) added for k = (" << -kx << "," << -ky << "," << kz << ")" << endl;	
	}
	
	else if ((basis_type == "SFF") && ((kz == 0) || (kz == N[3]/2)) ) {          
		localF = conj(Fx), conj(Fy), conj(Fz);
		universal->Assign_spectral_field(kx, -ky, kz, U.Force1, U.Force2, U.Force3, localF);
		if (master)
			cout << "Complex-conj(Force) added for k = (" << kx << "," << -ky << "," << kz << ")" << endl;	
	}	
}	

	//
	// Scalar 
	//

void FORCE::Add_complex_conj_force(FluidSF& T, int kx, int ky, int kz, Complex localG)	
{	
	
	// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FFF" || basis_type == "FFFW") && ((kz == 0) || (kz == N[3]/2)) )	{          
		universal->Assign_spectral_field(-kx, -ky, kz, T.Force, localG);
		if (master)
			cout << "Complex-conj(Force) added for k = (" << -kx << "," << -ky << "," << kz << ")" << endl;	
	}
	
	else if ((basis_type == "SFF")  && ((kz == 0) || (kz == N[3]/2)) ) 	{         
		universal->Assign_spectral_field(kx, -ky, kz, T.Force, localG);
		if (master)
			cout << "Complex-conj(Force) added for k = (" << kx << "," << -ky << "," << kz << ")" << endl;	
	}
}


//*********************************************************************************************

/** @brief D=3:  force field \f$ \vec{F}(kx, ky, kz) = (F_x, F_y, F_z) \f$ 
 *				and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, kz, Complex Fx, Complex Fy, Complex Fz)
 */
void FORCE::Assign_force_and_comp_conj(FluidVF& U, int kx, int ky, int kz, Complex Fx, Complex Fy, Complex Fz)
{
	TinyVector<Complex,3> localF(Fx, Fy, Fz);
	
	universal->Assign_spectral_field(kx, ky, kz, U.Force1, U.Force2, U.Force3, localF);	// in appropriate proc.
	
	Add_complex_conj_force(U, kx, ky, kz, Fx, Fy, Fz);	
}

	//
	//
	//*********************************************************************************************

/** @brief D=2:  force field for scalar \f$ F(kx, ky, kz) = FG \f$ and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, int kz, Complex ForceG)
 */
void  FORCE::Assign_force_and_comp_conj(FluidSF& T, int kx, int ky, int kz, Complex localG)
{
	universal->Assign_spectral_field(kx, ky, kz, T.Force, localG);
	
	Add_complex_conj_force(T, kx, ky, kz, localG);		
}


	//*********************************************************************************************

/** @brief D=3:  force field \f$ \vec{F}(kx, ky, kz) = (F_x, F_y, F_z) \f$ 
 *				and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, kz, Complex Fx, Complex Fy, Complex Fz)
 */
void FORCE::Assign_force(FluidVF& U, int kx, int ky, int kz, Real Fx, Real Fy, Real Fz)
{
	TinyVector<Real,3> localF(Fx, Fy, Fz);
	
	universal->Assign_spectral_field(kx, ky, kz, U.Force1, U.Force2, U.Force3, localF);	// in appropriate proc.
}

	//
	//
	//*********************************************************************************************

/** @brief D=2:  force field for scalar \f$ F(kx, ky, kz) = FG \f$ and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, int kz, Complex ForceG)
 */
void  FORCE::Assign_force(FluidSF& T, int kx, int ky, int kz, Real localG)
{
	universal->Assign_spectral_field(kx, ky, kz, T.Force, localG);
}

//*************************************************************************************

void FORCE::Model_force_spectrum(Real force_spectrum_amplitude, Real force_spectrum_exponent, Array<Real,1> Sk)
{
	Real k;
	
	for (int i=0; i < (Sk.length())(0); i++) {
		k = 1.0*i;
		Sk(i) = force_spectrum_amplitude * pow(k,force_spectrum_exponent)/(0.01);
    }
    
  //  cout << "Sk = " << Sk << endl;
}



//*******************************  End of force_model.cc **************************************

	
