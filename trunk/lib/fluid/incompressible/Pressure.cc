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

/*! \file compute_pressure.h 
 * 
 * @brief Functions for computing pressure and pressure spectrum.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date 30 August 2008
 * 
 * @bug  No known bugs
 */ 

//*********************************************************************************************

#include "Pressure.h"
#include "extern_var_incompress.h"

//********************************************************************************************* 

/** @brief Compute pressure and puts it in the CSF of NLIN.
 *
 *	Before entering this function:  
 *		nlin \f$ \vec{N} \f$ contains the nonlinear of udot equation (in Fourier space).
 *		Fluid & Passive scalar: nlin = FT[U.grad U] 
 *		MHD: nlin = T[U.grad U - B.grad B]
 *
 *	Procedure:  \f$ p(K) = \mathcal{F}(\nabla \cdot \vec{N}) \f$, then \f$ p(K) = p(K)/K^2 \f$.
 *
 *  @return  pressure that is stored in F of CSF of NLIN class.   
 */
void Pressure::Compute_pressure(FluidVF &U)
{
    if (global.program.kind != "KEPLERIAN") {
        U.Compute_divergence_nlin(F);  // div(nlin) -> F  (the proc uses temp arrays)

        Divide_ksqr();			// F = F/k^2
        
        if (!global.io.output_pressure_spectrum_done) {
            fluidIO_incompress.Output_pressure_spectrum(*this);	
            global.io.output_pressure_spectrum_done = true;
        }
        
        if (!global.io.output_pressure_done) {
            fluidIO_incompress.Output_pressure (*this);	
            global.io.output_pressure_done = true;
        }
    }
    
    // For Keplerian. ONLY FFF SLAB BASIS
    // -iq omega Ky ux + i N.mu = mu'^2 p(t)
    else {
        Real omega_keplerian = global.force.double_para(0);
        Real q_keplerian = global.force.double_para(1);
        
        Real q_omega_t = q_keplerian*omega_keplerian*global.time.keplerian;
        
        U.Compute_divergence_nlin(F);
        // div(nlin) = ik.N + iq omega t  Ky Nx -> F  (the proc uses temp arrays)
        
        Real Kx,Ky,Kz,Ksqr, mu_sqr;
        
        for (int lx=0; lx<local_Nx; lx++)
            for (int ly=0; ly<Ny; ly++)
                for (int lz=0; lz<=Nz/2; lz++) {
                    Kx = (universal->Get_kx(lx))*kfactor[1];
                    Ky = (universal->Get_ky(ly))*kfactor[2];
                    Kz = (universal->Get_kz(lz))*kfactor[3];
                    Ksqr = pow2(Kx)+pow2(Ky)+pow2(Kz);
                    mu_sqr = Ksqr+pow2(q_omega_t*Ky)+2*q_omega_t*Kx*Ky;
                    
                    F(lx,ly,lz) /= mu_sqr;
                }
        
        if (master)
            F(0,0,0) = 0.0;
    }
}



//*******************************  End of compute_pressure.cc  ********************************

