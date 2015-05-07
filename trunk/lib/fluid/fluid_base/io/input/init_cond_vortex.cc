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


/*! \file  init_cond_ABC.cc
 * 
 * @brief Initial conditions as ABC flow.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FFF basis; not for SCFT basis.
 *
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "FluidIO.h"


//*********************************************************************************************

void FluidIO::Init_cond_vortex(FluidVF& U)
{
/*	Real R0 = global.io.double_para(0);
    Real r0 = global.io.double_para(1);

    int rx,ry,rz;
    Real x,y,z;
    Real phi;  // really -minus phi, according to the definition
    
    for (int lx=0; lx<global.field.maxlx; lx++) {
        rx = universal->Get_rx_real_space(lx);
        x = rx*global.field.L[1]/Nx;
        
        for (int ly=0; ly<global.field.maxly; ly++) {
            ry = universal->Get_ry_real_space(ly);
            y = ry*global.field.L[2]/Ny;
            
            for (int lz=0; lz<global.field.maxlz; lz++) {
                rz = universal->Get_rz_real_space(lz);
                z = rz*global.field.L[3]/Nz;
                
                phi = atan((y-R0)/(sqrt((x-M_PI)*(x-M_PI)+(z-M_PI)*(z-M_PI))-r0)) - atan(y/(x-M_PI));
                
                if (universal->Probe_in_me_real_space(rx, ry, rz))
                  universal->Assign_real_field(rx, ry, rz, global.temp_array.Xr, phi);
                
                universal->Get_rz_real_space(lz, rz, 1);
                z = rz*global.field.L[3]/Nz;
                
                phi = atan((y-R0)/(sqrt((x-M_PI)*(x-M_PI)+(z-M_PI)*(z-M_PI))-r0)) - atan(y/(x-M_PI));
                
                if (universal->Probe_in_me_real_space(rx, ry, rz))
                    universal->Assign_real_field(rx, ry, rz, global.temp_array.Xr, phi);
            }
        }
    }
    
    universal->Forward_transform(global.temp_array.Xr, global.temp_array.X);
    
    universal->Xderiv(global.temp_array.X, U.cvf.V1);
    universal->Yderiv(global.temp_array.X, U.cvf.V2);
    universal->Zderiv(global.temp_array.X, U.cvf.V3); */
}

//*********************************************************************************************
// Scalar 
//

void FluidIO::Init_cond_vortex(FluidVF& U, FluidSF& T)
{
}


//*********************************************************************************************
// Vector
//

void FluidIO::Init_cond_vortex(FluidVF& U, FluidVF& W)
{
}

//*********************************************************************************************
//	Vector+scalar
//

void FluidIO::Init_cond_vortex(FluidVF& U, FluidVF& W, FluidSF& T)
{
}


// For GP
void FluidIO::Init_cond_vortex(FluidSF& T)
{
/*	Real R0 = global.io.double_para(0);
    Real r0 = global.io.double_para(1);
    Real A = global.io.double_para(2); */
    
    int rx,ry,rz;
    Real x,y,z;
    Real phi;
    Complex psi_r;
    
    /*if (!(global.field.transpose)) {
        for (int lx=0; lx<local_Nx; lx++) {
            rx = universal->Get_rx_real_space(lx);
            x = rx*global.field.L[1]/Nx;
            
            for (int ly=0; ly<Ny; ly++) {
                ry = universal->Get_ry_real_space(ly);
                y = ry*global.field.L[2]/Ny;
                
                for (int lz=0; lz<Nz; lz++) {
                    universal->Get_rz_real_space(lz, rz, 0);
                    z = rz*global.field.L[3]/Nz;
                    
                    phi = atan((y-R0)/(sqrt((x-M_PI)*(x-M_PI)+(z-M_PI)*(z-M_PI))-r0)); // - atan(y/(x-M_PI));
                    
                    psi_r = Complex(A*cos(phi), A*sin(phi));
                    
                    T.rsf.Fr(lx,ly,lz) = psi_r;
                }
            }
        }
        
        T.Forward_transform();
    }*/
}




//******************************** End of Init_cond_ABC.cc ************************************

