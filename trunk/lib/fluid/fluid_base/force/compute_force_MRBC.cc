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
/*! \file  compute_force_RB.cc
 * 
 * @brief Force for RB convection
 *
 *		For Pr_large, u_small:		 F1 = Ra*Pr*T(k);  
 *		For Pr_large, u_large:		 F1 = Ra*T(k);    	
 *		For Pr_small or 0, u_small:  F1 = Ra*T(k);	   
 *		For Pr_small or 0, u_large:  F1 = Pr*T(k)
 *
 *		For Pr_large				Ftheta = V1;    	
 *		For Pr_small				Ftheta = V1/Pr;	   
 *		For Pr = 0					Ftheta = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 
#include "FORCE.h"


//*********************************************************************************************

void FORCE::Compute_force_MRBC(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	/*int rx, ry, rz;
	Real x, Dxyz, Mxyz, Bxyz;
	bool is_real_part= true;	
	
	T1.Inverse_transform();
	T2.Inverse_transform();
	
	for (int lx=0; lx<global.field.maxlx; lx++) {
		rx = universal->Get_rx_real_space(lx);
		x = rx*global.field.L[1]/Nx;
		
		for (int ly=0; ly<global.field.maxly; ly++)
			for (int lz=0; lz<global.field.maxlz; lz++) {
				ry = universal->Get_ry_real_space(ly);
				universal->Get_rz_real_space(lz, rz, is_real_part);
					// even rz
				Dxyz = universal->Get_real_field(rx, ry, rz, T1.rsf.Fr);	// D=T1;
				Mxyz = universal->Get_real_field(rx, ry, rz, T2.rsf.Fr);	// M =T2
				
				Bxyz = max(Mxyz, Dxyz+global.MRBC.SSD+ x*(1-global.MRBC.CSA-global.MRBC.RaD/global.MRBC.RaM));						   
				universal->Assign_real_field(rx, ry, rz, global.temp_array.Xr, Bxyz);
				
					// odd rz's
				Dxyz = universal->Get_real_field(rx, ry, rz+1, T1.rsf.Fr);	// D=T1;
				Mxyz = universal->Get_real_field(rx, ry, rz+1, T2.rsf.Fr);	// M =T2
				
				Bxyz = max(Mxyz, Dxyz+global.MRBC.SSD+ x*(1-global.MRBC.CSA-global.MRBC.RaD/global.MRBC.RaM));						   
				universal->Assign_real_field(rx, ry, rz+1, global.temp_array.Xr, Bxyz);
			}
	}
	
    universal->Forward_transform(global.temp_array.Xr, U.Force1); 
	

	U.Force2 = 0.0; 
	U.Force3 = 0.0;
	
	// For the temperature field

	T1.Force =  complex<Real>(global.MRBC.RaD/global.MRBC.RaM, 0)*(U.cvf.V1); 
	T2.Force =  (U.cvf.V1); 
	 */
	
}


//************************ End of compute_force_RB.cc *****************************************

