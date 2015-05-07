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

/*! \file energy_tr_B0.cc
 * 
 * @brief  MHD: Computes \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$ over 
 *			a shell or ring. Effect of B0 to energy transfer..
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "EnergyTr.h"

//*********************************************************************************************


/// Compute for ET shells \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void EnergyTr::Compute_shell_ET_B0(FluidVF& U, FluidVF& W)
{
	
	TinyVector<Real,3> B0;
	
	if (my_id == master_id)
		B0 = real((W.cvf.V1)(0,0,0)), real((W.cvf.V2)(0,0,0)), real((W.cvf.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<Real*>(B0.data()), data_size, MPI_Real, master_id, MPI_COMM_WORLD);
	
	universal->Shell_mult_all_imagVW_B0(B0, U.cvf.V1, U.cvf.V2, U.cvf.V3, W.cvf.V1, W.cvf.V2, W.cvf.V3, energy_tr_shell_B0);	
												
}


//*********************************************************************************************

/// Compute for ET rings \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void EnergyTr::Compute_ring_ET_B0(FluidVF& U, FluidVF& W)
{
	TinyVector<Real,3> B0;
	
	if (my_id == master_id)
		B0 = real((W.cvf.V1)(0,0,0)), real((W.cvf.V2)(0,0,0)), real((W.cvf.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<Real*>(B0.data()), data_size, MPI_Real, master_id, MPI_COMM_WORLD);
	
	universal->Ring_mult_all_imagVW_B0(B0, U.cvf.V1, U.cvf.V2, U.cvf.V3, W.cvf.V1, W.cvf.V2, W.cvf.V3, energy_tr_ring_B0);	
												
}



//*********************************************************************************************

/// Compute for ET cylindrical rings 
///		\f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void EnergyTr::Compute_cylindrical_ring_ET_B0(FluidVF& U, FluidVF& W)
{

	TinyVector<Real,3> B0;
	
	if (my_id == master_id)
		B0 = real((W.cvf.V1)(0,0,0)), real((W.cvf.V2)(0,0,0)), real((W.cvf.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<Real*>(B0.data()), data_size, MPI_Real, master_id, MPI_COMM_WORLD);
	
	universal->Cyl_ring_mult_all_imagVW_B0(B0, U.cvf.V1, U.cvf.V2, U.cvf.V3, W.cvf.V1, W.cvf.V2, W.cvf.V3, energy_tr_cylindrical_ring_B0);
}

//******************************** End of energy_tr_B0.cc  ************************************



