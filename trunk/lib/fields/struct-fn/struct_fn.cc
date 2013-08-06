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
/*! \file struct_fn.cc
 * 
 * @sa struct_fn.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 
 
#include "struct_fn.h"


//*********************************************************************************************


// Along the periodic dirn, max dist is N[i]/2.
void Compute_dr
(
	string basis_type,
	int N[], 
	TinyVector<int,3> j1, TinyVector<int,3> j2, 
	TinyVector<DP,3> &dr,
	DP delta_x[]
) 
{

	
	// FOR PERIODIC BOX.. except for SCFT (top,bottom plates) 
	dr(0) = j2(0) - j1(0);
	dr(1) = j2(1) - j1(1);
	dr(2) = j2(2) - j1(2);
	
	if (basis_type == "FOUR")					// Periodic dirn
		if (abs(dr(0)) > N[1]/2) 
			dr(0) = (dr(0)>0) ? (dr(0) - N[1]) : (dr(0) + N[1]);
	
	dr(0) = dr(0)* delta_x[1];  // Straightaway for SCFT
	
	// Y, Z direction- periodic for both FOUR & SCFT
	if (abs(dr(1)) > N[2]/2)
		dr(1) = (dr(1)>0) ? (dr(1) - N[2]) : (dr(1) + N[2]);
	
	if (abs(dr(2)) > N[3]/2)
		dr(2) = (dr(2)>0) ? (dr(2) - N[3]) : (dr(2) + N[3]);
	
	dr(1) = dr(1) * delta_x[2];
	dr(2) = dr(2) * delta_x[3]; 
	 
	
	
	/*
	// For Box geometry with walls.. no periodicity
	dr(0) = (j2(0) - j1(0))*delta_x[1];
	dr(1) = (j2(1) - j1(1))*delta_x[2];
	dr(2) = (j2(2) - j1(2))*delta_x[3];
	 */
	 
}



//**************************************************************************************************
		

void RVF::Compute_dSt
(
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	TinyVector<DP,3> V1, TinyVector<DP,3> V2
)
{
	TinyVector<DP,3>	dr, dr_hat;
	TinyVector<DP,3>	dV;	
	DP					abs_dr, dV_pll, dV_perp;
	
	
	Compute_dr(RV_basis_type, Nrv, j1, j2, dr, RV_Delta_x);
	
	abs_dr = sqrt(dot(dr,dr));
	dr_hat = dr(0)/abs_dr, dr(1)/abs_dr, dr(2)/abs_dr;
	
	dV = V2(0)-V1(0), V2(1)-V1(1), V2(2)-V1(2);
	
	dV_pll = dot(dV, dr_hat);
	
	if (abs(dot(dV,dV)) > pow2(dV_pll))
		dV_perp = sqrt(dot(dV,dV) - pow2(dV_pll));
	else
		dV_perp = 0.0;
		
	int  rindex = (int) ceil(abs_dr*RV_structurefn_r_arraysize/RV_structurefn_r_max);
	
	// Skip when P and Q overlap, 
	// or when the distance is larger than max_radius_inside_real_space.
	if ((rindex > 0) && (rindex <= RV_structurefn_r_arraysize)) {
		(*RV_St_count)(rindex) = (*RV_St_count)(rindex)+1.0;

		for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++) {
				(*RV_St)(rindex,q-RV_structurefn_q_min,0) += my_pow(dV_pll, q);
				(*RV_St)(rindex,q-RV_structurefn_q_min,1) += my_pow(dV_perp, q);
		}
	}
}


//*************************************************************RV_structurefn_neighbour_index*************************************

void RVF::Compute_structure_function_unit
(
	int otherproc, 
	Array<complx,3> otherV1,
	Array<complx,3> otherV2,
	Array<complx,3> otherV3
)
{	
	
	TinyVector<int,3>	j1, j2;
	TinyVector<DP,3>	v_1, v_2;
	
	for (int Pi1=0; Pi1<=(local_N1-1); Pi1=Pi1+RV_structurefn_neighbour_index)
		for (int Pi2=0; Pi2<=(Nrv[2]-1); Pi2=Pi2+RV_structurefn_neighbour_index) 
			for (int Pi3=0; Pi3<=(Nrv[3]/2-1); Pi3=Pi3+RV_structurefn_neighbour_index) {
				
				j1(0) = local_N1_start + Pi1;
				j1(1) = Pi2;
				
				RV_structure_fn_Q_ends(Pi1, Pi2, Pi3);
				// Gives Qi_start and Qi_ends
				
				for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RV_structurefn_neighbour_index)
					for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RV_structurefn_neighbour_index)
						for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RV_structurefn_neighbour_index) 
						{
							 
							j2(0) = otherproc*local_N1 + Qi1;
							j2(1) = Qi2;
														
							// first set
							v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
								  real((*V3r)(Pi1,Pi2,Pi3)); 
							
							v_2 = real((otherV1)(Qi1,Qi2,Qi3)), real((otherV2)(Qi1,Qi2,Qi3)), 
								  real((otherV3)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3;
							j2(2) = 2*Qi3;
							Compute_dSt(j1, j2, v_1, v_2);
							
							// Second set 		
							v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
								  real((*V3r)(Pi1,Pi2,Pi3)); 
							
							v_2 = imag((otherV1)(Qi1,Qi2,Qi3)), imag((otherV2)(Qi1,Qi2,Qi3)), 
								  imag((otherV3)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3;
							j2(2) = 2*Qi3+1;
							Compute_dSt(j1, j2, v_1, v_2);
							
							// Third set 	
							v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
								  imag((*V3r)(Pi1,Pi2,Pi3)); 
							
							v_2 = real((otherV1)(Qi1,Qi2,Qi3)), real((otherV2)(Qi1,Qi2,Qi3)), 
								  real((otherV3)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3+1;
							j2(2) = 2*Qi3;
							Compute_dSt(j1, j2, v_1, v_2);
							
							// Fourth set 		
							v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
								  imag((*V3r)(Pi1,Pi2,Pi3)); 
							
							v_2 = imag((otherV1)(Qi1,Qi2,Qi3)), imag((otherV2)(Qi1,Qi2,Qi3)), 
								  imag((otherV3)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3+1;
							j2(2) = 2*Qi3+1;
							Compute_dSt(j1, j2, v_1, v_2);  
						}		// end of  loops
			}
	
}



void RVF::Compute_local_structure_function()
{	
	int data_size = 2* local_N1 * Nrv[2] * (Nrv[3]/2 + 1);
	int tag = 123;
	int receive_from, send_to;
	MPI_Request request, request2;
	
	// allocate otherprocVi vars
	Array<complx,3> otherprocV1(local_N1,Nrv[2],(Nrv[3]/2)+1); 
	Array<complx,3> otherprocV2(local_N1,Nrv[2],(Nrv[3]/2)+1); 
	Array<complx,3> otherprocV3(local_N1,Nrv[2],(Nrv[3]/2)+1); 
	
	*RV_St = 0.0;
	*RV_St_count = 0;
	
	// with itself
	otherprocV1 = *V1r;
	otherprocV2 = *V2r;
	otherprocV3 = *V3r;
	Compute_structure_function_unit(my_id, otherprocV1, otherprocV2, otherprocV3);
	
	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i >= 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id + i) % numprocs;
		
		// first component
		MPI_Isend( reinterpret_cast<DP*>((*V1r).data()), data_size, 
				MPI_DP, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<DP*>((otherprocV1).data()), data_size, 
				MPI_DP, receive_from, tag, MPI_COMM_WORLD, &request2);

		MPI_Wait(&request2, &status);
		
		// Second component
		if (Nrv[2] > 1)   // if 2D, then ignore V2 transfer.
		{
			MPI_Isend( reinterpret_cast<DP*>((*V2r).data()), data_size, 
				 MPI_DP, send_to, tag, MPI_COMM_WORLD, &request);
	
			MPI_Irecv( reinterpret_cast<DP*>((otherprocV2).data()), data_size, 
				 MPI_DP, receive_from, tag, MPI_COMM_WORLD, &request2);

			MPI_Wait(&request2, &status);
		}
		
		
		// Third component
		MPI_Isend( reinterpret_cast<DP*>((*V3r).data()), data_size, 
				 MPI_DP, send_to, tag, MPI_COMM_WORLD, &request);
	
		MPI_Irecv( reinterpret_cast<DP*>((otherprocV3).data()), data_size, 
				 MPI_DP, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status); 
		
		Compute_structure_function_unit(receive_from, 
							otherprocV1, otherprocV2, otherprocV3);
	}

}


void RVF::RV_Compute_structure_function()
{	

	Compute_local_structure_function();
	
	int data_size = (*RV_St).size(); 
	
	MPI_Reduce(reinterpret_cast<DP*>((*RV_St).data()), 
			   reinterpret_cast<DP*>((*RV_St_final).data()), 
			   data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 

	data_size = (*RV_St_count).size();  
	
	MPI_Reduce(reinterpret_cast<int*>((*RV_St_count).data()), 
			   reinterpret_cast<int*>((*RV_St_count_final).data()), 
			   data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	

	// Normalize St(..)	
	
	if (my_id == master_id) {
		for (int rindex=1; rindex<=RV_structurefn_r_arraysize; rindex++)
			if ((*RV_St_count_final)(rindex) > MYEPS2)
				for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++) {
		
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,0) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,0)
							/(*RV_St_count_final)(rindex);
					
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,1) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,1)
							/(*RV_St_count_final)(rindex);
				}
	}
		 
}



//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************

void RSF::Compute_dSt
(
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	DP F1, DP F2
)
{	
	TinyVector<DP,3>	dr;	
	DP					abs_dr;
	
	Compute_dr(RS_basis_type, Nrs, j1, j2, dr, RS_Delta_x);
	abs_dr = sqrt(dot(dr,dr));
	
	int rindex = (int) ceil(abs_dr*RS_structurefn_r_arraysize/RS_structurefn_r_max);
	
	
	// Skip when P and Q overlap, 
	// or when the distance is larger than max_radius_inside_real_space.
	if ((rindex >=0) && (rindex <= RS_structurefn_r_arraysize)) {
		(*RS_St_count)(rindex) = (*RS_St_count)(rindex)+1.0;
		
		for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++)
			(*RS_St)(rindex,q-RS_structurefn_q_min) += my_pow((F2-F1), q);
	}
}



//*********************************************************************************************

void RSF::Compute_structure_function_unit
(
	int otherproc,
	Array<complx,3> otherprocF
)
{
	TinyVector<int,3>	j1, j2;
	DP					F1, F2;
	
	for (int Pi1=0; Pi1<=(local_N1-1); Pi1=Pi1+RS_structurefn_neighbour_index)
		for (int Pi2=0; Pi2<=(Nrs[2]-1); Pi2=Pi2+RS_structurefn_neighbour_index) 
			for (int Pi3=0; Pi3<=(Nrs[3]/2-1); Pi3=Pi3+RS_structurefn_neighbour_index) {
				
				j1(0) = local_N1_start + Pi1;
				j1(1) = Pi2; 
				
				RS_structure_fn_Q_ends(Pi1, Pi2, Pi3);
				
				for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RS_structurefn_neighbour_index)
					for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RS_structurefn_neighbour_index)
						for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RS_structurefn_neighbour_index)
						{
							
							j2(0) = otherproc*local_N1 + Qi1;
							j2(1) = Qi2;
							
							// first set
							F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
							F2 = real((otherprocF)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3;
							j2(2) = 2*Qi3;
							Compute_dSt(j1, j2, F1, F2);
							
							// Second set 		
							F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
							F2 = imag((otherprocF)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3;
							j2(2) = 2*Qi3+1;
							Compute_dSt(j1, j2, F1, F2);
							
							// Third set 		
							F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
							F2 = real((otherprocF)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3+1;
							j2(2) = 2*Qi3;
							Compute_dSt(j1, j2, F1, F2);
							
							// Fourth set 		
							F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
							F2 = imag((otherprocF)(Qi1,Qi2,Qi3));
							
							j1(2) = 2*Pi3+1;
							j2(2) = 2*Qi3+1;
							Compute_dSt(j1, j2, F1, F2);
						}	
			}	// end of  loops
}
	


void RSF::Compute_local_structure_function()
{
	
	int data_size = 2* local_N1 * Nrs[2] * (Nrs[3]/2 + 1);
	int tag = 123;
	int receive_from, send_to;	
	MPI_Request request, request2;
	
	// allocate otherprocVi vars
	Array<complx,3> otherprocF(local_N1,Nrs[2],(Nrs[3]/2)+1);
	
	*RS_St = 0.0;
	*RS_St_count = 0;
	
	// with itself
	otherprocF = *Fr;
	Compute_structure_function_unit(my_id, otherprocF);
	
	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i > 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id+i) % numprocs;
		
		MPI_Isend( reinterpret_cast<DP*>((*Fr).data()), data_size, 
				 MPI_DP, send_to, tag, MPI_COMM_WORLD, &request);
	
		MPI_Irecv( reinterpret_cast<DP*>((otherprocF).data()), data_size, 
				 MPI_DP, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		Compute_structure_function_unit(receive_from, otherprocF);
	}
}


void RSF::RS_Compute_structure_function()
{
	
	Compute_local_structure_function();
	
	int data_size = (*RS_St).size(); 
	
	MPI_Reduce(reinterpret_cast<DP*>((*RS_St).data()), 
			   reinterpret_cast<DP*>((*RS_St_final).data()), 
			   data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	data_size = (*RS_St_count).size();

	MPI_Reduce(reinterpret_cast<int*>((*RS_St_count).data()), 
			   reinterpret_cast<int*>((*RS_St_count_final).data()), 
			   data_size, MPI_DP, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	// Normalize St(..)	
	if (my_id == master_id)
		for (int rindex=1; rindex<=RS_structurefn_r_arraysize; rindex++) {
			if ((*RS_St_count_final)(rindex) > MYEPS2)
				for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++)
					(*RS_St_final)(rindex,q-RS_structurefn_q_min) 
					= (*RS_St_final)(rindex,q-RS_structurefn_q_min)/(*RS_St_count_final)(rindex);
		}
	
}




//*********************************   End of struct_fn.cc **********************


