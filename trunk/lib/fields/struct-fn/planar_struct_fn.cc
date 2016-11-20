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
/*! \file planar_struct_fn.cc
 * 
 * @sa planar_struct_fn.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 
 
// #include "planar_struct_fn.h"
#include "struct_fn.h"


/**********************************************************************************************

	Planar structure function S(z, r, q) if anisotropy along z.
	Similar scheme for other anisotropy dirn.

***********************************************************************************************/

// Computes Structure function for y-z plane indexed by pll_index, and 
// print in file_out. No need of any data transfer across processors.
void RVF::Compute_planar_structure_function_unit_anis1(ofstream& file_out, int pll_index)
{
	TinyVector<int,3>	j1, j2;
	TinyVector<Real,3>	v_1, v_2;
	
	int data_size_St, data_size_St_count;
	int tag_St = 111;
	int tag_St_count = 222;
	
	int Pi1 = pll_index;
	int Qi1 = pll_index;
	
	j1(0) = local_N1_start + Pi1;
	j2(0) = local_N1_start + Pi1;
	
	for (int Pi2=0; Pi2<Nrv[2]; Pi2=Pi2+RV_structurefn_neighbour_index) 
		for (int Pi3=0; Pi3<Nrv[3]/2; Pi3=Pi3+RV_structurefn_neighbour_index) {
			
			j1(1) = Pi2; 
			
			RV_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			// Gives Qi_start and Qi_ends
	
			for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RV_structurefn_neighbour_index)
				for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RV_structurefn_neighbour_index) {
					
					j2(1) = Qi2;
					
					// first set
					v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
							real((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = real((*V1r)(Qi1,Qi2,Qi3)), real((*V2r)(Qi1,Qi2,Qi3)), 
							real((*V3r)(Qi1,Qi2,Qi3));
				
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Second set 		
					v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
							real((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = imag((*V1r)(Qi1,Qi2,Qi3)), imag((*V2r)(Qi1,Qi2,Qi3)), 
							imag((*V3r)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Third set 	
					v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
							imag((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = real((*V1r)(Qi1,Qi2,Qi3)), real((*V2r)(Qi1,Qi2,Qi3)), 
							real((*V3r)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Fourth set 		
					v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
							imag((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = imag((*V1r)(Qi1,Qi2,Qi3)), imag((*V2r)(Qi1,Qi2,Qi3)), 
							imag((*V3r)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, v_1, v_2);  
				}		// end of  loops
	}
	
	
	// if result in master_id, print it.
	// else, send to the master_id.  Master_id will print it
	if (my_id == master_id) {
		file_out << "%% pll_index = "  << local_N1_start + pll_index << endl;
		
		for (int rindex=1; rindex<=RV_structurefn_r_arraysize; rindex++) {
			
			file_out	<< rindex << "  ";
			
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++) 
				if ((*RV_St_count)(rindex) > MYEPS2) {
					
					(*RV_St)(rindex,q-RV_structurefn_q_min,0) 
					= (*RV_St)(rindex,q-RV_structurefn_q_min,0)/(*RV_St_count)(rindex);
					
					(*RV_St)(rindex,q-RV_structurefn_q_min,1) 
					= (*RV_St)(rindex,q-RV_structurefn_q_min,1)/(*RV_St_count)(rindex);
					
					file_out	<< (*RV_St)(rindex,q-RV_structurefn_q_min,0)	<< " " 
								<< (*RV_St)(rindex,q-RV_structurefn_q_min,1)	<< " ";
				} 
				
			file_out << endl;	
		}
			
		
		for (int source = 1; source <= numprocs-1; source++) {
			
			data_size_St = (*RV_St).size();
			
			MPI_Recv( reinterpret_cast<Real*>((*RV_St).data()), data_size_St, 
					 MPI_Real, source, tag_St, MPI_COMM_WORLD, &status );
					 
			data_size_St_count = (*RV_St_count).size(); 		 
					 
			MPI_Recv( reinterpret_cast<Real*>((*RV_St_count).data()), data_size_St_count, 
					 MPI_Real, source, tag_St_count, MPI_COMM_WORLD, &status );	
					 
			file_out << "%% pll_index = " << local_N1*source + pll_index << endl;
			
			for (int rindex=1; rindex<=RV_structurefn_r_arraysize; rindex++) {
				
				file_out	<< rindex << "  ";
				
				for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
					if ((*RV_St_count)(rindex) > MYEPS2) {
						
						(*RV_St)(rindex,q-RV_structurefn_q_min,0) 
						= (*RV_St)(rindex,q-RV_structurefn_q_min,0)/(*RV_St_count)(rindex);
						
						(*RV_St)(rindex,q-RV_structurefn_q_min,1) 
						= (*RV_St)(rindex,q-RV_structurefn_q_min,1)/(*RV_St_count)(rindex);
						
						file_out	<< (*RV_St)(rindex,q-RV_structurefn_q_min,0)	<< " " 
									<< (*RV_St)(rindex,q-RV_structurefn_q_min,1)	<< " ";
					} 
				
				file_out << endl;	
			}								
			
			
		}			
	}
	
	else {
		data_size_St = (*RV_St).size();
		
		MPI_Send( reinterpret_cast<Real*>((*RV_St).data()), data_size_St, 
				 MPI_Real, master_id, tag_St, MPI_COMM_WORLD );
			
		data_size_St_count = (*RV_St_count).size(); 
				 	 
		MPI_Send( reinterpret_cast<Real*>((*RV_St_count).data()), data_size_St_count, 
				 MPI_Real, master_id, tag_St_count, MPI_COMM_WORLD );		 
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}


void RVF::Compute_planar_structure_function_unit_anis2
(
	int pll_index,
	int otherproc, 
	Array<Complex,2> otherprocV1,
	Array<Complex,2> otherprocV2,
	Array<Complex,2> otherprocV3
)
{
	
	TinyVector<int,3>	j1, j2;
	TinyVector<Real,3>	v_1, v_2;
	
	int Pi2 = pll_index;
	int Qi2 = pll_index;
	
	j1(1) = pll_index;
	j2(1) = pll_index;
		
	for (int Pi1=0; Pi1<local_N1; Pi1=Pi1+RV_structurefn_neighbour_index) 
		for (int Pi3=0; Pi3<Nrv[3]/2; Pi3=Pi3+RV_structurefn_neighbour_index) {
			
			j1(0) = local_N1_start + Pi1;
			
			RV_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			// Gives Qi_start and Qi_ends
			
			for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RV_structurefn_neighbour_index)
				for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RV_structurefn_neighbour_index) {
					
					j2(0) = otherproc*local_N1 + Qi1; 
	
					// first set
					v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
							real((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = real((otherprocV1)(Qi1,Qi3)), real((otherprocV2)(Qi1,Qi3)), 
							real((otherprocV3)(Qi1,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Second set 		
					v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
							real((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = imag((otherprocV1)(Qi1,Qi3)), imag((otherprocV2)(Qi1,Qi3)), 
							imag((otherprocV3)(Qi1,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Third set 	
					v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
							imag((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = real((otherprocV1)(Qi1,Qi3)), real((otherprocV2)(Qi1,Qi3)), 
							real((otherprocV3)(Qi1,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, v_1, v_2);
					
					// Fourth set 		
					v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
							imag((*V3r)(Pi1,Pi2,Pi3)); 
					
					v_2 = imag((otherprocV1)(Qi1,Qi3)), imag((otherprocV2)(Qi1,Qi3)), 
							imag((otherprocV3)(Qi1,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, v_1, v_2);  
				}		// end of  loops
		}		
}


void RVF::Compute_planar_structure_function_unit_anis3
(
	int pll_index,
	int otherproc, 
	Array<Real,2> otherprocV1,
	Array<Real,2> otherprocV2,
	Array<Real,2> otherprocV3,
	int imag_switch
)
{
	TinyVector<int,3>	j1, j2;
	TinyVector<Real,3>	v_1, v_2;
	
	int Pi3 = pll_index;
	int Qi3 = pll_index;
	
	for (int Pi1=0; Pi1<local_N1; Pi1=Pi1+RV_structurefn_neighbour_index) 
		for (int Pi2=0; Pi2<Nrv[2]; Pi2=Pi2+RV_structurefn_neighbour_index) {
			
			j1(0) = local_N1_start + Pi1;
			j1(1) = Pi2;
			
			RV_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			// Gives Qi_start and Qi_ends
			
			for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RV_structurefn_neighbour_index)
				for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RV_structurefn_neighbour_index) {
					
					j2(0) = otherproc*local_N1 + Qi1;
					j2(1) = Qi2;
					
					if (imag_switch == 0) {
						v_1 = real((*V1r)(Pi1,Pi2,Pi3)), real((*V2r)(Pi1,Pi2,Pi3)), 
								real((*V3r)(Pi1,Pi2,Pi3)); 
						
						v_2 = otherprocV1(Qi1,Qi2), otherprocV2(Qi1,Qi2), 
								otherprocV3(Qi1,Qi2);
						
						j1(2) = 2*Pi3;
						j2(2) = 2*Qi3;
						Compute_dSt(j1, j2, v_1, v_2);
					}
					
					else if (imag_switch == 1) {
						v_1 = imag((*V1r)(Pi1,Pi2,Pi3)), imag((*V2r)(Pi1,Pi2,Pi3)), 
								imag((*V3r)(Pi1,Pi2,Pi3)); 
						
						v_2 = otherprocV1(Qi1,Qi2), otherprocV2(Qi1,Qi2), 
								otherprocV3(Qi1,Qi2);
						
						j1(2) = 2*Pi3+1;
						j2(2) = 2*Qi3+1;
						Compute_dSt(j1, j2, v_1, v_2);  
					}
				} // end of Q loop	
		} // end of P loop	
}


//**************************************************************************************************

// Computes Structure function for x-z plane indexed by pll_index.  
// Each processor computes for a subpart of the plane; 
// these subparts are summed in the master node.
void RVF::Compute_local_planar_structure_function_anis2(int pll_index)
{

	int data_size = 2* local_N1 * (Nrv[3]/2 + 1);
	int tag = 123;
	int receive_from, send_to;
	MPI_Request request, request2;
		
	*RV_St = 0.0;
	*RV_St_count = 0;

	// allocate otherprocVi vars
	Array<Complex,2> otherprocV1(local_N1,(Nrv[3]/2)+1); 
	Array<Complex,2> otherprocV2(local_N1,(Nrv[3]/2)+1); 
	Array<Complex,2> otherprocV3(local_N1,(Nrv[3]/2)+1);
	Array<Complex,2> temp(local_N1,(Nrv[3]/2)+1); 
	
	// with itself
	otherprocV1 = (*V1r)(Range::all(),pll_index,Range::all());
	otherprocV2 = (*V2r)(Range::all(),pll_index,Range::all());
	otherprocV3 = (*V3r)(Range::all(),pll_index,Range::all());
	Compute_planar_structure_function_unit_anis2(pll_index, my_id, otherprocV1, otherprocV2, otherprocV3);

	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i >= 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id + i) % numprocs;
		
		// first component of plane=pll_index
		temp = (*V1r)(Range::all(),pll_index,Range::all());
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
					 MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV1).data()), data_size, 
					 MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		// second component of plane=pll_index
		temp = (*V2r)(Range::all(),pll_index,Range::all());
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				  MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV2).data()), data_size, 
				 MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		
		// thrid component of plane=pll_index
		temp = (*V3r)(Range::all(),pll_index,Range::all());
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				  MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV3).data()), data_size, 
				 MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		Compute_planar_structure_function_unit_anis2(pll_index, receive_from, 
											otherprocV1, otherprocV2, otherprocV3);			
	} 
}


void RVF::Compute_local_planar_structure_function_anis3
(
	int pll_index,
	int imag_switch
)
{			
	int data_size = local_N1 * Nrv[2];
	int tag = 123;
	int receive_from, send_to;
	MPI_Request request, request2;
	
	// allocate otherprocVi vars
	Array<Real,2> otherprocV1(local_N1,Nrv[2]); 
	Array<Real,2> otherprocV2(local_N1,Nrv[2]);
	Array<Real,2> otherprocV3(local_N1,Nrv[2]);
	Array<Real,2> temp(local_N1,Nrv[2]);
	
	*RV_St = 0.0;
	*RV_St_count = 0;
	
	// with itself
	if (imag_switch == 0) {
		otherprocV1 = real((*V1r)(Range::all(),Range::all(),pll_index));
		otherprocV2 = real((*V2r)(Range::all(),Range::all(),pll_index));
		otherprocV3 = real((*V3r)(Range::all(),Range::all(),pll_index));
		Compute_planar_structure_function_unit_anis3(pll_index, my_id, otherprocV1, otherprocV2, 
													otherprocV3, imag_switch);
	}
	
	else if (imag_switch == 1) {
		otherprocV1 = imag((*V1r)(Range::all(),Range::all(),pll_index));
		otherprocV2 = imag((*V2r)(Range::all(),Range::all(),pll_index));
		otherprocV3 = imag((*V3r)(Range::all(),Range::all(),pll_index));
		Compute_planar_structure_function_unit_anis3(pll_index, my_id, otherprocV1, otherprocV2, 
													 otherprocV3, imag_switch);
	}				
	
	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i >= 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id + i) % numprocs;
		
		// first component of plane=pll_index
		if (imag_switch == 0)
			temp = real((*V1r)(Range::all(),Range::all(),pll_index));
		else if (imag_switch == 1) 
			temp = imag((*V1r)(Range::all(),Range::all(),pll_index));
			
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				   MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
				
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV1).data()), data_size, 
				  MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		
		// second component of plane=pll_index
		if (imag_switch == 0)
			temp = real((*V2r)(Range::all(),Range::all(),pll_index));
		else if (imag_switch == 1) 
			temp = imag((*V2r)(Range::all(),Range::all(),pll_index));
		
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				  MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV2).data()), data_size, 
				  MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		
		// third component of plane=pll_index
		if (imag_switch == 0)
			temp = real((*V3r)(Range::all(),Range::all(),pll_index));
		else if (imag_switch == 1) 
			temp = imag((*V3r)(Range::all(),Range::all(),pll_index));
		
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				  MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocV3).data()), data_size, 
				  MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		Compute_planar_structure_function_unit_anis3(pll_index, receive_from, 
						otherprocV1, otherprocV2, otherprocV3, imag_switch);			
	}

}		
			
			
			
//**************************************************************************************************

// Compute Structure function for each y-z plane (indexed by pll_index), and print in file_out.
// No need of any data transfer across processors.
void RVF::RV_Compute_planar_structure_function_anis1(ofstream& file_out, int pll_index)
{
	*RV_St = 0.0;
	*RV_St_count = 0.0;
	
	Compute_planar_structure_function_unit_anis1(file_out, pll_index);
}

// Compute Structure function for each x-z plane (indexed by pll_index).  Each processor computes
// for a subpart of the plane; these subparts are summed in the master node.
void RVF::RV_Compute_planar_structure_function_anis2(int pll_index)
{
	*RV_St = 0.0;
	*RV_St_count = 0.0;
	
	Compute_local_planar_structure_function_anis2(pll_index);
	
			
	int data_size = (*RV_St).size(); 
	
	MPI_Reduce(reinterpret_cast<Real*>((*RV_St).data()), 
			   reinterpret_cast<Real*>((*RV_St_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	data_size = (*RV_St_count).size();  
	
	MPI_Reduce(reinterpret_cast<int*>((*RV_St_count).data()), 
			   reinterpret_cast<int*>((*RV_St_count_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	
	// Normalize St(..)	
	
	if (my_id == master_id)
		for (int rindex=1; rindex<=RV_structurefn_r_arraysize; rindex++)
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
				if ((*RV_St_count_final)(rindex) > MYEPS2) {
					
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,0) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,0)/(*RV_St_count_final)(rindex);
					
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,1) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,1)/(*RV_St_count_final)(rindex);
				} 
				
}


// Compute Structure function for each x-y plane (indexed by pll_index).  Each processor computes
// for a subpart of the plane; these subparts are summed in the master node.
// Real and imag parts planes are done according to the imag_switch.
void RVF::RV_Compute_planar_structure_function_anis3(int pll_index, int imag_switch)
{
	*RV_St = 0.0;
	*RV_St_count = 0.0;
	
	cout.precision(MY_PRECISION); 
	
	Compute_local_planar_structure_function_anis3(pll_index, imag_switch);
	
	int data_size = (*RV_St).size(); 
	
	MPI_Reduce(reinterpret_cast<Real*>((*RV_St).data()), 
			   reinterpret_cast<Real*>((*RV_St_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	data_size = (*RV_St_count).size();  
	
	MPI_Reduce(reinterpret_cast<int*>((*RV_St_count).data()), 
			   reinterpret_cast<int*>((*RV_St_count_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	
	// Normalize St(..)	
	if (my_id == master_id)
		for (int rindex=1; rindex<=RV_structurefn_r_arraysize; rindex++)
			if ((*RV_St_count_final)(rindex) > MYEPS2)
				for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++) {
					
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,0) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,0)/(*RV_St_count_final)(rindex);
					
					(*RV_St_final)(rindex,q-RV_structurefn_q_min,1) 
					= (*RV_St_final)(rindex,q-RV_structurefn_q_min,1)/(*RV_St_count_final)(rindex);
				}  
	
}


//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************


// For anis1, each plane is independent, and there is no need of a data transfer.
// No mpi_send recv for planar_struct_fn_calc.
void RSF::Compute_planar_structure_function_unit_anis1(ofstream& file_out, int pll_index)
{
	TinyVector<int,3>	j1, j2;
	Real					F1, F2;
	
	int data_size_St, data_size_St_count;
	int tag_St = 111;
	int tag_St_count = 222;
	
	int Pi1 = pll_index;
	int Qi1 = pll_index;
	
	j1(0) = local_N1_start + Pi1;
	j2(0) = local_N1_start + Pi1;
	
	for (int Pi2=0; Pi2<Nrs[2]; Pi2=Pi2+RS_structurefn_neighbour_index) 
		for (int Pi3=0; Pi3<Nrs[3]/2; Pi3=Pi3+RS_structurefn_neighbour_index) {
			
			j1(1) = Pi2; 
			
			RS_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			
			for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RS_structurefn_neighbour_index)
				for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RS_structurefn_neighbour_index) {
					
					j2(1) = Qi2;
					
					// first set
					F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = real((*Fr)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, F1, F2);
					
					// Second set 		
					F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = imag((*Fr)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, F1, F2);
					
					// Third set 		
					F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = real((*Fr)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, F1, F2);
					
					// Fourth set 		
					F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = imag((*Fr)(Qi1,Qi2,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, F1, F2); 
				}		// end of  loops
		}
		
	// if result in master_id, print it.
	// else, send to the master_id.  Master_id will print it
	if (my_id == master_id) {
		file_out << "%% T: pll_index = "  << local_N1_start + pll_index << endl;
		
		for (int rindex=1; rindex<=RS_structurefn_r_arraysize; rindex++) {
			
			file_out	<< rindex << "  ";
			
			for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++) 
				if ((*RS_St_count)(rindex) > MYEPS2) {
					
					(*RS_St)(rindex,q-RS_structurefn_q_min,0) 
					= (*RS_St)(rindex,q-RS_structurefn_q_min,0)/(*RS_St_count)(rindex);
					
					(*RS_St)(rindex,q-RS_structurefn_q_min,1) 
					= (*RS_St)(rindex,q-RS_structurefn_q_min,1)/(*RS_St_count)(rindex);
					
					file_out	<< (*RS_St)(rindex,q-RS_structurefn_q_min,0)	<< " " 
								<< (*RS_St)(rindex,q-RS_structurefn_q_min,1)	<< " ";
				} 
			
			file_out << endl;	
		}
		
		
		for (int source = 1; source <= numprocs-1; source++) {
			
			data_size_St = (*RS_St).size();
			
			MPI_Recv( reinterpret_cast<Real*>((*RS_St).data()), data_size_St, 
					 MPI_Real, source, tag_St, MPI_COMM_WORLD, &status );
			
			data_size_St_count = (*RS_St_count).size(); 		 
			
			MPI_Recv( reinterpret_cast<Real*>((*RS_St_count).data()), data_size_St_count, 
					 MPI_Real, source, tag_St_count, MPI_COMM_WORLD, &status );	
			
			file_out << "%% T: pll_index = " << local_N1*source + pll_index << endl;
			
			for (int rindex=1; rindex<=RS_structurefn_r_arraysize; rindex++) {
				
				file_out	<< rindex << "  ";
				
				for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++) 
					if ((*RS_St_count)(rindex) > MYEPS2) {
						
						(*RS_St)(rindex,q-RS_structurefn_q_min,0) 
						= (*RS_St)(rindex,q-RS_structurefn_q_min,0)/(*RS_St_count)(rindex);
						
						(*RS_St)(rindex,q-RS_structurefn_q_min,1) 
						= (*RS_St)(rindex,q-RS_structurefn_q_min,1)/(*RS_St_count)(rindex);
						
						file_out	<< (*RS_St)(rindex,q-RS_structurefn_q_min,0)	<< " " 
									<< (*RS_St)(rindex,q-RS_structurefn_q_min,1)	<< " ";
					} 
				
				file_out << endl;	
			}							

		}			
	}
	
	else {
		data_size_St = (*RS_St).size();
		
		MPI_Send( reinterpret_cast<Real*>((*RS_St).data()), data_size_St, 
				 MPI_Real, master_id, tag_St, MPI_COMM_WORLD );
		
		data_size_St_count = (*RS_St_count).size(); 
		
		MPI_Send( reinterpret_cast<Real*>((*RS_St_count).data()), data_size_St_count, 
				 MPI_Real, master_id, tag_St_count, MPI_COMM_WORLD );		 
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
		
	
}


void RSF::Compute_planar_structure_function_unit_anis2
(
	 int pll_index,
	 int otherproc, 
	 Array<Complex,2> otherprocF
)
{
	TinyVector<int,3>	j1, j2;
	Real					F1, F2;
	
	int Pi2 = pll_index;
	int Qi2 = pll_index;
	
	j1(1) = pll_index;
	j2(1) = pll_index;
	
	for (int Pi1=0; Pi1<local_N1; Pi1=Pi1+RS_structurefn_neighbour_index) 
		for (int Pi3=0; Pi3<Nrs[3]/2; Pi3=Pi3+RS_structurefn_neighbour_index) {
			
			j1(0) = local_N1_start + Pi1;
			
			RS_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			
			for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RS_structurefn_neighbour_index)
				for (int Qi3=Qi3_start; Qi3<=Qi3_end; Qi3=Qi3+RS_structurefn_neighbour_index) {
					
					j2(0) = otherproc*local_N1 + Qi1; 
					
					// first set
					F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = real((otherprocF)(Qi1,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, F1, F2);
					
					// Second set 		
					F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = imag((otherprocF)(Qi1,Qi3));
					
					j1(2) = 2*Pi3;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, F1, F2);
					
					// Third set 		
					F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = real((otherprocF)(Qi1,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3;
					Compute_dSt(j1, j2, F1, F2);
					
					// Fourth set 		
					F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
					F2 = imag((otherprocF)(Qi1,Qi3));
					
					j1(2) = 2*Pi3+1;
					j2(2) = 2*Qi3+1;
					Compute_dSt(j1, j2, F1, F2); 
				}		// end of  loops
		}
	
}


void RSF::Compute_planar_structure_function_unit_anis3
(
	 int pll_index,
	 int otherproc, 
	 Array<Real,2> otherprocF,
	 int imag_switch
)
{
	TinyVector<int,3>	j1, j2;
	Real					F1, F2;
	
	int Pi3 = pll_index;
	int Qi3 = pll_index;
	
	for (int Pi1=0; Pi1<local_N1; Pi1=Pi1+RS_structurefn_neighbour_index) 
		for (int Pi2=0; Pi2<Nrs[2]/2; Pi2=Pi2+RS_structurefn_neighbour_index) {
			
			j1(0) = local_N1_start + Pi1;
			j1(1) = Pi2;
			
			RS_structure_fn_Q_ends(Pi1, Pi2, Pi3);
			
			for (int Qi1=Qi1_start; Qi1<=Qi1_end; Qi1=Qi1+RS_structurefn_neighbour_index)
				for (int Qi2=Qi2_start; Qi2<=Qi2_end; Qi2=Qi2+RS_structurefn_neighbour_index) {
					
					j2(0) = otherproc*local_N1 + Qi1;
					j2(1) = Qi2;
					
					if (imag_switch == 0) {
						F1 = real((*Fr)(Pi1,Pi2,Pi3)); 
						F2 = otherprocF(Qi1,Qi2);
						
						j1(2) = 2*Pi3;
						j2(2) = 2*Qi3;
						Compute_dSt(j1, j2, F1, F2);
					}
					
					else if (imag_switch == 1) {
						F1 = imag((*Fr)(Pi1,Pi2,Pi3)); 
						F2 = otherprocF(Qi1,Qi2);
						
						j1(2) = 2*Pi3+1;
						j2(2) = 2*Qi3+1;
						Compute_dSt(j1, j2, F1, F2);  
					}
				} // end of Q loop	
		} // end of P loop	
	
}




//*********************************************************************************************


void RSF::Compute_local_planar_structure_function_anis2(int pll_index)
{
	
	int data_size = 2* local_N1 * (Nrs[3]/2 + 1);
	int tag = 123;
	int receive_from, send_to;
	MPI_Request request, request2;
	
	// allocate otherprocVi vars
	Array<Complex,2> otherprocF(local_N1,(Nrs[3]/2)+1);
	Array<Complex,2> temp(local_N1,(Nrs[3]/2)+1); 
	
	*RS_St = 0.0;
	*RS_St_count = 0;
	
	// with itself
	otherprocF = (*Fr)(Range::all(),pll_index,Range::all());
	Compute_planar_structure_function_unit_anis2(pll_index, my_id, otherprocF);
	
	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i > 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id+i) % numprocs;
		
		temp = (*Fr)(Range::all(),pll_index,Range::all());
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				  MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);			 
		 
		MPI_Irecv( reinterpret_cast<Real*>((otherprocF).data()), data_size, 
				  MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);

		Compute_planar_structure_function_unit_anis2(pll_index, receive_from, otherprocF);			
	}
}


void RSF::Compute_local_planar_structure_function_anis3
(
	 int pll_index,
	 int imag_switch
)
{			
	int data_size = local_N1 * Nrs[2];
	int tag = 123;
	int receive_from, send_to;
	MPI_Request request, request2;	
	
	// allocate otherprocVi vars
	Array<Real,2> otherprocF(local_N1,Nrs[2]); 
	Array<Real,2> temp(local_N1,Nrs[2]);
	
	*RS_St = 0.0;
	*RS_St_count = 0;
	
	// with itself
	if (imag_switch == 0) {
		otherprocF = real((*Fr)(Range::all(),Range::all(),pll_index));
		Compute_planar_structure_function_unit_anis3(pll_index, my_id, otherprocF, imag_switch);
	}
	
	else if (imag_switch == 1) {
		otherprocF = imag((*Fr)(Range::all(),Range::all(),pll_index));
		Compute_planar_structure_function_unit_anis3(pll_index, my_id, otherprocF, imag_switch);
	}				
	
	for (int i=1; i<=dist_farthest_proc; i++) {
		if (my_id -i > 0)
			send_to = my_id - i; 
		else
			send_to = my_id - i + numprocs;
		
		receive_from = (my_id+i) % numprocs;
		
		if (imag_switch == 0)
			temp = real((*Fr)(Range::all(),Range::all(),pll_index));
		else if (imag_switch == 1) 
			temp = imag((*Fr)(Range::all(),Range::all(),pll_index));
		
		MPI_Isend( reinterpret_cast<Real*>((temp).data()), data_size, 
				 MPI_Real, send_to, tag, MPI_COMM_WORLD, &request);
		
		MPI_Irecv( reinterpret_cast<Real*>((otherprocF).data()), data_size, 
					 MPI_Real, receive_from, tag, MPI_COMM_WORLD, &request2);
		
		MPI_Wait(&request2, &status);
		
		Compute_planar_structure_function_unit_anis3(pll_index, receive_from, otherprocF, imag_switch);			
	}
	
}		


//**************************************************************************************************

// Computes Structure function for each y-z plane (indexed by pll_index), and 
// print in file_out. No need of any data transfer across processors.
void RSF::RS_Compute_planar_structure_function_anis1(ofstream& file_out, int pll_index)
{
	*RS_St = 0.0;
	*RS_St_count = 0.0;
	
	Compute_planar_structure_function_unit_anis1(file_out, pll_index);
}


// Compute Structure function for each x-z plane (indexed by pll_index).  
// Each processor computes for a subpart of the plane; 
// these subparts are summed in the master node.
void RSF::RS_Compute_planar_structure_function_anis2(int pll_index)
{
	*RS_St = 0.0;
	*RS_St_count = 0.0;
	
	Compute_local_planar_structure_function_anis2(pll_index);
	
	int data_size = (*RS_St).size(); 
	
	MPI_Reduce(reinterpret_cast<Real*>((*RS_St).data()), 
			   reinterpret_cast<Real*>((*RS_St_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	data_size = (*RS_St_count).size();  
	
	MPI_Reduce(reinterpret_cast<int*>((*RS_St_count).data()), 
			   reinterpret_cast<int*>((*RS_St_count_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	
	// Normalize St(..)	
	
	if (my_id == master_id)
		for (int rindex=1; rindex<=RS_structurefn_r_arraysize; rindex++)
			if ((*RS_St_count_final)(rindex) > MYEPS2)
				for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++)
					
					(*RS_St_final)(rindex,q-RS_structurefn_q_min) 
					= (*RS_St_final)(rindex,q-RS_structurefn_q_min)/(*RS_St_count_final)(rindex);
	
}


// Compute Structure function for each x-y plane (indexed by pll_index).  
// Each processor computes for a subpart of the plane; 
// these subparts are summed in the master node.
// Real and imag parts planes are done according to the imag_switch.
void RSF::RS_Compute_planar_structure_function_anis3(int pll_index, int imag_switch)
{
	*RS_St = 0.0;
	*RS_St_count = 0.0;
	
	Compute_local_planar_structure_function_anis3(pll_index, imag_switch);
	
	int data_size = (*RS_St).size(); 
	
	MPI_Reduce(reinterpret_cast<Real*>((*RS_St).data()), 
			   reinterpret_cast<Real*>((*RS_St_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	data_size = (*RS_St_count).size();  
	
	MPI_Reduce(reinterpret_cast<int*>((*RS_St_count).data()), 
			   reinterpret_cast<int*>((*RS_St_count_final).data()), 
			   data_size, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
	// Normalize St(..)	
	
	if (my_id == master_id)
		for (int rindex=1; rindex<=RS_structurefn_r_arraysize; rindex++) {
			if ((*RS_St_count_final)(rindex) > MYEPS2)
				for (int q=RS_structurefn_q_min; q<=RS_structurefn_q_max; q++)
					(*RS_St_final)(rindex,q-RS_structurefn_q_min) 
					= (*RS_St_final)(rindex,q-RS_structurefn_q_min)/(*RS_St_count_final)(rindex);
		}
	
}




//*********************************   End of struct_fn.cc *************************************




