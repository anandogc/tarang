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

/*! \file  output_main.cc
 * 
 * @brief  Output_all_inloop & Output_field_k_inloop. 
 *
 *	In Output_all_inloop <BR>
 *	Output_global:	Total energy, dissipation etc. <BR>
 *	Output_cout:	Total energy, etc. on cout.	<BR>
 *	Output_field:	Field in k space. <BR?
 *	Output_field_reduced: Field in smaller box. <BR>
 *	Output_realfield: Field in r space <BR>
 *	Output_shell_spectrum: Shell spectrum <BR>
 *	Output_flux:	flux <BR>
 *	Output_shell_to_shell: shell-to-shell energy tr <BR>
 *	Output_field_frequent: frequent saving of field vars; Overwrites the field. <BR>
 *	Output_ring_spectrum: ring <BR>
 *	Output_ring_to_ring: ring-to-ring enregy transfer <BR>
 *	Output_cylindrical_ring_spectrum: cylinderical ring spectrum <BR>
 *	Output_cylindrical_ring_to_ring: Cylindrical ring-to-ring energy transfer <BR>
 *  Output_structire_fn: Structure function <BR>
 *	Output_planar_structure_fn: Planar structure function. <BR>
 *
 *	In Output_field_k_inloop: <BR>
 *	Output_field_k:	Outputs \f$ V(k) \f$ at specified k's.
 *	Output_field_r: Outputs \f$ V(r) \f$ at specified r's.
 *	Output_Skpq: Outputs S(k|p|q) for specified triads.
 *
 *	Output_pressure_spectrum_inloop: Outputs pressure spectrum. <BR>
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 


#include "IncIO.h"   


//*********************************************************************************************

/** @brief Outputs various things in main loop.
 *
 *	Outputs every save_next time. Saving takes place immediately after global.time.now
 *	crosses save_next. <BR>
 *	save_next is updated to save_next + save_interval.
 *
 *	save_interval etc. are saved in Time class.
 */

void FluidIO_incompress::Output_all_inloop(FluidVF& U, Pressure& P)
{

	if (global.time.now >= global.io.time.global_save_next) {	 
		Output_global(U); 
		global.io.time.global_save_next +=  global.io.time.global_save_interval; 
	}
	
	if (global.time.now >= global.io.time.complex_field_save_next) {  
		Output_complex_field(U); 
		global.io.time.complex_field_save_next += global.io.time.complex_field_save_interval;
	}	
	
	if (global.time.now >= global.io.time.field_frequent_save_next) {\
		Output_complex_field_frequent(U); 
		global.io.time.field_frequent_save_next += global.io.time.field_frequent_save_interval;
	}	
	
	if (global.time.now >= global.io.time.field_reduced_save_next) {
		Output_reduced_complex_field(U);
		global.io.time.field_reduced_save_next += global.io.time.field_reduced_save_interval;
	}	
	

	if (global.time.now >= global.io.time.spectrum_save_next) {
		Output_shell_spectrum(U);	
		global.io.time.spectrum_save_next += global.io.time.spectrum_save_interval;
	}		
					
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.flux_save_next)) {
		Output_flux(U, P);
		global.io.time.flux_save_next += global.io.time.flux_save_interval;
	}	
			
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.shell_to_shell_save_next)) {
		Output_shell_to_shell(U, P);
		global.io.time.shell_to_shell_save_next += global.io.time.shell_to_shell_save_interval;
	}	
	

	if ((global.time.now >= global.io.time.ring_spectrum_save_next) && (global.spectrum.ring.turnon)) {  
		Output_ring_spectrum(U);
		global.io.time.ring_spectrum_save_next += global.io.time.ring_spectrum_save_interval;
	}	
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.ring_to_ring_save_next)  && (global.energy_transfer.ring_to_ring.turnon))  {
		Output_ring_to_ring(U,P);
		global.io.time.ring_to_ring_save_next += global.io.time.ring_to_ring_save_interval;
	}	
	
	if ((global.time.now >= global.io.time.cylindrical_ring_spectrum_save_next) && (global.spectrum.cylindrical_ring.turnon)) {  
		Output_cylindrical_ring_spectrum(U); 
		global.io.time.cylindrical_ring_spectrum_save_next += global.io.time.cylindrical_ring_spectrum_save_interval;
	}
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.cylindrical_ring_to_ring_save_next) && (global.energy_transfer.cylindrical_ring_to_ring.turnon))  {  
		Output_cylindrical_ring_to_ring(U, P);
		global.io.time.cylindrical_ring_to_ring_save_next += global.io.time.cylindrical_ring_to_ring_save_interval;
	} 

	if (global.time.now >= global.io.time.structure_fn_save_next) {
	/*	if (global.structure_fn.box_switch == 1)
			Output_structure_fn();			
		
		if (global.structure_fn.planar_switch == 1)
			Output_planar_structure_fn(); */
			
		global.io.time.structure_fn_save_next += global.io.time.structure_fn_save_interval;
	} 
	
	if (global.time.now >= global.io.time.cout_save_next) {  
		Output_cout(U); 
		global.io.time.cout_save_next += global.io.time.cout_save_interval; 
	}
}



void FluidIO_incompress::Output_last(FluidVF& U, Pressure& P)
{
	
	if (global.io.time.global_save_last)
		Output_global(U);
	
	if (global.io.time.complex_field_save_last)
		Output_complex_field(U); 
	
	if (global.io.time.field_frequent_save_last) {
		Real total_abs_div;
		U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
		Output_complex_field_frequent(U);
	}
	
	if (global.io.time.field_reduced_save_last)
		Output_reduced_complex_field(U);

	if (global.io.time.real_field_save_last) {
		U.Inverse_transform();
		Output_real_field(U);
	}	
	if (global.io.time.spectrum_save_last)
		Output_shell_spectrum(U);	
	
	if ((global.energy_transfer.turnon) && (global.io.time.flux_save_last))
		Output_flux(U, P);
	
	if ((global.energy_transfer.turnon) && (global.io.time.shell_to_shell_save_last))
		Output_shell_to_shell(U, P);
	
	
	if ((global.spectrum.ring.turnon) && (global.io.time.ring_spectrum_save_last)) 
		Output_ring_spectrum(U); 
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.ring_to_ring.turnon) && (global.io.time.ring_to_ring_save_last))
		Output_ring_to_ring(U, P);
	
	if ((global.spectrum.cylindrical_ring.turnon) && (global.io.time.cylindrical_ring_spectrum_save_last))
		Output_cylindrical_ring_spectrum(U);
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.cylindrical_ring_to_ring.turnon) && (global.io.time.cylindrical_ring_to_ring_save_last))
		Output_cylindrical_ring_to_ring(U, P);
	
	if (global.io.time.structure_fn_save_last) {
		/*	if (global.structure_fn.box_switch == 1)
		 Output_structure_fn();			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(); */
	} 
	
	if (global.io.time.cout_save_last)
		Output_cout(U); 
}


//*********************************************************************************************
// with scalar
//
void FluidIO_incompress::Output_all_inloop(FluidVF& U, FluidSF& T, Pressure& P)
{

	if (global.time.now >= global.io.time.global_save_next) {
		Output_global(U, T); 
		global.io.time.global_save_next +=  global.io.time.global_save_interval; 
	}
		
	if (global.time.now >= global.io.time.complex_field_save_next) {
		Output_complex_field(U, T);  
		global.io.time.complex_field_save_next += global.io.time.complex_field_save_interval;
	}	
	
	if (global.time.now >= global.io.time.field_frequent_save_next) {  
		Output_complex_field_frequent(U, T); 
		global.io.time.field_frequent_save_next += global.io.time.field_frequent_save_interval;
	}

	if (global.time.now >= global.io.time.field_reduced_save_next) {
		Output_reduced_complex_field(U, T);
		global.io.time.field_reduced_save_next += global.io.time.field_reduced_save_interval;
	}
		
	if (global.time.now >= global.io.time.spectrum_save_next) {
		Output_shell_spectrum(U, T);	
		global.io.time.spectrum_save_next += global.io.time.spectrum_save_interval;
	}		

	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.flux_save_next)) {
		Output_flux(U, T, P);
		global.io.time.flux_save_next += global.io.time.flux_save_interval;
	}	

	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.shell_to_shell_save_next)) {
		Output_shell_to_shell(U, T, P);
		global.io.time.shell_to_shell_save_next += global.io.time.shell_to_shell_save_interval;
	}	
	   
	if ((global.time.now >= global.io.time.ring_spectrum_save_next)  && (global.spectrum.ring.turnon)) {  
		Output_ring_spectrum(U, T); 
		global.io.time.ring_spectrum_save_next += global.io.time.ring_spectrum_save_interval;
	}	

	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.ring_to_ring_save_next) 
				&& (global.energy_transfer.ring_to_ring.turnon))  {  
		Output_ring_to_ring(U, T, P); 
		global.io.time.ring_to_ring_save_next += global.io.time.ring_to_ring_save_interval;
	}

	if ((global.time.now >= global.io.time.cylindrical_ring_spectrum_save_next) && (global.spectrum.cylindrical_ring.turnon)) {  
		Output_cylindrical_ring_spectrum(U, T); 
		global.io.time.cylindrical_ring_spectrum_save_next += global.io.time.cylindrical_ring_spectrum_save_interval;
	}

	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.cylindrical_ring_to_ring_save_next) && (global.energy_transfer.cylindrical_ring_to_ring.turnon)) {  
		Output_cylindrical_ring_to_ring(U, T, P); 
		global.io.time.cylindrical_ring_to_ring_save_next += global.io.time.cylindrical_ring_to_ring_save_interval;
	} 
	
/*	if (global.time.now >= global.io.time.structure_fn_save_next) {
		if (Rglobal.structure_fn.box_switch == 1)
			Output_structure_fn(T);			
		
		if (global.structure_fn.planar_switch == 1)
			Output_planar_structure_fn(T); 
			
		global.io.time.structure_fn_save_next += global.io.time.structure_fn_save_interval;
	}  */
	
	if (global.time.now >= global.io.time.cout_save_next) {  
		Output_cout(U, T); 
		global.io.time.cout_save_next += global.io.time.cout_save_interval; 
	}
}


void FluidIO_incompress::Output_last(FluidVF& U, FluidSF& T, Pressure& P)
{
	
	if (global.io.time.global_save_last)
		Output_global(U, T);
	
	if (global.io.time.complex_field_save_last)
		Output_complex_field(U, T); 
	
	if (global.io.time.field_frequent_save_last) 
		Output_complex_field_frequent(U, T);
	
	if (global.io.time.field_reduced_save_last)
		Output_reduced_complex_field(U, T);

	if (global.io.time.real_field_save_last) {
		U.Inverse_transform();
		T.Inverse_transform();
		Output_real_field(U, T);
	}
	
	if (global.io.time.spectrum_save_last)
		Output_shell_spectrum(U, T);	
	
	if ((global.energy_transfer.turnon) && (global.io.time.flux_save_last))
		Output_flux(U, T, P);			
	
	if ((global.energy_transfer.turnon) && (global.io.time.shell_to_shell_save_last))
		Output_shell_to_shell(U, T, P);
	
	
	if ((global.spectrum.ring.turnon) && (global.io.time.ring_spectrum_save_last)) 
		Output_ring_spectrum(U, T); 
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.ring_to_ring.turnon) && (global.io.time.ring_to_ring_save_last))
		Output_ring_to_ring(U, T, P); 
	
	if ((global.spectrum.cylindrical_ring.turnon) && (global.io.time.cylindrical_ring_spectrum_save_last))
		Output_cylindrical_ring_spectrum(U, T);
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.cylindrical_ring_to_ring.turnon) && (global.io.time.cylindrical_ring_to_ring_save_last))
		Output_cylindrical_ring_to_ring(U, T, P);
	
	if (global.io.time.structure_fn_save_last) {
		/*	if (global.structure_fn.box_switch == 1)
		 Output_structure_fn();			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(); */
	} 
	
	if (global.io.time.cout_save_last)
		Output_cout(U, T); 
}


//*********************************************************************************************
// with vector
//
void FluidIO_incompress::Output_all_inloop(FluidVF& U, FluidVF& W, Pressure& P)
{

	if (global.time.now >= global.io.time.global_save_next) {	 
		Output_global(U, W);
		global.io.time.global_save_next +=  global.io.time.global_save_interval; 
	}
	
	if (global.time.now >= global.io.time.complex_field_save_next)  {  
		Output_complex_field(U, W);  
		global.io.time.complex_field_save_next += global.io.time.complex_field_save_interval;
	}	

	if (global.time.now >= global.io.time.field_frequent_save_next) {  
		Output_complex_field_frequent(U, W); 
		global.io.time.field_frequent_save_next += global.io.time.field_frequent_save_interval;
	}
	
	if (global.time.now >= global.io.time.field_reduced_save_next) {	
		Output_reduced_complex_field(U, W);
		global.io.time.field_reduced_save_next += global.io.time.field_reduced_save_interval;
	}	

				
	if (global.time.now >= global.io.time.spectrum_save_next) {
		Output_shell_spectrum(U, W);	
		global.io.time.spectrum_save_next += global.io.time.spectrum_save_interval;
	}		
	
				
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.flux_save_next)) {
		Output_flux(U, W, P);			
		global.io.time.flux_save_next += global.io.time.flux_save_interval;
	}	
			
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.shell_to_shell_save_next)) {
		Output_shell_to_shell(U, W, P);
		global.io.time.shell_to_shell_save_next += global.io.time.shell_to_shell_save_interval;
	}	
	
	if ((global.time.now >= global.io.time.ring_spectrum_save_next) && (global.spectrum.ring.turnon)) {  
		Output_ring_spectrum(U, W); 
		global.io.time.ring_spectrum_save_next += global.io.time.ring_spectrum_save_interval;
	}	
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.ring_to_ring_save_next) 
				&& (global.energy_transfer.ring_to_ring.turnon))  {  
		Output_ring_to_ring(U, W, P); 
		global.io.time.ring_to_ring_save_next += global.io.time.ring_to_ring_save_interval;
	}
	
	if ((global.time.now >= global.io.time.cylindrical_ring_spectrum_save_next) && (global.spectrum.cylindrical_ring.turnon)) {  
		Output_cylindrical_ring_spectrum(U, W); 
		global.io.time.cylindrical_ring_spectrum_save_next += global.io.time.cylindrical_ring_spectrum_save_interval;
	}
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.cylindrical_ring_to_ring_save_next) && (global.energy_transfer.cylindrical_ring_to_ring.turnon)) {  
		Output_cylindrical_ring_to_ring(U, W, P); 
		global.io.time.cylindrical_ring_to_ring_save_next += global.io.time.cylindrical_ring_to_ring_save_interval;
	}
		
	if (global.time.now >= global.io.time.structure_fn_save_next) {
	//	if (global.structure_fn.box_switch == 1)
	//		Output_structure_fn(W);			
		
	//	if (global.structure_fn.planar_switch == 1)
	//		Output_planar_structure_fn(W);
			
		global.io.time.structure_fn_save_next += global.io.time.structure_fn_save_interval;
	}
	
		
	if (global.time.now >= global.io.time.cout_save_next) {  
		Output_cout(U, W); 
		global.io.time.cout_save_next += global.io.time.cout_save_interval; 
	}	
}


void FluidIO_incompress::Output_last(FluidVF& U, FluidVF& W, Pressure& P)
{
	
	if (global.io.time.global_save_last)
		Output_global(U, W);
	
	if (global.io.time.complex_field_save_last)
		Output_complex_field(U, W); 
	
	if (global.io.time.field_frequent_save_last) 
		Output_complex_field_frequent(U, W);
	
	if (global.io.time.field_reduced_save_last)
		Output_reduced_complex_field(U, W);

	if (global.io.time.real_field_save_last) {
		U.Inverse_transform();
		W.Inverse_transform();
		Output_real_field(U, W);
	}	
	if (global.io.time.spectrum_save_last)
		Output_shell_spectrum(U, W);	
	
	if ((global.energy_transfer.turnon) && (global.io.time.flux_save_last))
		Output_flux(U, W, P);			
	
	if ((global.energy_transfer.turnon) && (global.io.time.shell_to_shell_save_last))
		Output_shell_to_shell(U, W, P);
	
	
	if ((global.spectrum.ring.turnon) && (global.io.time.ring_spectrum_save_last)) 
		Output_ring_spectrum(U, W); 
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.ring_to_ring.turnon) && (global.io.time.ring_to_ring_save_last))
		Output_ring_to_ring(U, W, P); 
	
	if ((global.spectrum.cylindrical_ring.turnon) && (global.io.time.cylindrical_ring_spectrum_save_last))
		Output_cylindrical_ring_spectrum(U, W);
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.cylindrical_ring_to_ring.turnon) && (global.io.time.cylindrical_ring_to_ring_save_last))
		Output_cylindrical_ring_to_ring(U, W, P);
	
	if (global.io.time.structure_fn_save_last) {
		/*	if (global.structure_fn.box_switch == 1)
		 Output_structure_fn();			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(); */
	} 
	
	if (global.io.time.cout_save_last)
		Output_cout(U, W); 
}

//*********************************************************************************************
// W and T
//
void FluidIO_incompress::Output_all_inloop(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{

	if (global.time.now >= global.io.time.global_save_next) {	 
		Output_global(U, W, T); 
		global.io.time.global_save_next +=  global.io.time.global_save_interval; 
	}
			
	if (global.time.now >= global.io.time.complex_field_save_next) {  
		Output_complex_field(U, W, T); 
		global.io.time.complex_field_save_next += global.io.time.complex_field_save_interval;
	}	
	
	if (global.time.now >= global.io.time.field_frequent_save_next) {  
		Output_complex_field_frequent(U, W, T); 
		global.io.time.field_frequent_save_next += global.io.time.field_frequent_save_interval;
	}
	
	
	if (global.time.now >= global.io.time.field_reduced_save_next)  {	
		Output_reduced_complex_field(U, W, T);
		global.io.time.field_reduced_save_next += global.io.time.field_reduced_save_interval;
	}	
	
/*	if (global.time.now >= global.io.time.real_field_save_next) {
		Output_real_field(U, W, T);
		global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}	*/
			
	if (global.time.now >= global.io.time.spectrum_save_next) {
		Output_shell_spectrum(U, W, T);	
		global.io.time.spectrum_save_next += global.io.time.spectrum_save_interval;
	}		
			
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.flux_save_next)) {
		Output_flux(U, W, T, P);			
		global.io.time.flux_save_next += global.io.time.flux_save_interval;
	}	
			
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.shell_to_shell_save_next))
	{
		Output_shell_to_shell(U, W, T, P);
		global.io.time.shell_to_shell_save_next += global.io.time.shell_to_shell_save_interval;
	}	

	if ((global.time.now >= global.io.time.ring_spectrum_save_next)  && (global.spectrum.ring.turnon)) {  
		Output_ring_spectrum(U, W, T); 
		global.io.time.ring_spectrum_save_next += global.io.time.ring_spectrum_save_interval;
	}	
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.ring_to_ring_save_next) 
				&& (global.energy_transfer.ring_to_ring.turnon)) {  
		Output_ring_to_ring(U, W, T, P); 
		global.io.time.ring_to_ring_save_next += global.io.time.ring_to_ring_save_interval;
	}
	
	if ((global.time.now >= global.io.time.cylindrical_ring_spectrum_save_next) && (global.spectrum.cylindrical_ring.turnon)) {  
		Output_cylindrical_ring_spectrum(U, W, T); 
		global.io.time.cylindrical_ring_spectrum_save_next += global.io.time.cylindrical_ring_spectrum_save_interval;
	}
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.cylindrical_ring_to_ring_save_next) && (global.energy_transfer.cylindrical_ring_to_ring.turnon)) {  
		Output_cylindrical_ring_to_ring(U, W, T, P); 
		global.io.time.cylindrical_ring_to_ring_save_next += global.io.time.cylindrical_ring_to_ring_save_interval;
	}
	
	if (global.time.now >= global.io.time.structure_fn_save_next)
	{
	/*	if (global.structure_fn.box_switch == 1)
			Output_structure_fn(W, T);			
		
		if (global.structure_fn.planar_switch == 1)
			Output_planar_structure_fn(W, T); */
			
		global.io.time.structure_fn_save_next += global.io.time.structure_fn_save_interval;
	}
	
	if (global.time.now >= global.io.time.cout_save_next) {  
		Output_cout(U, W, T); 
		global.io.time.cout_save_next += global.io.time.cout_save_interval; 
	}
}

void FluidIO_incompress::Output_last(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P)
{
	
	if (global.io.time.global_save_last)
		Output_global(U, W, T);
	
	if (global.io.time.complex_field_save_last)
		Output_complex_field(U, W, T); 
	
	if (global.io.time.field_frequent_save_last) 
		Output_complex_field_frequent(U, W, T);
	
	if (global.io.time.field_reduced_save_last)
		Output_reduced_complex_field(U, W, T);

	if (global.io.time.real_field_save_last) {
		U.Inverse_transform();
		W.Inverse_transform();
		T.Inverse_transform();
		Output_real_field(U, W, T);
	}	
	if (global.io.time.spectrum_save_last)
		Output_shell_spectrum(U, W, T);	
	
	if ((global.energy_transfer.turnon) && (global.io.time.flux_save_last))
		Output_flux(U, W, T, P);			
	
	if ((global.energy_transfer.turnon) && (global.io.time.shell_to_shell_save_last))
		Output_shell_to_shell(U, W, T, P);
	
	
	if ((global.spectrum.ring.turnon) && (global.io.time.ring_spectrum_save_last)) 
		Output_ring_spectrum(U, W, T); 
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.ring_to_ring.turnon) && (global.io.time.ring_to_ring_save_last))
		Output_ring_to_ring(U, W, T, P); 
	
	if ((global.spectrum.cylindrical_ring.turnon) && (global.io.time.cylindrical_ring_spectrum_save_last))
		Output_cylindrical_ring_spectrum(U, W, T);
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.cylindrical_ring_to_ring.turnon) && (global.io.time.cylindrical_ring_to_ring_save_last))
		Output_cylindrical_ring_to_ring(U, W, T, P);
	
	if (global.io.time.structure_fn_save_last) {
		/*	if (global.structure_fn.box_switch == 1)
		 Output_structure_fn();			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(); */
	} 
	
	if (global.io.time.cout_save_last)
		Output_cout(U, W, T); 
}


//*********************************************************************************************
// W and T
//
void FluidIO_incompress::Output_all_inloop(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P)
{
	
	if (global.time.now >= global.io.time.global_save_next) {
		Output_global(U, W, T, C);
		global.io.time.global_save_next +=  global.io.time.global_save_interval;
	}
	
	if (global.time.now >= global.io.time.cout_save_next) {
		Output_cout(U, W, T, C);
		global.io.time.cout_save_next += global.io.time.cout_save_interval;
	}
}

void FluidIO_incompress::Output_last(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P)
{
}


//*********************************************************************************************
// with scalar
//
void FluidIO_incompress::Output_all_inloop(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P)
{
	
	if (global.time.now >= global.io.time.global_save_next) {	 
		Output_global(U, T1, T2); 
		global.io.time.global_save_next +=  global.io.time.global_save_interval; 
	}
	
	if (global.time.now >= global.io.time.complex_field_save_next) {  
		Output_complex_field(U, T1, T2);  
		global.io.time.complex_field_save_next += global.io.time.complex_field_save_interval;
	}	
	
	if (global.time.now >= global.io.time.field_frequent_save_next) {  
		Output_complex_field_frequent(U, T1, T2); 
		global.io.time.field_frequent_save_next += global.io.time.field_frequent_save_interval;
	}
	
	if (global.time.now >= global.io.time.field_reduced_save_next) {	
		Output_reduced_complex_field(U, T1, T2);
		global.io.time.field_reduced_save_next += global.io.time.field_reduced_save_interval;
	}	
	
/*	if (global.time.now >= global.io.time.real_field_save_next) {
		Output_real_field(U, T1, T2);
		global.io.time.real_field_save_next += global.io.time.real_field_save_interval;
	}	*/
	
	if (global.time.now >= global.io.time.spectrum_save_next) {
		Output_shell_spectrum(U, T1, T2);	
		global.io.time.spectrum_save_next += global.io.time.spectrum_save_interval;
	}		
	/*
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.flux_save_next)) {
		Output_flux(U, T1, T2);			
		global.io.time.flux_save_next += global.io.time.flux_save_interval;
	}	
	
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.shell_to_shell_save_next)) {
		Output_shell_to_shell(U, T);
		global.io.time.shell_to_shell_save_next += global.io.time.shell_to_shell_save_interval;
	}	
	
	
	if ((global.time.now >= global.io.time.ring_spectrum_save_next)  && (global.spectrum.ring.turnon)) {  
		Output_ring_spectrum(U, T); 
		global.io.time.ring_spectrum_save_next += global.io.time.ring_spectrum_save_interval;
	}	
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.ring_to_ring_save_next) 
		&& (global.energy_transfer.ring_to_ring.turnon))  {  
		Output_ring_to_ring(U, T); 
		global.io.time.ring_to_ring_save_next += global.io.time.ring_to_ring_save_interval;
	}
	
	if ((global.time.now >= global.io.time.cylindrical_ring_spectrum_save_next) && (global.spectrum.cylindrical_ring.turnon)) {  
		Output_cylindrical_ring_spectrum(U, T); 
		global.io.time.cylindrical_ring_spectrum_save_next += global.io.time.cylindrical_ring_spectrum_save_interval;
	}
	
	if ((global.energy_transfer.turnon) && (global.time.now >= global.io.time.cylindrical_ring_to_ring_save_next) && (global.energy_transfer.cylindrical_ring_to_ring.turnon)) {  
		Output_cylindrical_ring_to_ring(U, T); 
		global.io.time.cylindrical_ring_to_ring_save_next += global.io.time.cylindrical_ring_to_ring_save_interval;
	}
	
	if (global.time.now >= global.io.time.structure_fn_save_next) {
		if (Rglobal.structure_fn.box_switch == 1)
		 Output_structure_fn(T);			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(T); 
		
		global.io.time.structure_fn_save_next += global.io.time.structure_fn_save_interval;
	}
	*/
	
	if (global.time.now >= global.io.time.cout_save_next) {  
		Output_cout(U, T1, T2); 
		global.io.time.cout_save_next += global.io.time.cout_save_interval; 
	}
}


void FluidIO_incompress::Output_last(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P)
{
	
	if (global.io.time.global_save_last)
		Output_global(U, T1, T2);
	
	if (global.io.time.complex_field_save_last)
		Output_complex_field(U, T1, T2);
	
	if (global.io.time.field_frequent_save_last) 
		Output_complex_field_frequent(U, T1, T2);
	
	if (global.io.time.field_reduced_save_last)
		Output_reduced_complex_field(U, T1, T2);

	if (global.io.time.real_field_save_last) {
		U.Inverse_transform();
		T1.Inverse_transform();
		T2.Inverse_transform();
		Output_real_field(U, T1, T2);
	}	
	if (global.io.time.spectrum_save_last)
		Output_shell_spectrum(U, T1, T2);	
	
/*	if ((global.energy_transfer.turnon) && (global.io.time.flux_save_last))
		Output_flux(U, T1, T2);			
	
	if ((global.energy_transfer.turnon) && (global.io.time.shell_to_shell_save_last))
		Output_shell_to_shell(U, T1, T2);
	
	
	if ((global.spectrum.ring.turnon) && (global.io.time.ring_spectrum_save_last)) 
		Output_ring_spectrum(U, T1, T2); 
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.ring_to_ring.turnon) && (global.io.time.ring_to_ring_save_last))
		Output_ring_to_ring(U, T1, T2);
	
	if ((global.spectrum.cylindrical_ring.turnon) && (global.io.time.cylindrical_ring_spectrum_save_last))
		Output_cylindrical_ring_spectrum(U, T1, T2);
	
	if ((global.energy_transfer.turnon) && (global.energy_transfer.cylindrical_ring_to_ring.turnon) && (global.io.time.cylindrical_ring_to_ring_save_last))
		Output_cylindrical_ring_to_ring(U, T1, T2);
	
	if (global.io.time.structure_fn_save_last) {
		if (global.structure_fn.box_switch == 1)
		 Output_structure_fn();			
		 
		 if (global.structure_fn.planar_switch == 1)
		 Output_planar_structure_fn(); 
	} */
	
	if (global.io.time.cout_save_last)
		Output_cout(U, T1, T2);
}





//********************************  End of output_main.cc ************************************



