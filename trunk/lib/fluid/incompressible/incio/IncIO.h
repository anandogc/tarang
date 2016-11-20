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


/*! \file  IO.h
 * 
 *	@brief  Class declaration of input/output function.
 * 
 *
 *	@author  M. K. Verma
 *	@version 2.0 MPI
 *	@date feb 2012
 *
 * @bug  No known bugs
 */
 
#ifndef _H_FluidIO_incompres
#define _H_FluidIO_incompres

#include "fluid_base.h"
#include "Pressure.h"
#include "EnergyTr.h"
#include "FluidIO.h"

//*********************************************************************************************

	// Array read write
	// write each component separately, time.V1.hdf, time.V2.hdf... etc.

class FluidIO_incompress : public FluidIO
{
    
public:
	EnergyTr *energyTr;
	
	ofstream		pressure_file;
	ofstream		pressure_spectrum_file;
	
	ofstream		flux_file;
	ofstream		shell_to_shell_file;		
	ofstream		ring_to_ring_file;
	ofstream		cylindrical_ring_to_ring_file;	
	
public:
	~FluidIO_incompress(); 
	
	void Init_energy_transfer();
	void Open_files();
    void Close_files();
	
	void Output_all_inloop(FluidVF& U, Pressure& P, FluidVF& helicalU);
	void Output_all_inloop(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_all_inloop(FluidVF& U, FluidVF& W, Pressure& P, FluidVF& helicalU, FluidVF& helicalW);
	void Output_all_inloop(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);
	void Output_all_inloop(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P);
	void Output_all_inloop(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P);
	
	void Output_last(FluidVF& U, Pressure& P, FluidVF& helicalU);
	void Output_last(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_last(FluidVF& U, FluidVF& W, Pressure& P, FluidVF& helicalU, FluidVF& helicalW);
	void Output_last(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);
	void Output_last(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C, Pressure& P);
	void Output_last(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P);
	
	
	void Output_pressure(Pressure& P);
	void Output_pressure_spectrum(Pressure& P);
		
	void Output_flux(FluidVF& U, Pressure& P, FluidVF& helicalU);
	void Output_flux(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_flux(FluidVF& U, FluidVF& W, Pressure& P, FluidVF& helicalU, FluidVF& helicalW);
	void Output_flux(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);
	
  	void Output_shell_to_shell(FluidVF& U, Pressure& P, FluidVF& helicalU);
	void Output_shell_to_shell(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_shell_to_shell(FluidVF& U, FluidVF& W, Pressure& P, FluidVF& helicalU, FluidVF& helicalW);
	void Output_shell_to_shell(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);

  	void Output_ring_to_ring(FluidVF& U, Pressure& P);
	void Output_ring_to_ring(FluidVF& U, Pressure& P, FluidVF& helicalU);
	void Output_ring_to_ring(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_ring_to_ring_scalar(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_ring_to_ring_RBC(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_ring_to_ring(FluidVF& U, FluidVF& W, Pressure& P, FluidVF& helicalU, FluidVF& helicalW);
	void Output_ring_to_ring(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);
	
	void Output_cylindrical_ring_to_ring(FluidVF& U, Pressure& P);
	void Output_cylindrical_ring_to_ring(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_cylindrical_ring_to_ring_scalar(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_cylindrical_ring_to_ring_RBC(FluidVF& U, FluidSF& T, Pressure& P);
	void Output_cylindrical_ring_to_ring(FluidVF& U, FluidVF& W, Pressure& P);
	void Output_cylindrical_ring_to_ring(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P);
};


#endif
//========================= Class declaration of IO ends ============================== 
 
