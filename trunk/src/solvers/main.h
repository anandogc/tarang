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


/*! \file  main.h
 * 
 *	@brief  Main file for Tarang (MPI). Calls various modules.
 *
 *	@author  M. K. Verma
 *	@version 4.0  MPI
 *	@date Sept 2008
 */

//********************************************************************************************* 
#include "fluid.h"   
// Contains IncFlow declarations (which contains fields +fourier+sincosfour defs)
														
int Ifluid_main();

	// int Iscalar_main();

int Iscalar_main();

int MRBC_main();

int IMHD_main();

int IMHDscalar_main();

int IMHDastro_main();

int GP_main();

void process_mem_usage(Real& vm_usage, Real& resident_set);	
//********************************** End of main.h ********************************************


