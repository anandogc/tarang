
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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 


#include "gpu.h"

//****************************************************************************************					
 
void Forward_transform()
{	
}


void Inverse_transform()
{
}

void Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}

void Add_Xderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}


void Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}

void Add_Yderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}


void Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}

void Add_Zderiv(Array<complx,3> A, Array<complx,3> B)
{
	
}

void Real_space_multiply_gpu(Array<complx,3> A, Array<complx,3> B, Array<complx,3> C)
{
	
}

void Print_large_Fourier_elements(Array<complx,3> A, string array_name)
{
	
}

void Array_divide_ksqr(Array<complx,3> A)
{
	
}

void Array_exp_ksqr(Array<complx,3> A, DP factor)
{
	
}

void Array_exp_ksqr(Array<complx,3> A, DP factor, DP hyper_factor, int hyper_exponent)
{
	
}

void Compute_divergence(Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, Array<complx,3> div, string field_or_nlin, DP &total_abs_div, bool print_switch)
{
	
}

void Compute_total_energy(V1, V2, V3)
{
	
}




//********************************** End of Ifluid_main.cc ************************************



