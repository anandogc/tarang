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

/*! \file field_basic.h
 * 
 * @brief Some basic Array operations: Array_real_mult, Output_asreal, Model_initial_shell_spectrum
 *
 * @author  M. K. Verma
 * @date Sept 2008
 *
 * @bugs Out of date functions Read_data_split_arrays_MPI, Write_data_split_arrays_MPI contain bugs.
 */ 
 

#ifndef _ARRAY_OPS_H
#define _ARRAY_OPS_H

#include "def_vars.h"
#include "Global_extern_vars.h"
#include "basic_extern_vars.h"

//*********************************************************************************************

/** @brief Multiply real arrays A and B term by term, and put the result in C.  
 *
 *  The multiplication is done in each process independently 
 *	without any requirement of communications.
 *
 * @param  A  Real array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$. 
 * @param  B  Real array \f$ B [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$. 
 * @param  N[]  The size of the arrays A, B.
 * 
 * @return   C(i,j,k) = A(i,j,k)*B(i,j,k).  
 */

class ArrayOps 
{
public:

	static void Real_space_multiply
	(
		Array<DP,3> A, Array<DP,3> B, Array<DP,3> C
	);

	static void Real_space_divide
	(
		Array<DP,3> A, Array<DP,3> B, Array<DP,3> C
	);

	template<class T, int rank>
	static void Print_array_all_procs(Array<T, rank> A){
		for (int i=0; i<numprocs; i++){
			if (my_id == i)
				cout << "my_id, Array = " << my_id << ", " << A << endl;
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
};

#endif

//*********************************   End of array_basic.h  ***********************************
