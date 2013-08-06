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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _CONFIGURE_TRANSPOSE_H
#define _CONFIGURE_TRANSPOSE_H

#include "spectral_def_vars.h"


//*********************************************************************************************	

struct TransposeConfig{
	int sendcount;
	MPI_Datatype sendtype;

	int recvcount;
	MPI_Datatype recvtype;

	MPI_Comm communicator;

	int numprocs;

	TransposeConfig();
	TransposeConfig(int send_count, MPI_Datatype send_type, int recv_count, MPI_Datatype recv_type, MPI_Comm communicator, int numprocs);
	void configure(int send_count, MPI_Datatype send_type, int recv_count, MPI_Datatype recv_type, MPI_Comm communicator, int numprocs);
	TransposeConfig conj();
};

#endif


//******************************** End of field_basic.h  **************************************


