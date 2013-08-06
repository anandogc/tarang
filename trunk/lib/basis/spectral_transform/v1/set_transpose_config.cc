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

template <class T>
 void SpectralTransform::Set_transpose_config(int send_config_id, int recv_config_id, TransposeConfig& config, MPI_Comm communicator){

	int num_blocks;
	int extent;
	int stride;
	int block_lengths[2];
	
	MPI_Aint displacements[2];
	MPI_Datatype types[2];

	int unit = sizeof(T)/sizeof(DP);
	int Z_length;

	config.communicator = communicator;
	MPI_Comm_size(communicator, &config.numprocs);

	MPI_Datatype MPI_Vector_unit;	// DP or complx

	MPI_Datatype MPI_Vector_block;	//slab
	MPI_Datatype MPI_Struct_block;	//slab

	MPI_Datatype MPI_Vector_block_x;  //pencil
	MPI_Datatype MPI_Struct_block_x;  //pencil

    
	MPI_Datatype MPI_Vector_block_z;  //pencil
	MPI_Datatype MPI_Struct_block_z;  //pencil

	num_blocks = 1;
	extent = unit;
	stride = unit;

	MPI_Type_vector(num_blocks,extent,stride,MPI_DP,&MPI_Vector_unit);
	MPI_Type_commit(&MPI_Vector_unit);

	for (int i=0; i<2; i++){

		int& count = (i==0) ? config.sendcount : config.recvcount;
		MPI_Datatype& type = (i==0) ? config.sendtype : config.recvtype;
		int config_id = (i==0) ? send_config_id : recv_config_id;

		switch (config_id){
			case 0:	//3D slab horizontal
				Z_length = 2*shape_horizontal_array_3d(1);
				count = local_Nx*local_Ny*Z_length;
				type = MPI_DP;

				break;

			case 1: //3D slab vertical
				Z_length = 2*shape_vertical_array_3d(1)/unit;

				num_blocks = local_Ny*Z_length;
				extent = local_Nx;
				stride = Nx;

				MPI_Type_vector(num_blocks,extent,stride,MPI_Vector_unit,&MPI_Vector_block);
				MPI_Type_commit(&MPI_Vector_block);


				block_lengths[0]=1;
				block_lengths[1]=1;
				
				displacements[0]=0;
				displacements[1]=local_Nx*sizeof(T);

				types[0]=MPI_Vector_block;
				types[1]=MPI_UB;

				MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block);
				MPI_Type_commit(&MPI_Struct_block);

				count = 1;
				type = MPI_Struct_block;

				break;

			case 2:	//2D slab horizontal
				count = 2*local_Nx*local_Nz;
				type = MPI_DP;

				break;

			case 3: //2D slab vertical

				num_blocks = 2*local_Nz/unit;
				extent = local_Nx;
				stride = Nx;

				MPI_Type_vector(num_blocks,extent,stride,MPI_Vector_unit,&MPI_Vector_block);
				MPI_Type_commit(&MPI_Vector_block);


				block_lengths[0]=1;
				block_lengths[1]=1;
				
				displacements[0]=0;
				displacements[1]=local_Nx*sizeof(T);

				types[0]=MPI_Vector_block;
				types[1]=MPI_UB;

				MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block);
				MPI_Type_commit(&MPI_Struct_block);

				count = 1;
				type = MPI_Struct_block;

				break;


			case 4: //pencil along X in VERT communicator

			    num_blocks = 2*local_Ny_vert*local_Nz_hor/unit;
			    extent = local_Nx_vert;
			    stride = Nx;

			    MPI_Type_vector(num_blocks,extent,stride,MPI_Vector_unit,&MPI_Vector_block_x);
			    MPI_Type_commit(&MPI_Vector_block_x);


			    block_lengths[0]=1;
			    block_lengths[1]=1;
			    
			    displacements[0]=0;
			    displacements[1]=local_Nx_vert*sizeof(T);

			    types[0]=MPI_Vector_block_x;
			    types[1]=MPI_UB;
				
			    MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block_x);
			    MPI_Type_commit(&MPI_Struct_block_x);

			    count = 1;
			    type = MPI_Struct_block_x;

			    break;

			case 5: //pencil along X in HOR communicator
				//not used
				break;

			case 6: //pencil along Y in VERT communicator
				count = 2*local_Nx_vert*local_Ny_vert*local_Nz_hor;
				type = MPI_DP;

				break;

			case 7: //pencil along Y in HOR communicator
				count = 2*local_Nx_vert*local_Ny_hor*local_Nz_hor;
				type = MPI_DP;

				break;

			case 8: //pencil along Z in VERT communicator
				Z_length = 2*shape_z_array_3d(1)/unit;
				num_blocks = local_Ny_hor;
				extent = 2*local_Nx_vert*local_Nz_vert/unit;
				stride = local_Nx_vert*Z_length;
				
				MPI_Type_vector(num_blocks,extent,stride,MPI_DP,&MPI_Vector_block_z);
				MPI_Type_commit(&MPI_Vector_block_z);
				
				block_lengths[0]=1;
				block_lengths[1]=1;
				
				displacements[0]=0;
				displacements[1]=2*local_Nx_vert*local_Nz_vert*sizeof(T)/unit;
				
				types[0]=MPI_Vector_block_z;
				types[1]=MPI_UB;
				
				MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block_z);
				MPI_Type_commit(&MPI_Struct_block_z);

				count = 1;
				type = MPI_Struct_block_z;

				break;

			case 9: //pencil aling Z in HOR communicator
			    Z_length = 2*shape_z_array_3d(1)/unit;
			    num_blocks = local_Ny_hor;
				extent = 2*local_Nx_vert*local_Nz_hor/unit;
				stride = local_Nx_vert*Z_length;
				
				MPI_Type_vector(num_blocks,extent,stride,MPI_Vector_unit,&MPI_Vector_block_z);
				MPI_Type_commit(&MPI_Vector_block_z);
				
				block_lengths[0]=1;
				block_lengths[1]=1;
				
				displacements[0]=0;
				displacements[1]=2*local_Nx_vert*local_Nz_hor*sizeof(T)/unit;
				
				types[0]=MPI_Vector_block_z;
				types[1]=MPI_UB;
				
				MPI_Type_struct(2, block_lengths, displacements, types, &MPI_Struct_block_z);
				MPI_Type_commit(&MPI_Struct_block_z);

				count = 1;
				type = MPI_Struct_block_z;

				break;

		}
	}
 }