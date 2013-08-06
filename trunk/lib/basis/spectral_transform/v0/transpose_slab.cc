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


/*! \file	basis_basicfn.cc
 * 
 * @brief basic functions definitions that are common to all basis functions.
 *
 * @sa basis_basicfn.h
 *
 * @version 4.0 Parallel version (v0)
 * @author  M. K. Verma
 * @date	Sept 2008
 * @bug		The Transpose for 2D is not working at present.. it must transport 
 *				real array.
 */ 


#include "spectral_transform.h"


//
//
//******************************************************************************************

void SpectralTransform::Transpose_array_SLAB(Array<complx,3> A, Array<complx,3> B, char X1, char X2, char X3, char fixedaxis)
{
    
    Array<complx,2> *A2dTR, *A2dTR_recv;
    int A_first_dim, A_second_dim, ATR_local_size;

    // A(X,Y,Z) or A(Y,X<Z).. 
    if (fixedaxis == 'Z') {
        if (numprocs == 1)
            B = A.transpose(1,0,2);
        //    B = A.transpose(2,1,0);
        
        else {
            
            if (X1 == 'X') {
                A_first_dim = local_Nx;
                A_second_dim = Ny;
                ATR_local_size = local_Ny;
                
                A2dTR = &(slab_transpose.AxyTR);
                A2dTR_recv = &(slab_transpose.AxyTR_recv);
            }
            
            else if (X1 == 'Y') {
                A_first_dim = local_Ny;
                A_second_dim = Nx;
                ATR_local_size = local_Nx;
                
                A2dTR = &(slab_transpose.AyxTR);
                A2dTR_recv = &(slab_transpose.AyxTR_recv);
            }
            
           int data_size = 2* local_Nx * local_Ny;								
            // 2 for complex to double
        
        
            // upto i3=Nz/2 for FFF, SFF, SSF;  i3=Nz/2-1 for SSS
            for (int i3=0; i3<A.extent(2); i3++) {
               *A2dTR = A(Range::all(), Range::all(), i3).transpose(1,0);

                // The jth block sent from process i is received by process j and 
                // is placed in the ith block of recvbuf. 
                MPI_Alltoall(reinterpret_cast<DP*>(A2dTR->data()), data_size,  MPI_DP, 
                             reinterpret_cast<DP*>(A2dTR_recv->data()),
                             data_size,   MPI_DP, MPI_COMM_WORLD);
                
                for (int source = 0; source < numprocs; source++) 
                    B(Range::all(), Range(source*A_first_dim,(source+1)*A_first_dim-1), i3)
                    = (*A2dTR_recv)(Range(source*ATR_local_size,(source+1)*ATR_local_size-1), Range::all());
            } 
        }
    }
    
    
    // A(Z,Y,X) or A(X,Y,Z)
    else if (fixedaxis == 'Y') {
        if (numprocs == 1)
            B = A.transpose(2,1,0);
        
        else { // CHANGE.. Here for ch
            if (X1 == 'X') {
                A_first_dim = local_Nx;
                A_second_dim = Ny;
                ATR_local_size = local_Ny;
                
                A2dTR = &(slab_transpose.AxyTR);
                A2dTR_recv = &(slab_transpose.AxyTR_recv);
            }
            
            else if (X1 == 'Y') {
                A_first_dim = local_Ny;
                A_second_dim = Nx;
                ATR_local_size = local_Nx;
                
                A2dTR = &(slab_transpose.AyxTR);
                A2dTR_recv = &(slab_transpose.AyxTR_recv);
            }
            
            int data_size = 2* local_Nx * local_Ny;
            // 2 for complex to double
            
            
            // upto i3=Nz/2 for FFF, SFF, SSF;  i3=Nz/2-1 for SSS
            for (int i3=0; i3<A.extent(2); i3++) {
                *A2dTR = A(Range::all(), Range::all(), i3).transpose(1,0);
                
                // The jth block sent from process i is received by process j and
                // is placed in the ith block of recvbuf.
                MPI_Alltoall(reinterpret_cast<DP*>(A2dTR->data()), data_size,  MPI_DP,
                             reinterpret_cast<DP*>(A2dTR_recv->data()),
                             data_size,   MPI_DP, MPI_COMM_WORLD);
                
                for (int source = 0; source < numprocs; source++)
                    B(Range::all(), Range(source*A_first_dim,(source+1)*A_first_dim-1), i3)
                    = (*A2dTR_recv)(Range(source*ATR_local_size,(source+1)*ATR_local_size-1), Range::all());
            }
        }
    }
  
}

//****************************************************************************************

// 2D case: Axz to Azx; here N[2]=1
// N[3]/2+1 is completely divisible by numproc;  Choose N[3]=2^n-2.
// local_N3 = (N[3]/2+1)/numprocs
void SpectralTransform::Transpose_array_SLAB(Array<complx,2> A, Array<complx,2> B, char X1, char X2)
{
		
    Array<complx,2> *A2dTR, *A2dTR_recv;
    int A_first_dim, A_second_dim, ATR_local_size;
    
    if (numprocs == 1) 
        B = A.transpose(2,1,0);
    
    else {
        if (X1 == 'X') {
            A_first_dim = local_Nx;
            ATR_local_size = local_Nz;
            
            A2dTR = &(slab_transpose.AxzTR);
            A2dTR_recv = &(slab_transpose.AxzTR_recv);
        }
        
        else if (X1 == 'Z') {
            A_first_dim = local_Nz;
            A_second_dim = Nx;
            ATR_local_size = local_Nx;
            
			A2dTR = &(slab_transpose.AzxTR);
            A2dTR_recv = &(slab_transpose.AzxTR_recv);
        }
        
        int data_size = 2* local_Nx * local_Nz;								
        // 2 for complex to double
		
        *A2dTR = A(Range::all(),Range::all()).transpose(1,0);
        
        MPI_Alltoall(reinterpret_cast<DP*>(A2dTR->data()), data_size, MPI_DP, 
                     reinterpret_cast<DP*>(A2dTR_recv->data()),
                     data_size,  MPI_DP, MPI_COMM_WORLD);
        
        for (int source = 0; source < numprocs; source++)									
            B(Range::all(),Range(source*A_first_dim,(source+1)*A_first_dim-1))
            = (*A2dTR_recv)(Range(source*ATR_local_size,(source+1)*ATR_local_size-1), Range::all());

    }
}

//******************************************************************************************

/*
void SpectralTransform::Transpose_array_SLAB_new(Array<complx,3> A, Array<complx,3> B, char X1, char X2, char X3, char fixedaxis)
{
    
    Array<complx,2> *A2d, *A2d_recv;
    
    if (numprocs == 1)
        Atr = A.transpose(1,0,2);
    
    else {
        if (X1 == 'X') {
            A_first_dim = local_Nx;
            ATR_first_dim = local_Ny;
            A2d = &(global.temp_array.transpose.Axy);
            A2d_recv = &(global.temp_array.transpose.Axy_recv);
        }
        
        else if (X1 == 'Y') {
            A_first_dim = local_Ny;
            ATR_first_dim = local_Nx;
            A2d = &(global.temp_array.transpose.Ayx);
            A2d_recv = &(global.temp_array.transpose.Ayx_recv);
        } 
        
        int data_size = 2* local_Nx*local_Ny;								
        // 2 for complex to double
        
        for (int i3=0; i3<=N[3]/2; i3++) {
            A2d = A(Range::all(), Range::all(), i3);
            
            MPI_Alltoall(reinterpret_cast<DP*>(*A2d.data()), data_size, MPI_DP, reinterpret_cast<DP*>(*A2d_recv.data()), data_size, MPI_DP, MPI_COMM_WORLD);
            
            for (int source = 0; source < numprocs; source++) 
                B(Range::all(), Range(source*A_first_dim,(source+1)*A_first_dim-1),Range::all())
                = (*A2d_recv)(Range::all(),Range(source*ATR_first_dim,(source+1)*ATR_first_dim-1),Range::all()).transpose(1,0,2);
        }
    }
}   


/*****************************************************************************************

// 2D case: Axz to Azx; here N[2]=1
// N[3]/2+1 is completely divisible by numproc;  Choose N[3]=2^n-2.
// local_N3 = (N[3]/2+1)/numprocs
void SpectralTransform::Transpose_array2d_new(Array<complx,3> A, Array<complx,3> B, char X1, char X2)
{
    
    Array<complx,2> *A2d_recv, *A2dTR, *A2dTR_afterAll2all;
    
    if (numprocs == 1) 
        B = A.transpose(2,1,0);
    
    else {
        if (X1 == 'X') {
            A_first_dim = local_Nx;
            ATR_first_dim = local_Nz;
            A2d_recv = &global.temp_array.transpose.Axz;
        }
        
        else if (X1 == 'Z') {
            A_first_dim = local_Nz;
            ATR_first_dim = local_Nx;
            A2d_recv = &global.temp_array.transpose.Azx;
        }
        
        int data_size = 2* local_Nx * local_Nz;								
        // 2 for complex to double
        
        MPI_Alltoall(reinterpret_cast<DP*>(A.data()), data_size, MPI_DP, reinterpret_cast<DP*>(A2d_recv.data()), data_size,  MPI_DP, MPI_COMM_WORLD);
        
        for (int source = 0; source < numprocs; source++)									
            B(Range::all(), 0, Range(source*A_first_dim,(source+1)*A_first_dim-1)) = A2d_recv(Range::all(),Range(source*ATR_first_dim,(source+1)*ATR_first_dim-1)).transpose(1,0);
    }
}
*/

//***********************************  End of transpose.cc ********************************



