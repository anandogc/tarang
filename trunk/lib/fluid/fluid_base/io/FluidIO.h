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
 
//*********************************************************************************************
#ifndef _H_FLUIDIO
#define _H_FLUIDIO

#include "fluid_base.h" 
	// Array read write
	// write each component separately, time.V1.hdf, time.V2.hdf... etc.
class FluidIO 
{	
public:
	ifstream		field_in_file;
	ifstream		force_field_in_file;
	
	
	ofstream		field_out_file;
	ofstream		field_frequent_out_file;
	ofstream		field_out_reduced_file;
	ofstream		realfield_out_file;
	
	ofstream		field_k_out_file;
	ofstream		field_r_out_file;
	
	ofstream		global_file;
    ofstream        global_real_space_file;
	ofstream		spectrum_file;
	ofstream		ring_spectrum_file;
	ofstream		cylindrical_ring_spectrum_file;
	ofstream		structure_fn_file;
	ofstream		planar_structure_fn_file;
	ofstream		pressure_file;
    
    ofstream        Tk_spectrum_file;
    ofstream        Tk_ring_spectrum_file;
    ofstream        Tk_cylindrical_ring_spectrum_file;
    
    ofstream        misc_file;
	
	
public:
	
	void Open_base_files();
    void Close_base_files();
	
	void Init_cond_complex_field(FluidVF& U);
	void Init_cond_complex_field(FluidVF& U, FluidSF& T);
	void Init_cond_scalar_complex_field(FluidVF& U, FluidSF& T);
	void Init_cond_RBC_complex_field(FluidVF& U, FluidSF& T);
	void Init_cond_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Init_cond_complex_field(FluidVF& U, FluidVF& W);
	void Init_cond_complex_field(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Init_cond_reduced_complex_field(FluidVF& U);
	void Init_cond_reduced_complex_field(FluidVF& U, FluidSF& T);
	void Init_cond_reduced_complex_field(FluidVF& U, FluidVF& W);
	void Init_cond_reduced_complex_field(FluidVF& U, FluidVF& W, FluidSF& T);
	void Init_cond_reduced_complex_field_scalar(FluidVF& U, FluidSF& T);
	void Init_cond_reduced_complex_field_RBC(FluidVF& U, FluidSF& T);
	void Init_cond_reduced_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	
	void Init_cond_real_field(FluidVF& U);
	void Init_cond_real_field(FluidVF& U, FluidSF& T);
	void Init_cond_real_field_scalar(FluidVF& U, FluidSF& T);
	void Init_cond_real_field_RBC(FluidVF& U, FluidSF& T);
	void Init_cond_real_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Init_cond_real_field(FluidVF& U, FluidVF& W);
	void Init_cond_real_field(FluidVF& U, FluidVF& W, FluidSF& T);
	void Read_init_cond(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Init_cond_modes_Chebyshev(FluidVF&U);
	void Init_cond_modes_scalar_Chebyshev(FluidVF&U, FluidSF& T);
	void Init_cond_modes_RBC_Chebyshev(FluidVF&U, FluidSF& T);
	
	void Init_cond_modes(FluidVF &U);
	void Init_cond_modes(FluidVF&U, FluidSF& T);
	void Init_cond_modes_scalar(FluidVF&U, FluidSF& T);
	void Init_cond_modes_RBC(FluidVF&U, FluidSF& T);
	void Init_cond_modes(FluidVF& U, FluidVF& W);
	void Init_cond_modes(FluidVF& U, FluidVF& W, FluidSF& T);
	void Init_cond_modes(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	void Init_cond_modes(FluidVF&U, FluidSF& T1, FluidSF& T2);

	void Read_init_cond(FluidVF& U);	
	void Read_init_cond(FluidVF& U, FluidSF& T);
	void Read_init_cond(FluidVF& U, FluidVF& W);
	void Read_init_cond(FluidVF& U, FluidVF& W, FluidSF& T);
	void Read_init_cond(FluidVF& U, FluidSF& T1, FluidSF& T2);
    void Read_init_cond(FluidSF& T);

	void Setup_Taylor_Green_field(FluidVF& U, int k0, Real amp);
	void Setup_ABC_field(FluidVF& U, int k0, Real amp, Real A, Real B, Real C);
	void Setup_SIX_MODE_field(FluidVF& U, int k0, Real amp101, Real amp011, Real amp112, Real h);
	void Model_initial_using_shell_spectrum_Pope(Real dissipation_coefficient, Real epsilon, Array<Real,1> Sk);
	void Model_initial_using_shell_spectrum_Corrsin(Real a, Array<Real,1> Sk);
	
	void Init_cond_Taylor_Green(FluidVF& U);
	void Init_cond_Taylor_Green(FluidVF& U, FluidSF& T);
	void Init_cond_Taylor_Green_scalar(FluidVF& U, FluidSF& T);
	void Init_cond_Taylor_Green_RBC(FluidVF& U, FluidSF& T);
	void Init_cond_Taylor_Green(FluidVF& U, FluidVF& W);
	void Init_cond_Taylor_Green(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Init_cond_ABC(FluidVF& U);
	void Init_cond_ABC(FluidVF& U, FluidSF& T);
	void Init_cond_ABC(FluidVF& U, FluidVF& W);
	void Init_cond_ABC(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Init_cond_Rayleigh_Taylor(FluidVF&U, FluidSF& T);
	
	void Init_cond_dynamo_full_velocity_field(FluidVF& U, FluidVF& W);
	
	void Init_cond_DYNAMO_SIX_MODE(FluidVF& U, FluidVF& W);
	
	void Put_vector_amp_phase_comp_conj(FluidVF U, int lx, int ly, int lz,  Real amp, Real phase1, Real phase2, Real phase3);
	
	void Put_scalar_amp_phase_comp_conj(FluidSF& T, int lx, int ly, int lz, Real amp, Real phase);
	
	void Initialize_using_energy_helicity_spectrum(FluidVF& U, Real spectrum_amp, Real hk_by_kek);
	void Initialize_using_energy_helicity_spectrum(FluidSF& T, Real spectrum_amp);
	
	void Init_cond_energy_helicity_spectrum(FluidVF& U);
	void Init_cond_energy_helicity_spectrum(FluidVF& U, FluidSF& T);
	void Init_cond_energy_helicity_spectrum_scalar(FluidVF& U, FluidSF& T);
	void Init_cond_energy_helicity_spectrum_RBC(FluidVF& U, FluidSF& T);
	void Init_cond_energy_helicity_spectrum(FluidVF& U, FluidVF& W);
	void Init_cond_energy_helicity_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
    
    void Init_cond_channel_flow(FluidVF& U);
    
    void Init_cond_vortex(FluidVF& U);
    void Init_cond_vortex(FluidVF& U, FluidSF& T);
	void Init_cond_vortex(FluidVF& U, FluidVF& W);
    void Init_cond_vortex(FluidVF& U, FluidVF& W, FluidSF& T);
    void Init_cond_vortex(FluidSF& T);
	
		///
	
	bool Global_data_buffer_full();
	bool Spectral_space_buffer_full();
	bool Real_space_buffer_full();
	
	void Output_complex_field(FluidVF& U);
	void Output_complex_field(FluidVF& U, FluidSF& T);
	void Output_complex_field_scalar(FluidVF& U, FluidSF& T);
	void Output_complex_field_RBC(FluidVF& U, FluidSF& T);
	void Output_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_complex_field(FluidVF& U, FluidVF& W);
	void Output_complex_field(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_complex_field_frequent(FluidVF& U);
	void Output_complex_field_frequent(FluidVF& U, FluidSF& T);
	void Output_complex_field_scalar_frequent(FluidVF& U, FluidSF& T);
	void Output_complex_field_RBC_frequent(FluidVF& U, FluidSF& T);
	void Output_complex_field_frequent(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_complex_field_frequent(FluidVF& U, FluidVF& W);
	void Output_complex_field_frequent(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_reduced_complex_field(FluidVF& U);
	void Output_reduced_complex_field(FluidVF& U, FluidSF& T);
	void Output_reduced_complex_field_scalar(FluidVF& U, FluidSF& T);
	void Output_reduced_complex_field_RBC(FluidVF& U, FluidSF& T);
	void Output_reduced_complex_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_reduced_complex_field(FluidVF& U, FluidVF& W);
	void Output_reduced_complex_field(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_real_field(FluidVF& U);
	void Output_real_field(FluidVF& U, FluidSF& T);
	void Output_real_field(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_real_field(FluidVF& U, FluidVF& W);
	void Output_real_field(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_field_k(FluidVF& U);
	void Output_field_k(FluidVF& U, FluidSF& T);
	void Output_field_k(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_field_k(FluidVF& U, FluidVF& W);
	void Output_field_k(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_field_r(FluidVF& U);
	void Output_field_r(FluidVF& U, FluidSF& T);
	void Output_field_r(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_field_r(FluidVF& U, FluidVF& W);
	void Output_field_r(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_global(FluidVF& U);
	void Output_global(FluidVF& U, FluidSF& T);
	void Output_global_scalar(FluidVF &U, FluidSF& T);
	void Output_global_RBC(FluidVF& U, FluidSF& T);
	void Output_global(FluidVF& U, FluidVF& W);
	void Output_global(FluidVF& U, FluidVF& W, FluidSF& T);
	void Output_global(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	void Output_global(FluidVF &U, FluidSF& T1, FluidSF& T2);
	
	void Output_cout(FluidVF& U);
	void Output_cout(FluidVF& U, FluidSF& T);
	void Output_cout(FluidVF &U, FluidSF& T1, FluidSF& T2);
	void Output_cout(FluidVF& U, FluidVF& W);
	void Output_cout(FluidVF& U, FluidVF& W, FluidSF& T);
	void Output_cout(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C);
	
	void Output_shell_spectrum(FluidVF& U);
	void Output_shell_spectrum(FluidVF& U, FluidSF& T);
	void Output_shell_spectrum_scalar(FluidVF& U, FluidSF& T);
	void Output_shell_spectrum(FluidVF& U, FluidSF& T1, FluidSF& T2);
	void Output_shell_spectrum_RBC(FluidVF& U, FluidSF& T);
	void Output_shell_spectrum(FluidVF& U, FluidVF& W);
	void Output_shell_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_ring_spectrum(FluidVF& U);
	void Output_ring_spectrum(FluidVF& U, FluidSF& T);
	void Output_ring_spectrum_scalar(FluidVF& U, FluidSF& T);
	void Output_ring_spectrum_RBC(FluidVF& U, FluidSF& T);
	void Output_ring_spectrum(FluidVF& U, FluidVF& W);
	void Output_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_cylindrical_ring_spectrum(FluidVF& U);
	void Output_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T);
	void Output_cylindrical_ring_spectrum_scalar(FluidVF& U, FluidSF& T);
	void Output_cylindrical_ring_spectrum_RBC(FluidVF& U, FluidSF& T);
	void Output_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W);
	void Output_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_structure_fn(FluidVF& U);
	void Output_structure_fn(FluidVF& U, FluidSF& T);
	void Output_structure_fn(FluidVF& U, FluidVF& W);
	void Output_structure_fn(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_planar_structure_fn(FluidVF& U);
	void Output_planar_structure_fn(FluidVF& U, FluidSF& T);
	void Output_planar_structure_fn(FluidVF& U, FluidVF& W);
	void Output_planar_structure_fn(FluidVF& U, FluidVF& W, FluidSF& T);
    
    void Output_Tk_shell_spectrum(FluidVF& U);
	void Output_Tk_shell_spectrum(FluidVF& U, FluidSF& T);;
	void Output_Tk_shell_spectrum(FluidVF& U, FluidVF& W);
	void Output_Tk_shell_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_Tk_ring_spectrum(FluidVF& U);
	void Output_Tk_ring_spectrum(FluidVF& U, FluidSF& T);
	void Output_Tk_ring_spectrum(FluidVF& U, FluidVF& W);
	void Output_Tk_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_Tk_cylindrical_ring_spectrum(FluidVF& U);
	void Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidSF& T);
	void Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W);
	void Output_Tk_cylindrical_ring_spectrum(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_global_real_space(FluidVF& U);
	void Output_global_real_space(FluidVF& U, FluidSF& T);
	void Output_global_real_space_scalar(FluidVF& U, FluidSF& T);
	void Output_global_real_space_RBC(FluidVF& U, FluidSF& T);
	void Output_global_real_space(FluidVF& U, FluidVF& W);
	void Output_global_real_space(FluidVF& U, FluidVF& W, FluidSF& T);
	
	void Output_cout_real_space(FluidVF& U);
	void Output_cout_real_space(FluidVF& U, FluidSF& T);
	void Output_cout_real_space(FluidVF& U, FluidVF& W);
	void Output_cout_real_space(FluidVF& U, FluidVF& W, FluidSF& T);


	//Interprets the given array as 2d and prints the fastest(last) index in a row
	template <int N_rank>
	void Print_array(ofstream& file, string label, Array<Real,N_rank> B)
	{
		if (master){
			
			Array<Real,N_rank> A(B.shape());
			A=B;
			
			int dim1=A.extent(A.rank()-1);	//extent of last index
			int dim2=A.size()/dim1;

			std::ostream_iterator<Real> file_it (file,"\t");

			file << "%%" << label << "\n";

			for (int i=0; i<dim2; i++){	//print dim2 rows
				std::copy ( A.data()+i*dim1, A.data()+(i+1)*dim1, file_it );	//print dim1 cols
				file << '\n';
			}
		}
	}
	
	//Interprets the given array as 2d and prints the fastest(last) index in a row
	template <int N_rank>
	void Print_array(ofstream& file, string label, Array<Real,N_rank> B1, Array<Real,N_rank> B2)
	{
		if (master){
			Array<Real,N_rank> A1(B1.shape());
			Array<Real,N_rank> A2(B2.shape());
			
			A1=B1;
			A2=B2;

			int dim1=A1.extent(A1.rank()-1);	//extent of last index
			int dim2=A1.size()/dim1;

			std::ostream_iterator<Real> file_it (file,"\t");

			file << "%%" << label << "\n";

			for (int i=0; i<dim2; i++){	//print dim2 rows
				std::copy ( A1.data()+i*dim1, A1.data()+(i+1)*dim1, file_it );	//print dim1 cols
				file << '\n';
			}

			if (!global.program.two_dimension){
				for (int i=0; i<dim2; i++){	//print dim2 rows
					std::copy ( A2.data()+i*dim1, A2.data()+(i+1)*dim1, file_it );	//print dim1 cols
					file << '\n';
				}
			}
		}
	}

	template <int N_rank>
	void Print_array(ofstream& file, string label, Array<Real,N_rank> B1, Array<Real,N_rank> B2, Array<Real,N_rank> B3)
	{
		if (master){
			Array<Real,N_rank> A1(B1.shape());
			Array<Real,N_rank> A2(B2.shape());
			Array<Real,N_rank> A3(B3.shape());
			
			A1=B1;
			A2=B2;
			A3=B3;

			int dim1=A1.extent(A1.rank()-1);	//extent of last index (num cols)
			int dim2=A1.size()/dim1;			//num rows

			std::ostream_iterator<Real> file_it (file,"\t");

			file << "%%" << label << "\n";

			for (int i=0; i<dim2; i++){	//print dim2 rows
				std::copy ( A1.data()+i*dim1, A1.data()+(i+1)*dim1, file_it );	//print dim1 cols
				file << '\n';
			}


			if (!global.program.two_dimension){
				for (int i=0; i<dim2; i++){	//print dim2 rows
					std::copy ( A2.data()+i*dim1, A2.data()+(i+1)*dim1, file_it );	//print dim1 cols
					file << '\n';
				}
			}

			for (int i=0; i<dim2; i++){	//print dim2 rows
				std::copy ( A3.data()+i*dim1, A3.data()+(i+1)*dim1, file_it );	//print dim1 cols
				file << '\n';
			}
		}
	}

	
};


#endif
//========================= Class declaration of IO ends ============================== 
 
