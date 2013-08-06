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

/*! \file  Global.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
			  

#ifndef _H_Global
#define _H_Global

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>

#include "def_vars.h"

//*********************************************************************************************

class Global
{
public:

	YAML::Node para;            //stores all the para in the para.yaml file

	map<string, map<string, map<string, string> > > basis_table;
	map<string, int> no_components_table;

	map<string, int> global_data_packet_size_table;
	map<string, map<string, int> > spectral_probe_packet_size_table;
	
	map<string, int> real_probe_packet_size_table;
	

	struct program {
		const string version;
		
		string kind;
		string iter_or_diag;
		string alias_option;
		string integration_scheme;
		string basis_type;
		string decomposition;				// PENCIL, SLAB
		bool LES_switch;
		bool T_exists;
		bool W_exists;
		bool apply_strong_realitycond_alltime_switch;
		bool apply_weak_realitycond_alltime_switch;
		bool low_dimensional_switch;
		bool two_and_half_dimension;
		bool two_dimension;
		bool helicity_switch;

		int dt_option;  // 0=fixed; 1=dx/umax; 
		string sincostr_switch;             // SINX, COSX, SCC, CSC, CCS, SSC, CSS, SCS
		string sincostr_switch_Vx, sincostr_switch_Vy, sincostr_switch_Vz;
		string sincostr_switch_F;
		string sincostr_switch_Visqr;
		string sincostr_switch_VxVy, sincostr_switch_VxVz, sincostr_switch_VyVz;
		string sincostr_switch_FVx, sincostr_switch_FVy, sincostr_switch_FVz;
		string sincostr_switch_divergence;	
		
		

		program();
		void Print(int my_level);
	} program;

	struct PHYSICS {
		string	Pr_option;					// Prandtl number switch (PRLARGNE..) for RB
		string	Uscaling;					// UBscaling (ULARGE... ) for RB
		DP		Rayleigh;							// Rayleigh number
		DP		Prandtl;							// Prandtl number
		int		temperature_grad;			// +1 for convection; -1 for stratification;
		
		DP Chandrasekhar; // Chandrashekhar number
		DP Prandtl_mag; //Magnetic Prandlt number
		DP Prandtl_c;
		DP Reynolds;
		DP Reynolds_mag;
		DP Peclet;
		DP Peclet_c;
		// factor for u3 in temperature eqn
		void Print(int my_level);
	} PHYSICS;

	struct MRBC {
		string	Pr_option;					// Prandtl number switch (PRLARGNE..) for RB
		string	Uscaling;					// UBscaling (ULARGE... ) for RB
		DP		Pr;							// Prandtl number
		
		DP		RaD;
		DP		RaM;
		DP		SSD;						// Surface saturation deficit
		DP		CSA;						// Condensate in saturated ascent
		void Print(int my_level);
	} MRBC;
	
	struct field {
		//string SPACE;	// REAL or SPECTRAL
		bool incompressible;
		bool waveno_switch;
		int anisotropy_dirn;  //1,2,3
		
		int N[4];
		int Nx, Ny, Nz;
		int howmany;  // how many fields, e.g. RBC has 4.

		vector<DP> kfactor;
		DP xfactor[4];
		vector<DP> L;		// Box size
		DP Delta_x[4];	// Delta x along each direction
		vector<DP> diss_coefficients;
		vector<DP> hyper_diss_coefficients;
		vector<int> hyper_diss_exponents;
		
		TinyVector<int,3> shape_complex_array;
		TinyVector<int,3> shape_real_array;
		
		DP twobyL1;			// for Chebyshev basis
		
		int maxlx, maxly, maxlz;
		
		void Print(int my_level);
	} field;

	struct mpi {
		int my_id;
		int numprocs;
		int master_id;
		bool master;
		
			// for pencil
		int num_p_hor;              // processors along the horizontal direction  
		int num_p_vert;             // processors along the vertical direction
		
		int num_x_procs;
		int num_y_procs;
		int num_z_procs;

		int num_x_procs_real;
		int num_y_procs_real;
		int num_z_procs_real;

		int my_vert_pcoord;         // processor coordinate along vertical direction
		int my_hor_pcoord;          // processor coordinate along horizontal direction

		int my_x_pcoord;
		int my_y_pcoord;
		int my_z_pcoord;

		int my_x_pcoord_real;
		int my_y_pcoord_real;
		int my_z_pcoord_real;
		
		MPI_Comm MPI_COMM_HOR_SEGMENT;
		MPI_Comm MPI_COMM_VERT_SEGMENT;
		
		mpi();
		void Print(int my_level);
	} mpi;


	struct time{
		DP init; 
		DP final;
		DP dt_fixed;
		DP Courant_no;
		string job_time;
		
		DP dt;				// variable time including CFL
		DP now;
		DP previous;
		DP keplerian;
		bool dt_computation_done;
		
		clock_t job_time_final;
		void Print(int my_level);
	} time;

	struct force{
		bool U_switch;
		bool W_switch;
		bool T_switch;
		bool C_switch;
		bool configuration_done;

		int field_procedure;

		Array<int,1> int_para;
		Array<DP,1> double_para;
		Array<string,1> string_para;

		struct modes{
			bool read_done;
			int number;
			int number_components;  // 2 for fluid, 4 for U,W (incompressiblity)
			Array<int,2> coords;
			Array<DP,2> field_array_real;
			Array<complx,2> field_array_complex;
			void Print(int my_level);
		} modes;
		void Print(int my_level);
	} force;

	struct io {
		string data_dir;
		int input_field_procedure;
		bool input_vx_vy_switch;
		bool output_vx_vy_switch;
		const int output_precision;  // 6 for float and 12 for double in ascii format

		bool output_real_field_done;
		bool output_field_k_done;
		bool output_pressure_done;
		bool output_pressure_spectrum_done;
		bool real_space_field_available;
		
		bool output_nlin_magnitude_done;
		
		vector<int>	diagnostic_procedures;

		vector<int> N_in_reduced;
		vector<int> N_out_reduced;

		Array<int,1> int_para;
		Array<DP,1> double_para;
		Array<string,1> string_para;
		
		struct init_cond_modes {
			int number;
			int number_components;		// if incompressible: subtract 1 for V, W
			Array<int,2> coords;
			Array<DP,2> field_array_real;
			Array<complx,2> field_array_complex;
			void Print(int my_level);
			int In_table(int kx, int ky, int kz);
		} init_cond_modes;
		
		struct global_data {
			unsigned long buffer_size;
			unsigned long buffer_index;		// initialize to zero
			unsigned long packet_size;

			Array<DP,1> buffer;
		} global_data;

		struct probes{
			struct spectral_space{
				int number;
				unsigned long buffer_size;
				Array<int,2> coords;
				unsigned long buffer_index;		// initialize to zero
				int packet_size;		// for each time frame (define)
				// (global.io.probes.spectral_space.number * 6) +1 + Tk for each field; 

				Array<DP,1> field_buffer;  // in all nodes
				void Print(int my_level);
			} spectral_space;

			struct real_space {
				int number;
				unsigned long buffer_size;
				Array<int,2> coords;
				unsigned long buffer_index;		// initialize to zero
				int packet_size;

				Array<DP,1> field_buffer;
				void Print(int my_level);
			} real_space;

			void Print(int my_level);
		} probes;

		struct time{
			DP global_save_next;
			DP complex_field_save_next;
			DP field_frequent_save_next;
			DP field_reduced_save_next;
			DP real_field_save_next;
			DP field_k_save_next;
			DP field_r_save_next;
			DP spectrum_save_next;
			DP pressure_save_next;
			DP pressure_spectrum_save_next;
			DP flux_save_next;
			DP shell_to_shell_save_next;
			DP ring_spectrum_save_next;
			DP ring_to_ring_save_next;
			DP cylindrical_ring_spectrum_save_next;
			DP cylindrical_ring_to_ring_save_next;
			DP structure_fn_save_next;
			DP Tk_shell_spectrum_save_next;
			DP Tk_ring_spectrum_save_next;
			DP Tk_cylindrical_ring_spectrum_save_next;
			DP cout_save_next;


			DP	global_save_interval;	
			DP	complex_field_save_interval; 
			DP	field_frequent_save_interval;
			DP	field_reduced_save_interval;
			DP	real_field_save_interval;
			DP	field_k_save_interval;				
			DP	field_r_save_interval;	
			DP  pressure_save_interval;
			DP	spectrum_save_interval;
			DP	pressure_spectrum_save_interval;
			DP	flux_save_interval;
			DP	shell_to_shell_save_interval;
			DP	ring_spectrum_save_interval;
			DP	ring_to_ring_save_interval;
			DP	cylindrical_ring_spectrum_save_interval;
			DP	cylindrical_ring_to_ring_save_interval;
			DP	structure_fn_save_interval;
			DP  Tk_shell_spectrum_save_interval;
			DP  Tk_ring_spectrum_save_interval;
			DP  Tk_cylindrical_ring_spectrum_save_interval;
			DP	cout_save_interval;
			
			bool	global_save_last;	
			bool	complex_field_save_last; 
			bool	field_frequent_save_last;
			bool	field_reduced_save_last;
			bool	real_field_save_last;
			bool	field_k_save_last;				
			bool	field_r_save_last;	
			bool pressure_save_last;
			bool	spectrum_save_last;
			bool	pressure_spectrum_save_last;
			bool	flux_save_last;
			bool	shell_to_shell_save_last;
			bool	ring_spectrum_save_last;
			bool	ring_to_ring_save_last;
			bool	cylindrical_ring_spectrum_save_last;
			bool	cylindrical_ring_to_ring_save_last;
			bool	structure_fn_save_last;
			bool Tk_shell_spectrum_save_last;
			bool Tk_ring_spectrum_save_last;
			bool Tk_cylindrical_ring_spectrum_save_last;
			bool	cout_save_last;
			
			void Print(int my_level);
		} time;

		io();
		void Print(int my_level);
		void Dump_buffer(ofstream& file, Array<DP,1> buffer, unsigned long& buffer_index, unsigned long packet_size);
	} io;


	struct  spectrum {
		struct shell {
			bool turnon;
			int no_shells;
			void Print(int my_level);
		} shell;

		struct ring {
			bool turnon;
			int no_shells;     // prog defined
			int no_sectors;   // to be read
			string sector_option;
			Array<DP,1> sector_angles;
			void Print(int my_level);
		} ring;

		struct cylindrical_ring {
			bool turnon;
			int no_shells;		// prog define
			int no_slabs;		 // to be read
			string kpll_option;
			Array<DP,1> kpll_array;  // z coords for the sections
			void Print(int my_level);
		} cylindrical_ring;

		void Print(int my_level);
	} spectrum;

	struct energy_transfer {
		bool turnon;
		bool helicity_flux_switch;
		bool helicity_shell_to_shell_switch;
		bool Elsasser;
		bool Vpll_switch;  // for energy transfer of parallel diagnostics

		struct flux {
			bool turnon;
			int no_spheres;
			Array<DP,1> radii;
			void Print(int my_level);
		} flux;

		struct shell_to_shell {
			bool turnon;
			int no_shells;
			Array<DP,1> radii;
			void Print(int my_level);
		} shell_to_shell;

		struct ring_to_ring {
			bool turnon;
			int no_shells;
			int no_sectors;
			string sector_option;
			Array<DP,1> radii;
			Array<DP,1> sector_angles;
			void Print(int my_level);
		} ring_to_ring;

		struct cylindrical_ring_to_ring {
			bool turnon;
			int no_shells;
			int no_slabs;
			Array<DP,1> radii;
			string kpll_option;
			Array<DP,1> kpll_array;  // kpll coords for the sections
			void Print(int my_level);
		} cylindrical_ring_to_ring;
		void Print(int my_level);
	} energy_transfer;

	struct structure_fn {
		bool turnon;
		bool box_switch;
		bool planar_switch;
		bool approx_switch;

		int qmin;
		int qmax;
		int array_size;
		int rmax_div_st_fn_r_max;
		int nearest_neighbour_index;  // points_to_skip+1
		int dist_farthest_proc;			// for i1 index...
	} structure_fn;


	struct temp_array {
		Array<complx,3> X;	

		// Used in div calc & in MHD (Compute_nlin_offdiag(W))
		Array<complx,3> X2;					

		//!  \f$ (local_{N2}, N_1, N_3/2+1) \f$. 	
		Array<DP,3> Xr;

		// Ar not reqd NOW.. CHANGE to Xr2
		// BOTH Ar and Xr are reqd. Xr is the temp array..Ar is the real space array.
		//!  temp array \f$ (local_{Ny}, N_x, N_z/2+1) \f$.
		// For transforms.. MUST
		// Also in Compute_nlin_offdiag(W)
		Array<DP,3> Xr2;


		// for the function satisfy reality condition
		Array<complx,2> plane_xy;  // Ny,Nx 	
		Array<complx,2> plane_xy_inproc;  // local_Nx,Ny orlocal_Nx_hor,Ny 
		
		// for Helmholtz solver: Chebyshev basis
		Array<complx,1> in_helm_complex;    // Nx
		Array<complx,1> out_helm_complex;    // Nx
   
		Array<complx,1> V1_x;
		
		// pl remove to
		Array<DP,1> in_helm_real;    // Nx
		Array<DP,1> out_helm_real;    // Nx

		
		Array<DP,3> pressure_plus,pressure_minus;   // local_Nz,Ny,Nx
		Array<DP,1> pressure_derivative_real_1d;
		
		Array<DP,3> vx_plus, vx_minus;      // local_Nz,Ny, Nx
		Array<DP,3> Xreal;
		
		Array<DP,4> influence_matrix;       // local_Nz,Ny,2,2
		//ChFF end
		
	} temp_array;


	struct myconstant{
		const complx  I;

		/// minusI = -sqrt(-1).			
		const complx  minusI;

		/// minus2I = -2*sqrt(-1).			
		const complx  minus2I;	

		// Number of digits for output files
		// NOTE:  for double only.. For float put MY_PRECISION = 6
		//const int MY_PRECISION = 12;


		const DP  MYEPS;
		const DP  MYEPS2;

		const int MY_MAX_INT;
		// cut off while reading diagnostic_procedure() array from input file and similar ops

		/// Infinite radius.. All the modes outside -- for flux and shelltr calc
		const DP INF_RADIUS;

		const DP INF_TIME; 
		
		const int MAX_NO_GLOB_BUFFER_PACKETS;
		const int MAX_NO_PROBE_PACKETS;

		myconstant();
	} myconstant;

	Global();

	void Global_Parse(int argc, char** argv, bool is_test_module=false);
	void Global_Read();
	void Global_init_default();
	void Process_global_vars_basic();
	void Process_global_vars_advanced();
	void Assign_parameters();
	void Alias_some_global_vars();
	bool Input_provided(const YAML::Node& node, const string parameter);
	
	template <typename T> 
	void Assign_if_input_provided(const YAML::Node& node, const string parameter, T &var, T default_value);
	
	void Show_error(string message);

	void Print();
};

#endif

//******************************* End of Global.h ********************************
