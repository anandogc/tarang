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
 * Tarang-2 is distributed in the hope that it will be useful,Get_
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file  universea_fn_names.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

	// KEEP FNS NAMES THAT ARE IN "ALL" BASIS DIRS
	// basic

	void Last_component(int kx, int ky, int kz, DP &Vx, DP &Vy, DP &Vz);
	void Last_component(int kx, int ky, int kz, complx& Vx, complx& Vy, complx& Vz);
	void Dealias(Array<complx,3> A);
	bool Is_dealiasing_necessary(Array<complx,3> A, DP outer_radius);
    void Satisfy_strong_reality_condition_in_Array(Array<complx,3> A);
    void Satisfy_weak_reality_condition_in_Array(Array<complx,3> A);
    void Test_reality_condition_in_Array(Array<complx,3> A);

    // transform
    
	void Forward_transform(Array<DP,3> Ar, Array<complx,3> A);
    
	void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar);
    
	void Xderiv(Array<complx,3> A, Array<complx,3> B);
	void Yderiv(Array<complx,3> A, Array<complx,3> B);
	void Zderiv(Array<complx,3> A, Array<complx,3> B);

	void Add_Xderiv(Array<complx,3> A, Array<complx,3> B);
	void Add_Yderiv(Array<complx,3> A, Array<complx,3> B);
	void Add_Zderiv(Array<complx,3> A, Array<complx,3> B);

	void  Xderiv(Array<DP,3> A, Array<DP,3> B);

	void Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B);
	void Subtract_Laplacian(DP factor, Array<complx,3> A, Array<complx,3> B);

	// Energy
    
    /*DP Get_local_energy_XZ_plane(Array<complx,3> A, int ny);
    
    DP Get_local_energy_XZ_plane
	(
		Array<complx,3> A, Array<complx,3> B, 
		int ny
	);*/
    
	DP Get_local_energy_real_space(Array<DP,3> Ar);
	DP Get_local_energy(Array<complx,3> A);
	DP Get_local_energy_real_space(Array<DP,3> Ar, Array<DP,3> Br);
    DP Get_local_energy(Array<complx,3> A, Array<complx,3> B);

    void Compute_local_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        DP &local_helicity1, DP &local_helicity2, 
        DP &local_dissipation_H1, DP &local_dissipation_H2
    );
    
    void Compute_total_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        DP &total_helicity1, DP &total_helicity2, 
        DP &total_dissipation_H1, DP &total_dissipation_H2
     );
    
   void Compute_local_shell_spectrum_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
        Array<DP,1> local_H1k_count
    );
    
    void Compute_shell_spectrum_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3 
    );

    void Compute_local_ring_spectrum_helicity
    ( 
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
    );
    
    void Compute_ring_spectrum_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,2> local_H1k1, Array<DP,2> local_H1k2, Array<DP,2> local_H1k3
    );

    void Compute_local_cylindrical_ring_spectrum_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,2> local_H1k1, Array<DP,2> local_H1k2
    );
    
    void Compute_cylindrical_ring_spectrum_helicity
    (
        Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
        Array<DP,2> local_H1k1, Array<DP,2> local_H1k2
    );
    
	
		// Inline functions							 
	int Get_kx(int lx);
	int Get_lx(int kx);
	int Get_ix(int kx);
	
	int Get_ky(int ly);
	int Get_ly(int ky);
	int Get_iy(int ky);

	int Get_kz(int lz);
	int Get_lz(int kz);
	int Get_iz(int kz);
	
	bool Probe_in_me(int kx, int ky, int kz);
	
	complx Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A);
	TinyVector<complx,3> Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
//	DP Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A);
//	TinyVector<DP,3> Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
	
	void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A,complx field);
	void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V);
	void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A, DP field);
	void Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V);

    void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A,complx field);
    void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V);
    void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A, DP field);
    void Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V);
	
	complx Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A);
	TinyVector<complx,3> Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
//	DP Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A);
//	TinyVector<DP,3> Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az);
	
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A,complx field);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V);
	
	void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A,complx field);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V);
	
	int Get_lx_real_space(int rx);
	int Get_ly_real_space(int ry);
	int Get_lz_real_space(int rz);

	int Get_rx_real_space(int lx);
	int Get_ry_real_space(int ly);	
	int Get_rz_real_space(int lz);


	bool Probe_in_me_real_space(int rx, int ry, int rz);
	DP Get_real_field(int rx, int ry, int rz, Array<DP,3> A);
	TinyVector<DP,3> Get_real_field(int rx, int ry, int rz, Array<DP,3> Ax, Array<DP,3> Ay, Array<DP,3> Az);
	void Assign_real_field(int rx, int ry, int rz, Array<DP,3> A, DP field);
	void Assign_real_field(int rx, int ry, int rz, Array<DP,3> Ax, Array<DP,3> Ay, Array<DP,3> Az, TinyVector<DP,3> V);
	
	
	void Wavenumber(int lx, int ly, int lz, TinyVector<DP,3> & K);
	void Wavenumber(int lx, int ly, int lz, TinyVector<complx,3> & K);
	
	DP Kmagnitude(int lx, int ly, int lz);
	
	int Min_radius_outside();
	int Max_radius_inside();

	DP Approx_number_modes_in_shell(int radius);
	
	DP Multiplicity_factor(int lx, int ly, int lz);
	
	DP Modal_energy(int lx, int ly, int lz, Array<complx,3> A);
	DP Get_Modal_helicity
	(
		int lx, int ly, int lz,
		Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
	);
	
	void Compute_Modal_vorticity
	(
		 int lx, int ly, int lz, 
		 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
		 TinyVector<complx,3> &vorticity
	);
	
	void Compute_Modal_vorticity_y_component
	(
		 int lx, int ly, int lz, 
		 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
		 complx &vort_y
	 );
	
	DP AnisKpll(int lx, int ly, int lz);
	DP AnisKperp(int lx, int ly, int lz);
	DP AnisKh1(int lx, int ly, int lz);
	DP AnisKh2(int lx, int ly, int lz);
	DP Anis_min_Kpll();
	DP Anis_max_Kpll() ;
	
	int Anis_max_Krho_radius_inside();
	DP Get_max_polar_angle() ;
	
	DP AnisKvect_polar_angle(int lx, int ly, int lz);
	DP AnisKvect_azimuthal_angle(int lx, int ly, int lz);

    int Read(Array<complx,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name="");
    int Read(Array<DP,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name="");

    int Write(Array<complx,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name="");
    int Write(Array<DP,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name="");
    