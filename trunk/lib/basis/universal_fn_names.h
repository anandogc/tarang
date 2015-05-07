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

	void Last_component(int kx, int ky, int kz, Real &Vx, Real &Vy, Real &Vz);
	void Last_component(int kx, int ky, int kz, Complex& Vx, Complex& Vy, Complex& Vz);
	void Dealias(Array<Complex,3> A);
	bool Is_dealiasing_necessary(Array<Complex,3> A, Real outer_radius);
    void Satisfy_strong_reality_condition_in_Array(Array<Complex,3> A);
    void Satisfy_weak_reality_condition_in_Array(Array<Complex,3> A);
    void Test_reality_condition_in_Array(Array<Complex,3> A);

    // transform
    
	void Forward_transform(Array<Real,3> Ar, Array<Complex,3> A);
    
	void Inverse_transform(Array<Complex,3> A, Array<Real,3> Ar);
    
	void Xderiv(Array<Complex,3> A, Array<Complex,3> B);
	void Yderiv(Array<Complex,3> A, Array<Complex,3> B);
	void Zderiv(Array<Complex,3> A, Array<Complex,3> B);

	void Add_Xderiv(Array<Complex,3> A, Array<Complex,3> B);
	void Add_Yderiv(Array<Complex,3> A, Array<Complex,3> B);
	void Add_Zderiv(Array<Complex,3> A, Array<Complex,3> B);

	void  Xderiv(Array<Real,3> A, Array<Real,3> B);

	void Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B);
	void Subtract_Laplacian(Real factor, Array<Complex,3> A, Array<Complex,3> B);

	// Energy
    
    /*Real Get_local_energy_XZ_plane(Array<Complex,3> A, int ny);
    
    Real Get_local_energy_XZ_plane
	(
		Array<Complex,3> A, Array<Complex,3> B, 
		int ny
	);*/
    
	Real Get_local_energy_real_space(Array<Real,3> Ar);
	Real Get_local_energy(Array<Complex,3> A);
	Real Get_local_energy_real_space(Array<Real,3> Ar, Array<Real,3> Br);
    Real Get_local_energy(Array<Complex,3> A, Array<Complex,3> B);

    void Compute_local_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Real &local_helicity1, Real &local_helicity2, 
        Real &local_dissipation_H1, Real &local_dissipation_H2
    );
    
    void Compute_total_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Real &total_helicity1, Real &total_helicity2, 
        Real &total_dissipation_H1, Real &total_dissipation_H2
     );
    
   void Compute_local_shell_spectrum_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,1> local_H1k1, Array<Real,1> local_H1k2, Array<Real,1> local_H1k3, 
        Array<Real,1> local_H1k_count
    );
    
    void Compute_shell_spectrum_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,1> H1k1, Array<Real,1> H1k2, Array<Real,1> H1k3 
    );

    void Compute_local_ring_spectrum_helicity
    ( 
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,2> local_H1k1, Array<Real,2> local_H1k2, Array<Real,2> local_H1k3
    );
    
    void Compute_ring_spectrum_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,2> local_H1k1, Array<Real,2> local_H1k2, Array<Real,2> local_H1k3
    );

    void Compute_local_cylindrical_ring_spectrum_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,2> local_H1k1, Array<Real,2> local_H1k2
    );
    
    void Compute_cylindrical_ring_spectrum_helicity
    (
        Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
        Array<Real,2> local_H1k1, Array<Real,2> local_H1k2
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
	
	Complex Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A);
	TinyVector<Complex,3> Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
//	Real Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A);
//	TinyVector<Real,3> Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	
	void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Complex field);
	void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V);
	void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Real field);
	void Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V);

    void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Complex field);
    void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V);
    void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Real field);
    void Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V);
	
	Complex Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A);
	TinyVector<Complex,3> Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
//	Real Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A);
//	TinyVector<Real,3> Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Complex field);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Real field);
	void Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V);
	
	void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Complex field);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Real field);
	void Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V);
	
	int Get_lx_real_space(int rx);
	int Get_ly_real_space(int ry);
	int Get_lz_real_space(int rz);

	int Get_rx_real_space(int lx);
	int Get_ry_real_space(int ly);	
	int Get_rz_real_space(int lz);


	bool Probe_in_me_real_space(int rx, int ry, int rz);
	Real Get_real_field(int rx, int ry, int rz, Array<Real,3> A);
	TinyVector<Real,3> Get_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az);
	void Assign_real_field(int rx, int ry, int rz, Array<Real,3> A, Real field);
	void Assign_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az, TinyVector<Real,3> V);
	
	
	void Wavenumber(int lx, int ly, int lz, TinyVector<Real,3> & K);
	void Wavenumber(int lx, int ly, int lz, TinyVector<Complex,3> & K);
	
	Real Kmagnitude(int lx, int ly, int lz);
	
	int Min_radius_outside();
	int Max_radius_inside();

	Real Approx_number_modes_in_shell(int radius);
	
	Real Multiplicity_factor(int lx, int ly, int lz);
	
	Real Modal_energy(int lx, int ly, int lz, Array<Complex,3> A);
	Real Get_Modal_helicity
	(
		int lx, int ly, int lz,
		Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az
	);
	
	void Compute_Modal_vorticity
	(
		 int lx, int ly, int lz, 
		 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		 TinyVector<Complex,3> &vorticity
	);
	
	void Compute_Modal_vorticity_y_component
	(
		 int lx, int ly, int lz, 
		 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
		 Complex &vort_y
	 );
	
	Real AnisKpll(int lx, int ly, int lz);
	Real AnisKperp(int lx, int ly, int lz);
	Real AnisKh1(int lx, int ly, int lz);
	Real AnisKh2(int lx, int ly, int lz);
	Real Anis_min_Kpll();
	Real Anis_max_Kpll() ;
	
	int Anis_max_Krho_radius_inside();
	Real Get_max_polar_angle() ;
	
	Real AnisKvect_polar_angle(int lx, int ly, int lz);
	Real AnisKvect_azimuthal_angle(int lx, int ly, int lz);

    int Read(Array<Complex,3> A, BasicIO::H5_plan plan, string file_name, string dataset_name="");
    int Read(Array<Real,3> Ar, BasicIO::H5_plan plan, string file_name, string dataset_name="");

    int Write(Array<Complex,3> A, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name="");
    int Write(Array<Real,3> Ar, BasicIO::H5_plan plan, string folder_name, string file_name, string dataset_name="");

