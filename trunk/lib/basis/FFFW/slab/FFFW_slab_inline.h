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


/*! \file four_inline.h 
 * 
 * @brief Inline functions to compute array indices \f$ \vec{i} \f$ given wavenumber  
 * \f$ \vec{k} \f$ and viceversa.
 * Also contains other useful functions like Kmagnitude, Max radius etc.
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 * \f$ \vec{i} \f$.
 * Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber \f$ \vec{k} \f$ 
 *		using \f$ K_i = k_i * f_i \f$
 * where \f$ f_i \f$ is the kfactor[i]. 
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 *	or gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.   
 *	Actual wavenumber is preferred.
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	August 2008
 * @bug		No known bugs
 *
 */ 
 

#include "FFFW_slab.h"

using namespace blitz;


//*********************************************************************************************

/*! @brief	Get grid waveno kx given first local array index lx.
 * 
 *  if i1 <=Ni/2, ki=i1; else kx=i1-Ni.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param lx  first local index of an array
 * \return kx corresponding to lx
 */
inline int FFFW_SLAB::Get_kx(int lx)  
{ 			
    return  ((local_Nx_start + lx) <= N[1]/2) ? 
            (local_Nx_start + lx) : (local_Nx_start + lx-N[1]); 
} 
 

/*! @brief	Get local array index lx given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param	kx  grid wavenumber along x
 * \return	lx  local array index along x
 */
inline int FFFW_SLAB::Get_lx(int kx)  
{ return  (kx >= 0) ? (kx-local_Nx_start) : (kx + N[1]-local_Nx_start); }	


/*! @brief	Get array index ix given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 * 
 * \param	kx  grid wavenumber along x
 * \return	i1  array index along x
 */
inline int FFFW_SLAB::Get_ix(int kx)			
{ 	return  (kx >= 0) ? kx : (kx + N[1]);  }


/*! @brief	Get grid waveno ky given first local array index ly.
 * 
 *  If ly<=N2/2, k2=ly; else ky=ly-N2.
 * 
 * \param ly  second local index of an array
 * \return kx corresponding to lx
 */
 inline int FFFW_SLAB::Get_ky(int ly)  { return  (ly <= N[2]/2) ? ly : (ly-N[2]); } 

/*! @brief	Get local array index ly given grid waveno ky.
 * 
 *  If ky>=0, ly=ky; else ly=ky+N2.
 * 
 * \param	ky  grid wavenumber along y
 * \return	ly  local array index along y
 */
inline int FFFW_SLAB::Get_ly(int ky) { return  (ky >= 0) ? ky : (ky + N[2]); } 

inline int FFFW_SLAB::Get_iy(int ky) { return  0; } //This is not implemented for this basis

inline int FFFW_SLAB::Get_kz(int lz)  { return lz; } 

inline int FFFW_SLAB::Get_lz(int kz)  { return kz; }

	// array index
inline int FFFW_SLAB::Get_iz(int kz)  { return kz; }



inline bool FFFW_SLAB::Probe_in_me(int kx, int ky, int kz) 
{
	int lx = Get_lx(kx);
	return ((lx >= 0) && (lx < local_Nx));
}

inline void FFFW_SLAB::Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A, complx &field)
{ 
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if  ((lx >= 0) && (lx < local_Nx))
		field = A(lx, ly, kz); 
}

inline void FFFW_SLAB::Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> &V)
{
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if  ((lx >= 0) && (lx < local_Nx)) 
		V =  Ax(lx, ly, kz), Ay(lx, ly, kz), Az(lx, ly, kz);

}


inline void FFFW_SLAB::Get_spectral_field(int kx, int ky, int kz, Array<complx,3> A,DP &field)
{ 
	cout << "MYERROR: FFFW_SLAB::Get_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Get_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> &V)
{
	
	cout << "MYERROR: FFFW_SLAB::Get_spectral_field(); Use complex data type " << endl;
}

	// Assign....

inline void FFFW_SLAB::Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A, complx field)
{ 
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if ((lx >= 0) && (lx < local_Nx))
		A(lx, ly, kz) = field;
}

inline void FFFW_SLAB::Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V)
{
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if ((lx >= 0) && (lx < local_Nx)) {
		Ax(lx, ly, kz) = V(0);
		Ay(lx, ly, kz) = V(1);
		Az(lx, ly, kz) = V(2);
	}
}

inline void FFFW_SLAB::Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> A,DP field)
{ 
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Assign_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V)
{
	
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl;
}


inline void FFFW_SLAB::Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A, complx field)
{ 
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if ((lx >= 0) && (lx < local_Nx))
		A(lx, ly, kz) += field;
}

inline void FFFW_SLAB::Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V)
{
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
	
	if ((lx >= 0) && (lx < local_Nx)) {
		Ax(lx, ly, kz) += V(0);
		Ay(lx, ly, kz) += V(1);
		Az(lx, ly, kz) += V(2);
	}
}

inline void FFFW_SLAB::Add_spectral_field(int kx, int ky, int kz, Array<complx,3> A,DP field)
{ 
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Add_spectral_field(int kx, int ky, int kz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V)
{
	
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl;
}

	// Local get

inline void FFFW_SLAB::Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, complx &field)
{ 
	field = A(lx, ly, lz); 
}

inline void FFFW_SLAB::Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> &V)
{
	
	V =  Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
}


inline void FFFW_SLAB::Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A,DP &field)
{ 
	cerr << "MYERROR: FFFW_SLAB::Get_local_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Get_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> &V)
{
	
	 cerr << "MYERROR: FFFW_SLAB::Get_local_spectral_field(); Use complex data type " << endl;
}


inline void FFFW_SLAB::Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, complx field)
{ 	
	if ((lx >= 0) && (lx < local_Nx))
		A(lx, ly, lz) = field;
}


inline void FFFW_SLAB::Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V)
{	

	if ((lx >= 0) && (lx < local_Nx)) {
		Ax(lx, ly, lz) = V(0);
		Ay(lx, ly, lz) = V(1);
		Az(lx, ly, lz) = V(2);
	}
}

inline void FFFW_SLAB::Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V)
{
	
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl;
}

inline void FFFW_SLAB::Assign_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field)
{ 
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, complx field)
{ 	
	if ((lx >= 0) && (lx < local_Nx))
		A(lx, ly, lz) += field;
}

inline void FFFW_SLAB::Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<complx,3> V)
{	

	if ((lx >= 0) && (lx < local_Nx)) {
		Ax(lx, ly, lz) += V(0);
		Ay(lx, ly, lz) += V(1);
		Az(lx, ly, lz) += V(2);
	}
}

inline void FFFW_SLAB::Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> A, DP field)
{ 
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void FFFW_SLAB::Add_local_spectral_field(int lx, int ly, int lz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V)
{
	
	cout << "MYERROR: FFFW_SLAB::Assign_spectral_field(); Use complex data type " << endl;
}

				 
// REAL - SPACE FNS
inline int FFFW_SLAB::Get_lx_real_space(int rx)  
{ return  global.fft.transpose ? rx : (rx-local_Nx_start); }	

inline int FFFW_SLAB::Get_ly_real_space(int ry)  
{ return  global.fft.transpose ? (ry-local_Nx_start) : ry; }	

inline void FFFW_SLAB::Get_lz_real_space(int rz, int &lz, bool &is_real_part)  
{	lz = rz/2;  
	is_real_part = rz % 2;
}

inline int FFFW_SLAB::Get_rx_real_space(int lx)  
{ return  global.fft.transpose ? lx : (lx+local_Nx_start); }

inline int FFFW_SLAB::Get_ry_real_space(int ly)  
{ return  global.fft.transpose ? (ly+local_Ny_start) : ly; }	

inline void FFFW_SLAB::Get_rz_real_space(int lz, int& rz, bool is_real_part)  
{	
if (is_real_part)
	rz = 2*lz;  
else
	rz = 2*lz+1;
}	

inline bool FFFW_SLAB::Probe_in_me_real_space(int rx, int ry, int rz) 
{
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
	
	if (global.fft.transpose)
		return ((ly >= 0) && (ly < local_Ny));
	else
		return ((lx >= 0) && (lx < local_Nx));
}


inline void FFFW_SLAB::Get_real_field(int rx, int ry, int rz, Array<complx,3> A, DP &field)
{	

	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
    
    int lz;
	bool is_real_part;

    Get_lz_real_space(rz, lz, is_real_part);
	
	
	if (global.fft.transpose) {
		if ((ly >= 0) && (ly < local_Ny))  {
			if (is_real_part)
				field = real(A(ly, lx, lz));
			else
				field = imag(A(ly, lx, lz));
		}	
	}
	
	else {
		if ((lx >= 0) && (lx < local_Nx))  {
			if (is_real_part)
				field = real(A(lx, ly, lz));
			else
				field = imag(A(lx, ly, lz));
		}
	}
}

inline void FFFW_SLAB::Get_real_field(int rx, int ry, int rz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> &V)
{
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
    
    int lz;  bool is_real_part;
    Get_lz_real_space(rz, lz, is_real_part);
	
	if (global.fft.transpose) {
		if ((ly >= 0) && (ly < local_Ny))  {
			if (is_real_part)
				V =  real(Ax(ly, lx, lz)), real(Ay(ly, lx, lz)), real(Az(ly, lx, lz));
			else
				V =  imag(Ax(ly, lx, lz)), imag(Ay(ly, lx, lz)), imag(Az(ly, lx, lz));
		}	
	}
	
	else {
		if ((lx >= 0) && (lx < local_Nx))  {
			if (is_real_part)
				V =  real(Ax(lx, ly, lz)), real(Ay(lx, ly, lz)), real(Az(lx, ly, lz));
			else
				V =  imag(Ax(lx, ly, lz)), imag(Ay(lx, ly, lz)), imag(Az(lx, ly, lz));
		}
	}
}



inline void FFFW_SLAB::Assign_real_field(int rx, int ry, int rz, Array<complx,3> A, DP field)
{	
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
	
	int lz;  bool is_real_part;
    Get_lz_real_space(rz, lz, is_real_part);
	
	
	if (global.fft.transpose) {
		if ((ly >= 0) && (ly < local_Ny))  {
			if (is_real_part)
				real(A(ly, lx, lz)) = field;
			else
				imag(A(ly, lx, lz)) = field;
		}	
	}
	
	else {
		if ((lx >= 0) && (lx < local_Nx))  {
			if (is_real_part)
				real(A(lx, ly, lz)) = field;
			else
				imag(A(lx, ly, lz)) = field;
		}
	}
}

inline void FFFW_SLAB::Assign_real_field(int rx, int ry, int rz, Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, TinyVector<DP,3> V)
{
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
    
    int lz;  bool is_real_part;
    Get_lz_real_space(rz, lz, is_real_part);
	
	if (global.fft.transpose) {
		if ((ly >= 0) && (ly < local_Ny))  {
			if (is_real_part) {
				real(Ax(ly, lx, lz)) = V(0);
				real(Ay(ly, lx, lz)) = V(1);
				real(Az(ly, lx, lz)) = V(2);
			}
			else {
				imag(Ax(ly, lx, lz)) = V(0);
				imag(Ay(ly, lx, lz)) = V(1);
				imag(Az(ly, lx, lz)) = V(2);
			}
		}	
	}
	
	else {
		if ((lx >= 0) && (lx < local_Nx))  {
			if (is_real_part) {
				real(Ax(lx, ly, lz)) = V(0);
				real(Ay(lx, ly, lz)) = V(1);
				real(Az(lx, ly, lz)) = V(2);
			}
			else {
				imag(Ax(lx, ly, lz)) = V(0);
				imag(Ay(lx, ly, lz)) = V(1);
				imag(Az(lx, ly, lz)) = V(2);
			}
		}
	}
}


/**********************************************************************************************

	Compute Wavenumber components for grid index (i1,i2).

***********************************************************************************************/

// Real K
inline void FFFW_SLAB::Wavenumber(int lx, int ly, int lz, TinyVector<DP,3> & K)
{
	 K = Get_kx(lx)*kfactor[1],  Get_ky(ly)*kfactor[2], lz*kfactor[3];
}


// Complex K; The imaginary part is zero.  Written to use cross function of blitz.
// Omega = cross(V,K).
inline void FFFW_SLAB::Wavenumber(int lx, int ly, int lz, TinyVector<complx,3> & K)
{
	 K = complx(Get_kx(lx)*kfactor[1], 0.0), complx(Get_ky(ly)*kfactor[2], 0.0), complx(lz*kfactor[3], 0.0);
}


/**********************************************************************************************

				 If wavenos computed using actual wavenumber:  
						global.field.waveno_switch = 0;  
						Ki = Kfactor[i]*grid[i]
				 
				 If wavenos computed using grid wavenumber:  Ki = grid[i]
						global.field.waveno_switch = 1;  
						Ki = grid[i]

***********************************************************************************************/


///  WAVENOACTUAL: \f$ K = \sqrt{K_x^2 + K_y^2 + K_z^2} \f$
inline DP FFFW_SLAB::Kmagnitude(int lx, int ly, int lz)
{
	if	(global.field.waveno_switch)
		return sqrt(  pow2(Get_kx(lx)*kfactor[1]) + pow2(Get_ky(ly)*kfactor[2]) + pow2(lz*kfactor[3]) );
	
	else 
		return sqrt(  pow2(Get_kx(lx)) + pow2(Get_ky(ly)) + pow2(lz) ); 
}

/// WAVENOACTUAL -- Radius of the smallest sphere that contains the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int FFFW_SLAB::Min_radius_outside() 
{
	if	(global.field.waveno_switch) 
		return (int) ceil(sqrt(  pow2((N[1]/2)*kfactor[1]) + pow2(N[2]/2*kfactor[2]) + pow2(N[3]/2*kfactor[3]) ));
	
	else 
		return (int) ceil(sqrt(pow2(N[1]/2) + pow2(N[2]/2) + pow2(N[3]/2)));
		
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int FFFW_SLAB::Max_radius_inside() 
{
	
	int ans = 1;
	DP Kmag;
	
	if	(global.field.waveno_switch) {
		if (N[2] > 1) {
			Kmag = min( (N[1]/2)*kfactor[1], (N[2]/2)*kfactor[2]);
			Kmag = min(Kmag, (N[3]/2)*kfactor[3]); 
		} 
		
		else if (N[2] == 1)
			Kmag = min((N[1]/2)*kfactor[1], (N[3]/2)*kfactor[3]);
		
		ans = ((int) Kmag);
	}
	
	else  {
		if (N[2] > 1) {
			ans = min(N[1]/2, N[2]/2);
			ans = min(ans, N[3]/2);
		}
		
		else if (N[2] == 1)
			ans = min(N[1]/2, (N[3]/2));
	}
	
	return ans;
	
}


/*! \brief WAVENOACTUAL -- Returns the approximate number of modes in a wavenumber K semishell 
 *							(upper hemisphere) of radius "radius".
 * 
 *  We divide the area in Fourier space by volume of unit lattice \f$ \Pi_i f_i \f$, 
 *	where \f$ f_i \f$ is the factor[i].
 *  In 1D, The area in NOT divided by kfactor.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 * or gridwavenumber depending on the switch
 * WAVENOACTUAL or WAVENOGRID.   Actual wavenumber is preferred.
 *
 * \param  radius
 * \return The number of modes in a hemispheric shell of radius. 
 *			In 3D, it is hemisphere sphere with (kz>=0).
 */
inline DP FFFW_SLAB::Approx_number_modes_in_shell(int radius)
{
	DP ans = 0.0;
	
	if (N[2] > 1)	{
		if	(global.field.waveno_switch)
			ans = (4*M_PI*radius*radius)/(kfactor[1]*kfactor[2]*kfactor[3]);	
	
		else 
			ans = (4*M_PI*radius*radius);
	}
    
	else if (N[2] == 1) 	{
		if	(global.field.waveno_switch)
			ans = (2*M_PI*radius)/(kfactor[1]*kfactor[3]);
		
		else 
			ans = (2*M_PI*radius);
	}
	
	return ans;
}

//*********************************************************************************************



/*! \brief Returns multiplication factor for computing enregy spectrum etc. 
 * 
 *  \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$, but in simulation we
 *  double the energy of most of the modes because complex conjugates modes with 
 * -ky are not stored in the simulation.  
 * The modes on the x-axis are not doubled because their c.c. are already counted.
 * 
 * \param  lx, ly
 * \return Multiplication factor to \f$ |\vect{u}(\vect{k})|^2  \f$ for computing 
 *		energy spectrum etc. factor = 1 implies that the modal energy is already doubled.
 */
 
inline DP FFFW_SLAB::Multiplicity_factor(int lx, int ly, int lz)
{
	DP factor;

	int kx = Get_kx(lx);
	int ky = Get_ky(ly);
	
	if (lz > 0)
		factor = 1.0;
		
	else
		factor = 0.5;
		
	if ( (kx == N[1]/2) || ((ky == N[2]/2)&& (N[2] > 1)) )
		return 2*factor;		// for kx = -N[1]/2 or ky = -N[2]/2 ; 
								// Ignoring corner which would have factor 4
	else
		return factor;
}



/// Modal energy -- \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$
inline DP FFFW_SLAB::Modal_energy(int lx, int ly, int lz, Array<complx,3> A)
{
	return pow2(abs(A(lx,ly,lz)))/2;	
}

/**********************************************************************************************

	Compute Modal helicity for grid index (i1, i2, i3) = K. [Vr x Vi].

***********************************************************************************************/

inline DP FFFW_SLAB::Get_Modal_helicity
(
	int lx, int ly, int lz,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
)
{

	TinyVector<DP,3> Vreal, Vimag, VrcrossVi,  K;
	
	Vreal = real(Ax(lx, ly, lz)), real(Ay(lx, ly, lz)), real(Az(lx, ly, lz));
	Vimag = imag(Ax(lx, ly, lz)), imag(Ay(lx, ly, lz)), imag(Az(lx, ly, lz));
		
	VrcrossVi = cross(Vreal, Vimag);
	Wavenumber(lx, ly, lz, K);
	
	return (dot( K, VrcrossVi));	

}



/**********************************************************************************************=

	Compute Modal Vorticity for grid index (i1, i2, i3) = I K x V(k).

***********************************************************************************************/


inline void FFFW_SLAB::Compute_Modal_vorticity
(
	int lx, int ly, int lz, 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	TinyVector<complx,3> &vorticity
)
{

	TinyVector<DP,3>  K;
	TinyVector<complx,3> Vi;
	
	Vi = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);

	Wavenumber(lx, ly, lz,  K);
	
	vorticity(0) = I *  ( K(1)*Vi(2) -  K(2)*Vi(1));
	vorticity(1) = I *  ( K(2)*Vi(0) -  K(0)*Vi(2));
	vorticity(2) = I *  ( K(0)*Vi(1) -  K(1)*Vi(0));
}




inline void FFFW_SLAB::Compute_Modal_vorticity_y_component
(
    int lx, int ly, int lz, 
    Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
    complx &vort_y
)
{
	TinyVector<DP,3>  K;
	TinyVector<complx,3> Vi;
	
	Vi = Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz);
	
	Wavenumber(lx, ly, lz,  K);
	
	vort_y = I *  ( K(2)*Vi(0) -  K(0)*Vi(2));
}



//*********************************************************************************************	

			
/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_{||} = K_1 \f$			
inline DP FFFW_SLAB::AnisKpll(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_kx(lx) * kfactor[1]); 
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_ky(ly) * kfactor[2]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (lz * kfactor[3]);
	
	else
		return 0;		// for -Wall
		
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_\perp =\sqrt{K_2^2 + K_3^2} \f$			
inline DP FFFW_SLAB::AnisKperp(int lx, int ly, int lz)
{
	if (global.field.anisotropy_dirn == 1)
		return sqrt( pow2(Get_ky(ly) * kfactor[2]) + pow2(lz* kfactor[3]) ); 
	
	else if (global.field.anisotropy_dirn == 2)
		return sqrt( pow2(Get_kx(lx) * kfactor[1]) + pow2(lz* kfactor[3]) );
		
	else if (global.field.anisotropy_dirn == 3)
		return sqrt( pow2(Get_kx(lx) * kfactor[1]) + pow2(Get_ky(ly) * kfactor[2]) );
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 1, \f$ K_{h1} = K_2 \f$										
inline DP FFFW_SLAB::AnisKh1(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_ky(ly) * kfactor[2]);  
	
	else if (global.field.anisotropy_dirn == 2)
		return (lz * kfactor[3]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_kx(lx) * kfactor[1]);
	
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 2, \f$ K_{h2} = K_3 \f$				
inline DP FFFW_SLAB::AnisKh2(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (lz * kfactor[3]);  
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_kx(lx)  * kfactor[1]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_ky(ly) * kfactor[2]); 
	
	else
		return 0;		// for -Wall
}
			
			
/// Cylindrical: Anis_min_Kpll
inline DP FFFW_SLAB::Anis_min_Kpll() 
{ 
	return 0.0;
}
				
/// Cylindrical: Anis_max_Kpll
inline DP FFFW_SLAB::Anis_max_Kpll() 
{ 
	DP maxKpll = 0;
		
    if (global.field.anisotropy_dirn == 1)
        maxKpll = ((N[1]/2) * kfactor[1]); 
    
    else if (global.field.anisotropy_dirn == 2)
        maxKpll = ((N[2]/2) * kfactor[2]); 
        
    else if (global.field.anisotropy_dirn == 3)
        maxKpll = ((N[3]/2) * kfactor[3]);
	
	return maxKpll;
}

/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box
inline int FFFW_SLAB::Anis_max_Krho_radius_inside() 	
{
	DP Kmag = 0.0;
	
    if (global.field.anisotropy_dirn == 1)
        Kmag = min( (N[2]/2)*kfactor[2], (N[3]/2)*kfactor[3] );
    
    else if (global.field.anisotropy_dirn == 2)
        Kmag = min( (N[1]/2)*kfactor[1], (N[3]/2)*kfactor[3] );  
    
    else if (global.field.anisotropy_dirn == 3)
        Kmag = min( (N[1]/2)*kfactor[1], (N[2]/2)*kfactor[2] ); 

    return ((int) Kmag);
}

// Max polar angle
inline DP FFFW_SLAB::Get_max_polar_angle() 
{ 
	
	return M_PI/2;
	
}			

//*********************************************************************************************	


/*! \brief Returns the angle K vector makes with the anisotropic axis 
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(K_{\perp}/K_{||}) \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline DP FFFW_SLAB::AnisKvect_polar_angle(int lx, int ly, int lz)
{
	DP kkpll, kkperp;
	
	kkpll = AnisKpll(lx, ly, lz);
	kkperp = AnisKperp(lx, ly, lz);
	
	return Get_polar_angle(kkperp, kkpll);
}


/*! \brief 3D: Returns the azimutal angle.
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(Ky}/Kx \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline DP FFFW_SLAB::AnisKvect_azimuthal_angle(int lx, int ly, int lz)
{
	
	DP kkh1 = AnisKh1(lx, ly, lz);
	DP kkh2 = AnisKh2(lx, ly, lz);
	
	return Get_azimuthal_angle(kkh1, kkh2);
}			


//********************************** END of four_inline.h *************************************



