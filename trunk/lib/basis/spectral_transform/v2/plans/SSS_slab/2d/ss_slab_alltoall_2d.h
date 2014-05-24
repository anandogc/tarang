/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "spectral_plan_slab_2d.h"

class SS_slab_Alltoall_2D : public SpectralPlan_Slab_2D {
private:

	FFTW_PLAN plan_sintr_x, plan_isintr_x;
	FFTW_PLAN plan_costr_x, plan_icostr_x;

	FFTW_PLAN plan_sintr_z, plan_isintr_z;
	FFTW_PLAN plan_costr_z, plan_icostr_z;

	DP f(string sincostr_option, int rx, int rz);
	void Init_array();
	
	void SinCostr_x(char sincostr_option, Array<DP,2> Ar);
	void ISinCostr_x(char sincostr_option, Array<DP,2> Ar);
	void SinCostr_z(char sincostr_option, Array<complx,2> Ar);
	void ISinCostr_z(char sincostr_option, Array<complx,2> Ar);


	void Normalize(Array<complx,2> A);
public:

	SS_slab_Alltoall_2D(int my_id, int numprocs, int num_iter, int Nx, int Nz);

	void Forward_transform(string sincostr_switch, Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform(string sincostr_switch, Array<complx,2> A, Array<DP,2> Ar);

	void Transpose(Array<DP,2> Ar, Array<complx,2> A);
	void Transpose(Array<complx,2> A, Array<DP,2> Ar);
};

