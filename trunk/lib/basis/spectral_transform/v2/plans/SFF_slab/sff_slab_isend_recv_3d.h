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


#include "spectral_plan.h"

class SFF_slab_Isend_Recv_3D : public SpectralPlan {
private:
	FFTW_PLAN plan_sintr_x;
	FFTW_PLAN plan_costr_x;
	FFTW_PLAN plan_isintr_x;
	FFTW_PLAN plan_icostr_x;
	FFTW_PLAN plan_ft_r2c_z;
	FFTW_PLAN plan_ift_c2r_z;
	FFTW_PLAN plan_ft_c2c_y;
	FFTW_PLAN plan_ift_c2c_y;

	void Normalize(Array<complx,3> A);

	void SinCostr_x(char sincostr_option, Array<DP,3> Ar);
	void ISinCostr_x(char sincostr_option, Array<DP,3> Ar);
	void FT_r2c_z_and_Isend_yz(Array<DP,3> Ar);
	void Recv_yz_and_FT_c2c_y(Array<complx,3> A);
	void IFT_c2c_y_and_Isend_yz(Array<complx,3> A);
	void Recv_yz_and_IFT_c2r_z(Array<DP,3> Ar);

public:

	SFF_slab_Isend_Recv_3D(int my_id, int numprocs, int N0, int N1, int N2);

	void Forward_transform(string sincostr_option, Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform(string sincostr_option, Array<complx,3> A, Array<DP,3> Ar);

	void Transpose(Array<DP,3> Ar, Array<complx,3> A);
	void Transpose(Array<complx,3> A, Array<DP,3> Ar);
};

