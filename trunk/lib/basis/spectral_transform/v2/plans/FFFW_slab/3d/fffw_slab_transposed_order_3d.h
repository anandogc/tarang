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


#include "spectral_plan_slab_3d.h"

class FFFW_slab_transposed_order_3D : public SpectralPlan_Slab_3D {
private:

	FFTW_PLAN plan_ft_r2c_xyz;
	FFTW_PLAN plan_ift_c2r_xyz;

	DP f(int rx, int ry, int rz);
	void Init_array();


	void Normalize(Array<complx,3> A);
public:

	FFFW_slab_transposed_order_3D(int my_id, int numprocs, int num_iter, int Nx, int Ny, int Nz);

	void Forward_transform(Array<DP,3> Ar, Array<complx,3> A);
	void Inverse_transform(Array<complx,3> A, Array<DP,3> Ar);

	void Transpose(Array<DP,3> Ar, Array<complx,3> A);
	void Transpose(Array<complx,3> A, Array<DP,3> Ar);
};

