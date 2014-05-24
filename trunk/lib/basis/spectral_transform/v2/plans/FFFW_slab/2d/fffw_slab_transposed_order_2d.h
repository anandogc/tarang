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

class FFFW_slab_transposed_order_2D : public SpectralPlan_Slab_2D {
private:

	FFTW_PLAN plan_ft_r2c_xz;
	FFTW_PLAN plan_ift_c2r_xz;

	DP f(int rx, int rz);
	void Init_array();

	void Normalize(Array<complx,2> A);
public:

	FFFW_slab_transposed_order_2D(int my_id, int numprocs, int num_iter, int Nx, int Nz);

	void Forward_transform(Array<DP,2> Ar, Array<complx,2> A);
	void Inverse_transform(Array<complx,2> A, Array<DP,2> Ar);

	void Transpose(Array<DP,2> Ar, Array<complx,2> A);
	void Transpose(Array<complx,2> A, Array<DP,2> Ar);
};

