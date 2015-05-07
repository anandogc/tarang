#ifndef _SPECTRAL_TRANSFORM_H
#define _SPECTRAL_TRANSFORM_H

#include "spectral_def_vars.h"

class SpectralTransform {
private:
	string basis;

	int my_id;
	int numprocs;

	MPI_Comm MPI_COMM_ROW;
	MPI_Comm MPI_COMM_COL;

	int my_row_id;
	int my_col_id;

	int Nx,Ny,Nz;

	int Fx,Fy,Fz;
	int Ix,Iy,Iz;
	int Rx,Ry,Rz;

	int Fz_unit, Iz_unit, Rz_unit;

	int num_p_rows;
	int num_p_cols;

	//Set Array size
	TinyVector<int,3> FA_shape, IA_shape, RA_shape;

	int maxfx, maxfy, maxfz;
	int maxix, maxiy, maxiz;
	int maxrx, maxry, maxrz;

	int fx_start, fy_start, fz_start;
	int ix_start, iy_start, iz_start;
	int rx_start, ry_start, rz_start;

	FFTW_PLAN fft_plan_Fourier_forward_x;
	FFTW_PLAN fft_plan_Fourier_inverse_x;
	FFTW_PLAN fft_plan_Fourier_forward_y;
	FFTW_PLAN fft_plan_Fourier_inverse_y;
	FFTW_PLAN fft_plan_Fourier_forward_z;
	FFTW_PLAN fft_plan_Fourier_inverse_z;

	FFTW_PLAN fft_plan_Sine_forward_x;
	FFTW_PLAN fft_plan_Sine_inverse_x;
	FFTW_PLAN fft_plan_Sine_forward_y;
	FFTW_PLAN fft_plan_Sine_inverse_y;
	FFTW_PLAN fft_plan_Sine_forward_z;
	FFTW_PLAN fft_plan_Sine_inverse_z;

	FFTW_PLAN fft_plan_Cosine_forward_x;
	FFTW_PLAN fft_plan_Cosine_inverse_x;
	FFTW_PLAN fft_plan_Cosine_forward_y;
	FFTW_PLAN fft_plan_Cosine_inverse_y;
	FFTW_PLAN fft_plan_Cosine_forward_z;
	FFTW_PLAN fft_plan_Cosine_inverse_z;

	//MPI
	MPI_Datatype MPI_Vector_ZY_planes;
	MPI_Datatype MPI_Vector_resized_ZY_planes;

	Real *buffer_Y;
	Real *buffer_Z;

	size_t plane_size_YX;
	size_t copy_size_YX;
	size_t block_size_YX;

	size_t plane_size_YZ;
	size_t block_size_YZ;
	size_t copy_size_YZ;
	//


	Array<Complex,3> FA;
	Array<Complex,3> IA;
	Array<Real,3> RA;

	Real normalizing_factor;

	char orientation;

	void set_vars();
	void init_transpose();
	void init_transform();
	void init_timers();
	void Transpose_XY(Array<Complex,3> A1, Array<Complex,3> A2);
	void Transpose_YX(Array<Complex,3> A1, Array<Complex,3> A2);
	void Transpose_YZ(Array<Complex,3> A, Array<Real,3> Ar);
	void Transpose_ZY(Array<Real,3> Ar, Array<Complex,3> A);

	void ArrayShiftRight_basic(void *data, TinyVector<int,3> shape, int axis);
	void ArrayShiftRight(Array<Real,3> Ar, int axis);
	void ArrayShiftRight(Array<Complex,3> A, int axis);
	void ArrayShiftLeft_basic(void *data, TinyVector<int,3> shape, int axis);
	void ArrayShiftLeft(Array<Real,3> Ar, int axis);
	void ArrayShiftLeft(Array<Complex,3> A, int axis);

public:
	void Init(string basis, int Nx, int Ny, int Nz, int num_p_rows);
	void Finalize();
	
	void Forward_transform(string basis_option, Array<Real,3> RA, Array<Complex,3> FA);
	void Inverse_transform(string basis_option, Array<Complex,3> FA, Array<Real,3> RA);
	
	void Transpose(Array<Real,3> Ar, Array<Complex,3> A);
	void Transpose(Array<Complex,3> A, Array<Real,3> Ar);
	
	void Zero_pad_last_plane(Array<Real,3> Ar);
	void Normalize(Array<Complex,3> A);

	void To_slab(Array<Real,3> ArPencil, Array<Real,3> ArSlab);
	void To_pencil(Array<Real,3> ArSlab, Array<Real,3> ArPencil);

	MPI_Comm Get_communicator(string which);

	TinyVector<int,3> Get_FA_shape();
	TinyVector<int,3> Get_IA_shape();
	TinyVector<int,3> Get_RA_shape();

	TinyVector<int,3> Get_FA_start();
	TinyVector<int,3> Get_IA_start();
	TinyVector<int,3> Get_RA_start();

	int Get_row_id();
	int Get_col_id();

	int Get_num_p_rows();
	int Get_num_p_cols();

};
#endif
