#include "spectral_transform.h"

void SpectralTransform::Init(string basis, int Nx, int Ny, int Nz, int num_p_rows)
{
	this->basis=basis;
	this->Nx=Nx;
	this->Ny=Ny;
	this->Nz=Nz;
	this->num_p_rows=num_p_rows;
	this->orientation='Z';	//Required by To_slab and To_pencil. By default assume data is along Z

	if (basis[0] == 'F' || basis[0] == 'S') {
		Fx=Nx;
		Ix=Nx;
		Rx=Nx;
	}
	else {
		cerr << "basis[0] can be S or F";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	if (basis[1] == 'F' || basis[1] == 'S') {
		Fy=Ny;
		Iy=Ny;
		Ry=Ny;
	}
	else {
		cerr << "basis[1] can be S or F";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (basis[2] == 'F') {
		Fz=Nz/2+1;
		Iz=Nz/2+1;
		Rz=Nz+2;

		Fz_unit=2;
		Iz_unit=2;
		Rz_unit=2;
	}
	else if (basis[2] == 'S') {
		Fz=Nz/2;
		Iz=Nz/2;
		Rz=Nz;

		Fz_unit=2;
		Iz_unit=2;
		Rz_unit=1;
	}
	else {
		cerr << "basis[] can be S or F";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	set_vars();
	init_transpose();
	init_transform();
}

void SpectralTransform::set_vars() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	num_p_cols=numprocs/num_p_rows;

	MPI_Comm_split(MPI_COMM_WORLD, my_id/num_p_cols, 0, &MPI_COMM_ROW);
	MPI_Comm_split(MPI_COMM_WORLD, my_id%num_p_cols, 0, &MPI_COMM_COL);

	MPI_Comm_rank(MPI_COMM_ROW, &my_col_id);
	MPI_Comm_rank(MPI_COMM_COL, &my_row_id);

	maxfx=Fx;
	maxfy=Fy/num_p_cols;
	maxfz=Fz/num_p_rows;

	maxix=Ix/num_p_cols;
	maxiy=Iy;
	maxiz=Iz/num_p_rows;

	maxrx=Rx/num_p_cols;
	maxry=Ry/num_p_rows;
	maxrz=Rz;


	fx_start=0;
	fy_start=maxfy*my_col_id;
	fz_start=maxfz*my_row_id;

	ix_start=maxix*my_col_id;
	iy_start=0;
	iz_start=maxiz*my_row_id;

	rx_start=maxrx*my_col_id;
	ry_start=maxry*my_row_id;
	rz_start=0;

	//Set Array size
	FA_shape=shape(Fx,Fy/num_p_cols,Fz/num_p_rows);
	IA_shape=shape(Ix/num_p_cols,Iy,Iz/num_p_rows);
	RA_shape=shape(Rx/num_p_cols,Ry/num_p_rows,Rz);

	//FA.resize(FA_shape);
	IA.resize(IA_shape);
	RA.resize(RA_shape);

	fft_plan_Fourier_forward_x=NULL;
	fft_plan_Fourier_inverse_x=NULL;
	fft_plan_Fourier_forward_y=NULL;
	fft_plan_Fourier_inverse_y=NULL;
	fft_plan_Fourier_forward_z=NULL;
	fft_plan_Fourier_inverse_z=NULL;

	fft_plan_Sine_forward_x=NULL;
	fft_plan_Sine_inverse_x=NULL;
	fft_plan_Sine_forward_y=NULL;
	fft_plan_Sine_inverse_y=NULL;
	fft_plan_Sine_forward_z=NULL;
	fft_plan_Sine_inverse_z=NULL;

	fft_plan_Cosine_forward_x=NULL;
	fft_plan_Cosine_inverse_x=NULL;
	fft_plan_Cosine_forward_y=NULL;
	fft_plan_Cosine_inverse_y=NULL;
	fft_plan_Cosine_forward_z=NULL;
	fft_plan_Cosine_inverse_z=NULL;

}


void SpectralTransform::init_transpose()
{
	//To revise
	buffer_Y = new Real[size_t(Nx/num_p_cols)*size_t(Ny/num_p_rows)*size_t(Nz+2)];
	buffer_Z = new Real[size_t(Nx/num_p_cols)*size_t(Ny/num_p_rows)*size_t(Nz+2)];

	plane_size_YX = Iz_unit*(Iz/num_p_rows)*Iy;
	copy_size_YX  = Iz_unit*(Iz/num_p_rows)*(Iy/num_p_cols);
	block_size_YX = Iz_unit*(Iz/num_p_rows)*(Iy/num_p_cols)*(Ix/num_p_cols);

	plane_size_YZ = Iz_unit*(Iz/num_p_rows)*Iy;
	copy_size_YZ  = Iz_unit*(Iz/num_p_rows)*(Iy/num_p_rows);
	block_size_YZ = Iz_unit*(Iz/num_p_rows)*(Iy/num_p_rows)*(Ix/num_p_cols);


	//Configure pencil to slab and vice versa
	int count=maxrx;
	int blocklength=Rz*Ry/num_p_rows;
	int stride=Rz*Ry;
	// int count=1;
	// int blocklength=750;
	// int stride=1;
	MPI_Type_vector(count, blocklength, stride, MPI_Real, &MPI_Vector_ZY_planes);
	MPI_Type_commit(&MPI_Vector_ZY_planes);

	int lower_bound=0;
	int extent=Rz*Ry/num_p_rows*sizeof(Real);
	// int extent=1280;
	MPI_Type_create_resized(MPI_Vector_ZY_planes, lower_bound, extent, &MPI_Vector_resized_ZY_planes);
	MPI_Type_commit(&MPI_Vector_resized_ZY_planes);
/*

	int count=Nx;
	int blocklength=2*maxly;
	int stride=2*Ny;
	MPI_Type_vector(count,blocklength,stride,MPI_Real,&MPI_Vector_y_plane_block);
	MPI_Type_commit(&MPI_Vector_y_plane_block);

	int lower_bound=0;
	int extent=2*maxly*sizeof(Real);
	MPI_Type_create_resized(MPI_Vector_y_plane_block, lower_bound, extent, &MPI_Vector_resized_y_plane_block);
	MPI_Type_commit(&MPI_Vector_resized_y_plane_block);	*/
}

void SpectralTransform::init_transform(){
	//Initialize plans
	int Nx_dims[]={Nx};
	int Ny_dims[]={Ny};
	int Nz_dims[]={Nz};

	fftw_r2r_kind kind[1];

	normalizing_factor=1;

 	//X transforms
	if (basis[0] == 'F') {
		fft_plan_Fourier_forward_x = FFTW_PLAN_MANY_DFT(1, Nx_dims, maxfy*maxfz,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxfy*maxfz, 1,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxfy*maxfz, 1,
			FFTW_FORWARD, FFTW_PLAN_FLAG);

		fft_plan_Fourier_inverse_x = FFTW_PLAN_MANY_DFT(1, Nx_dims, maxfy*maxfz,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxfy*maxfz, 1,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxfy*maxfz, 1,
			FFTW_BACKWARD, FFTW_PLAN_FLAG);

		normalizing_factor/=Nx;

	}
	else if (basis[0] == 'S') {
		kind[0]=FFTW_RODFT10;
		fft_plan_Sine_forward_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Fz_unit*maxfz*maxfy,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_RODFT01;
		fft_plan_Sine_inverse_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Fz_unit*maxfz*maxfy,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT10;
		fft_plan_Cosine_forward_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Fz_unit*maxfz*maxfy,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT01;
		fft_plan_Cosine_inverse_x = FFTW_PLAN_MANY_R2R(1, Nx_dims, Fz_unit*maxfz*maxfy,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Fz_unit*maxfz*maxfy, 1,
		    kind, FFTW_PLAN_FLAG);

		normalizing_factor/=2*Nx;
	}
	else {
		cerr << "SpectralTransform::init_transform() - invalid basis" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}



	//Y transforms
	if (basis[1] == 'F') {
		fft_plan_Fourier_forward_y = FFTW_PLAN_MANY_DFT(1, Ny_dims, maxiz,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxiz, 1,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxiz, 1,
			FFTW_FORWARD, FFTW_PLAN_FLAG);

		fft_plan_Fourier_inverse_y = FFTW_PLAN_MANY_DFT(1, Ny_dims, maxiz,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxiz, 1,
			reinterpret_cast<FFTW_Complex*>(IA.data()), NULL,
			maxiz, 1,
			FFTW_BACKWARD, FFTW_PLAN_FLAG);

		normalizing_factor/=Ny;
	}
	else if (basis[1] == 'S') {
		kind[0]=FFTW_RODFT10;
		fft_plan_Sine_forward_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Iz_unit*maxiz,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_RODFT01;
		fft_plan_Sine_inverse_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Iz_unit*maxiz,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT10;
		fft_plan_Cosine_forward_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Iz_unit*maxiz,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT01;
		fft_plan_Cosine_inverse_y = FFTW_PLAN_MANY_R2R(1, Ny_dims, Iz_unit*maxiz,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    reinterpret_cast<Real*>(IA.data()), NULL,
		    Iz_unit*maxiz, 1,
		    kind, FFTW_PLAN_FLAG);

		normalizing_factor/=2*Ny;
	}
	else {
		cerr << "SpectralTransform::init_transform() - invalid basis" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//Z transforms
	if (basis[2] == 'F') {
		fft_plan_Fourier_forward_z = FFTW_PLAN_MANY_DFT_R2C(1, Nz_dims, maxrx*maxry,
			reinterpret_cast<Real*>(RA.data()), NULL,
			1, Nz+2,
			reinterpret_cast<FFTW_Complex*>(RA.data()), NULL,
			1, Nz/2+1,
			FFTW_PLAN_FLAG);

		fft_plan_Fourier_inverse_z = FFTW_PLAN_MANY_DFT_C2R(1, Nz_dims, maxrx*maxry,
			reinterpret_cast<FFTW_Complex*>(RA.data()), NULL,
			1, Nz/2+1,
			reinterpret_cast<Real*>(RA.data()), NULL,
			1, Nz+2,
			FFTW_PLAN_FLAG);

		normalizing_factor/=Nz;
	}
	else if (basis[2] == 'S') {

		kind[0]=FFTW_RODFT10;
		fft_plan_Sine_forward_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, maxrx*maxry,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_RODFT01;
		fft_plan_Sine_inverse_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, maxrx*maxry,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT10;
		fft_plan_Cosine_forward_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, maxrx*maxry,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    kind, FFTW_PLAN_FLAG);

		kind[0]=FFTW_REDFT01;
		fft_plan_Cosine_inverse_z = FFTW_PLAN_MANY_R2R(1, Nz_dims, maxrx*maxry,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    reinterpret_cast<Real*>(RA.data()), NULL,
		    1, Rz,
		    kind, FFTW_PLAN_FLAG);

		normalizing_factor/=2*Nz;
	}
	else {
		cerr << "SpectralTransform::init_transform() - invalid basis" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}


/*void SpectralTransform::init_timers(){
	timers = new double[22];
	for (int i=0; i<22; i++)
		timers[i]=0;
	num_timers=22;
}*/



void SpectralTransform::Forward_transform(string basis_option, Array<Real,3> Ar, Array<Complex,3> A)
{
	if ( (basis_option[2] == 'F') && (fft_plan_Fourier_forward_z != NULL) )
		FFTW_EXECUTE_DFT_R2C(fft_plan_Fourier_forward_z, reinterpret_cast<Real*>(Ar.data()), reinterpret_cast<FFTW_Complex*>(Ar.data()));
	else if ( (basis_option[2] == 'S') && (fft_plan_Sine_forward_z != NULL) ) {
		FFTW_EXECUTE_R2R(fft_plan_Sine_forward_z, reinterpret_cast<Real*>(Ar.data()), reinterpret_cast<Real*>(Ar.data()));
		ArrayShiftRight(Ar, 'Z');
	}
	else if ( (basis_option[2] == 'C') && (fft_plan_Cosine_forward_z != NULL) )
		FFTW_EXECUTE_R2R(fft_plan_Cosine_forward_z, reinterpret_cast<Real*>(Ar.data()), reinterpret_cast<Real*>(Ar.data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::forward_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}


	Transpose_ZY(Ar, IA);


	if ( (basis_option[1] == 'F') && (fft_plan_Fourier_forward_y != NULL) )
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_DFT(fft_plan_Fourier_forward_y, reinterpret_cast<FFTW_Complex*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_Complex*>(IA(ix,Range::all(),Range::all()).data()));
	else if ( (basis_option[1] == 'S') && (fft_plan_Sine_forward_y != NULL) ) {
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_R2R(fft_plan_Sine_forward_y, reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()));
		ArrayShiftRight(IA, 'Y');
	}
	else if ( (basis_option[1] == 'C') && (fft_plan_Cosine_forward_y != NULL) )
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_R2R(fft_plan_Cosine_forward_y, reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::forward_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}


	Transpose_YX(IA, A);


	if ( (basis_option[0] == 'F') && (fft_plan_Fourier_forward_x != NULL) )
		FFTW_EXECUTE_DFT(fft_plan_Fourier_forward_x, reinterpret_cast<FFTW_Complex*>(A.data()), reinterpret_cast<FFTW_Complex*>(A.data()));
	else if ( (basis_option[0] == 'S') && (fft_plan_Sine_forward_x != NULL) ) {
		FFTW_EXECUTE_R2R(fft_plan_Sine_forward_x, reinterpret_cast<Real*>(A.data()), reinterpret_cast<Real*>(A.data()));
		ArrayShiftRight(A, 'X');
	}
	else if ( (basis_option[0] == 'C') && (fft_plan_Cosine_forward_x != NULL) )
		FFTW_EXECUTE_R2R(fft_plan_Cosine_forward_x, reinterpret_cast<Real*>(A.data()), reinterpret_cast<Real*>(A.data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::forward_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}

	Normalize(A);
}


void SpectralTransform::Inverse_transform(string basis_option, Array<Complex,3> A, Array<Real,3> Ar) 
{
	if ( (basis_option[0] == 'F') && (fft_plan_Fourier_inverse_x != NULL) )
		FFTW_EXECUTE_DFT(fft_plan_Fourier_inverse_x, reinterpret_cast<FFTW_Complex*>(A.data()), reinterpret_cast<FFTW_Complex*>(A.data()));
	else if ( (basis_option[0] == 'S') && (fft_plan_Sine_inverse_x != NULL) ) {
		ArrayShiftLeft(A, 'X');
		FFTW_EXECUTE_R2R(fft_plan_Sine_inverse_x, reinterpret_cast<Real*>(A.data()), reinterpret_cast<Real*>(A.data()));
	}
	else if ( (basis_option[0] == 'C') && (fft_plan_Cosine_inverse_x != NULL) )
		FFTW_EXECUTE_R2R(fft_plan_Cosine_inverse_x, reinterpret_cast<Real*>(A.data()), reinterpret_cast<Real*>(A.data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::inverse_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}


	Transpose_XY(A, IA);


	if ( (basis_option[1] == 'F') && (fft_plan_Fourier_inverse_y != NULL) )
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_DFT(fft_plan_Fourier_inverse_y, reinterpret_cast<FFTW_Complex*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<FFTW_Complex*>(IA(ix,Range::all(),Range::all()).data()));
	else if ( (basis_option[1] == 'S') && (fft_plan_Sine_inverse_y != NULL) ){
		ArrayShiftLeft(IA, 'Y');
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_R2R(fft_plan_Sine_inverse_y, reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()));
	}
	else if ( (basis_option[1] == 'C') && (fft_plan_Cosine_inverse_y != NULL) )
		for (int ix = 0; ix < maxix; ix++)
			FFTW_EXECUTE_R2R(fft_plan_Cosine_inverse_y, reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()), reinterpret_cast<Real*>(IA(ix,Range::all(),Range::all()).data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::inverse_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}


	Transpose_YZ(IA, Ar);


	if ( (basis_option[2] == 'F') && (fft_plan_Fourier_inverse_z != NULL) )
		FFTW_EXECUTE_DFT_C2R(fft_plan_Fourier_inverse_z, reinterpret_cast<FFTW_Complex*>(Ar.data()), reinterpret_cast<Real*>(Ar.data()));
	else if ( (basis_option[2] == 'S') && (fft_plan_Sine_inverse_z != NULL) ){
		ArrayShiftLeft(Ar, 'Z');
		FFTW_EXECUTE_R2R(fft_plan_Sine_inverse_z, reinterpret_cast<Real*>(Ar.data()), reinterpret_cast<Real*>(Ar.data()));
	}
	else if ( (basis_option[2] == 'C') && (fft_plan_Cosine_inverse_z != NULL) )
		FFTW_EXECUTE_R2R(fft_plan_Cosine_inverse_z, reinterpret_cast<Real*>(Ar.data()), reinterpret_cast<Real*>(Ar.data()));
	else {
		if (my_id==0) cerr << "SpectralTransform::forward_transform - Invalid basis_option '" << basis_option << "'" << endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
}


void SpectralTransform::Transpose(Array<Real,3> Ar, Array<Complex,3> A)
{
	Transpose_ZY(Ar,IA);
	Transpose_YX(IA,A);
}

void SpectralTransform::Transpose(Array<Complex,3> A, Array<Real,3> Ar)
{
	Transpose_XY(A,IA);
	Transpose_YZ(IA,Ar);
}


void SpectralTransform::Transpose_XY(Array<Complex,3> FA, Array<Complex,3> IA)
{

	if (num_p_cols>1) {
		Real *__restrict__ FA_data=reinterpret_cast<Real*>(FA.data());
		Real *__restrict__ IA_data=reinterpret_cast<Real*>(IA.data());

		// timers[12] -= MPI_Wtime();
		MPI_Alltoall(FA_data, block_size_YX, MPI_Real, buffer_Y, block_size_YX, MPI_Real, MPI_COMM_ROW);
		// timers[12] += MPI_Wtime();

		// timers[13] -= MPI_Wtime();
		for (int c=0; c<num_p_cols; c++)
			for (int ix=0; ix<maxix; ix++)
				std::copy(buffer_Y + c*block_size_YX + ix*copy_size_YX, buffer_Y + c*block_size_YX + (ix+1)*copy_size_YX, IA_data + ix*plane_size_YX + c*copy_size_YX);
		// timers[13] += MPI_Wtime();
	}
	else
		IA=Array<Complex,3>(reinterpret_cast<Complex*>(FA.data()), IA_shape, neverDeleteData);

	orientation='Y';
}

void SpectralTransform::Transpose_YX(Array<Complex,3> IA, Array<Complex,3> FA)
{
	if (num_p_cols>1) {
		Real *__restrict__ IA_data=reinterpret_cast<Real*>(IA.data());
		Real *__restrict__ FA_data=reinterpret_cast<Real*>(FA.data());


		// timers[14] -= MPI_Wtime();
		for (size_t c=0; c<num_p_cols; c++)
			for (size_t ix=0; ix<maxix; ix++)
				std::copy(IA_data + ix*plane_size_YX + c*copy_size_YX, IA_data + ix*plane_size_YX + (c+1)*copy_size_YX, buffer_Y + c*block_size_YX + ix*copy_size_YX);
		// timers[14] += MPI_Wtime();

		// timers[15] -= MPI_Wtime();
		MPI_Alltoall(buffer_Y, block_size_YX, MPI_Real, FA_data, block_size_YX, MPI_Real, MPI_COMM_ROW);
		// timers[15] += MPI_Wtime();
	}
	else
		FA=Array<Complex,3>(reinterpret_cast<Complex*>(IA.data()), FA_shape, neverDeleteData);

	orientation='X';
}

void SpectralTransform::Transpose_YZ(Array<Complex,3> IA, Array<Real,3> RA)
{
	if (num_p_rows>1) {
		Real *__restrict__ IA_data=reinterpret_cast<Real*>(IA.data());
		Real *__restrict__ Ar_data=reinterpret_cast<Real*>(RA.data());

		// timers[16] -= MPI_Wtime();
		for (size_t r=0; r<num_p_rows; r++)
			for (size_t ix=0; ix<maxix; ix++)
				std::copy(IA_data + ix*plane_size_YZ + r*copy_size_YZ, IA_data + ix*plane_size_YZ + (r+1)*copy_size_YZ, buffer_Y + r*block_size_YZ + ix*copy_size_YZ);
		// timers[16] += MPI_Wtime();


		// timers[17] -= MPI_Wtime();
		MPI_Alltoall(buffer_Y, block_size_YZ, MPI_Real, buffer_Z, block_size_YZ, MPI_Real, MPI_COMM_COL);
		// timers[17] += MPI_Wtime();


		// timers[18] -= MPI_Wtime();
		for (size_t r=0; r<num_p_rows; r++)
			for (size_t j=0; j<maxry*maxrx; j++)
				std::copy(buffer_Z + r*block_size_YZ + j*Iz_unit*maxiz,buffer_Z + r*block_size_YZ + (j+1)*Iz_unit*maxiz, Ar_data + j*Rz + r*Iz_unit*maxiz);
		// timers[18] += MPI_Wtime();
	}
	else 
		RA=Array<Real,3>(reinterpret_cast<Real*>(IA.data()), RA_shape, neverDeleteData);


 	if (Iz_unit*maxiz*num_p_rows < Rz) {
		// timers[11] -= MPI_Wtime();
		Zero_pad_last_plane(RA);
		// timers[11] += MPI_Wtime();
	}

	orientation='Z';
}

void SpectralTransform::Transpose_ZY(Array<Real,3> RA, Array<Complex,3> IA)
{
	if (num_p_rows>1) {
		Real *__restrict__ RA_data=reinterpret_cast<Real*>(RA.data());
		Real *__restrict__ IA_data=reinterpret_cast<Real*>(IA.data());

		// timers[19] -= MPI_Wtime();
		for (size_t r=0; r<num_p_rows; r++)
			for (size_t j=0; j<maxry*maxrx; j++)
				std::copy(RA_data + j*Rz + r*Iz_unit*maxiz, RA_data + j*Rz + (r+1)*Iz_unit*maxiz, buffer_Z + r*block_size_YZ + j*Iz_unit*maxiz);
		// timers[19] += MPI_Wtime();

		// timers[20] -= MPI_Wtime();
		MPI_Alltoall(buffer_Z, block_size_YZ, MPI_Real, buffer_Y, block_size_YZ, MPI_Real, MPI_COMM_COL);
		// timers[20] += MPI_Wtime();

		// timers[21] -= MPI_Wtime();
		for (size_t r=0; r<num_p_rows; r++)
			for (size_t j=0; j<maxix; j++)
				std::copy(buffer_Y + r*block_size_YZ + j*copy_size_YZ, buffer_Y + r*block_size_YZ + (j+1)*copy_size_YZ, IA_data + j*plane_size_YZ + r*copy_size_YZ);
		// timers[21] += MPI_Wtime();
	}
	else
		IA=Array<Complex,3>(reinterpret_cast<Complex*>(RA.data()), FA_shape, neverDeleteData);

	orientation='Y';
}

void SpectralTransform::Zero_pad_last_plane(Array<Real,3> Ar) {
	Ar(Range::all(),Range::all(),Range(Nz,Nz+1))=0;
}

void  SpectralTransform::Normalize(Array<Complex,3> A) {
	A*=normalizing_factor;
}


//After SinCos transform, shift right
void SpectralTransform::ArrayShiftRight_basic(void *data, TinyVector<int,3> shape, int axis)
{	
	Array<Real,3> Ar((Real*)(data), shape, neverDeleteData);

	if (axis == 'X') {
		Ar(Range(Ar.extent(0)-1,1,-1),Range::all(),Range::all()) = Ar(Range(Ar.extent(0)-2,0,-1),Range::all(),Range::all());
		Ar(0,Range::all(),Range::all()) = 0.0;
	}
	else if (axis == 'Y') {
		Ar(Range::all(),Range(Ar.extent(1)-1,1,-1),Range::all()) = Ar(Range::all(),Range(Ar.extent(1)-2,0,-1),Range::all());
		Ar(Range::all(),0,Range::all()) = 0.0;	
	}
	else if (axis == 'Z') {
		Ar(Range::all(),Range::all(),Range(Ar.extent(2)-1,1,-1)) = Ar(Range::all(),Range::all(),Range(Ar.extent(2)-2,0,-1));
		Ar(Range::all(),Range::all(),0) = 0.0;	
	}
}

void SpectralTransform::ArrayShiftRight(Array<Real,3> Ar, int axis)
{	
	ArrayShiftRight_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralTransform::ArrayShiftRight(Array<Complex,3> A, int axis)
{	
	ArrayShiftRight_basic(A.data(), A.shape()*shape(1,1,2), axis);
}

//Before inverse SinCos transform, shift left
void SpectralTransform::ArrayShiftLeft_basic(void *data, TinyVector<int,3> shape, int axis)
{
	Array<Real,3> Ar((Real*)(data), shape, neverDeleteData);

	if (axis == 'X') {
		Ar(Range(0,Ar.extent(0)-2),Range::all(),Range::all()) = Ar(Range(1,Ar.extent(0)-1),Range::all(),Range::all());
		Ar(Ar.extent(0)-1,Range::all(),Range::all()) = 0.0;	
	}
	else if (axis == 'Y') {
		Ar(Range::all(),Range(0,Ar.extent(1)-2),Range::all()) = Ar(Range::all(),Range(1,Ar.extent(1)-1),Range::all());
		Ar(Range::all(),Ar.extent(1)-1,Range::all()) = 0.0;

	}
	else if (axis == 'Z') {
		Ar(Range::all(),Range::all(),Range(0,Ar.extent(2)-2)) = Ar(Range::all(),Range::all(),Range(1,Ar.extent(2)-1));
		Ar(Range::all(),Range::all(),Ar.extent(2)-1) = 0.0;	
	}
}

void SpectralTransform::ArrayShiftLeft(Array<Real,3> Ar, int axis)
{	
	ArrayShiftLeft_basic(Ar.data(), Ar.shape(), axis);
}

void SpectralTransform::ArrayShiftLeft(Array<Complex,3> A, int axis)
{	
	ArrayShiftLeft_basic(A.data(), A.shape()*shape(1,1,2), axis);
}


void SpectralTransform::To_slab(Array<Real,3> ArPencil, Array<Real,3> ArSlab)
{
	/*MPI_Gather(ArPencil.data(), ArPencil.size(), MPI_Real,
           ArSlab.data(), 1, MPI_Vector_resized_ZY_planes, 0,
           MPI_COMM_COL);
	*/

	Array<Real,3> ArSlab_temp(num_p_rows*maxrx,Ry/num_p_rows,Rz);

	MPI_Gather(ArPencil.data(), ArPencil.size(), MPI_Real,
           ArSlab_temp.data(), ArPencil.size(), MPI_Real, 0,
           MPI_COMM_COL);


	if (my_row_id==0) {
		for (int p=0; p<num_p_rows; p++) {
			ArSlab(Range::all(), Range(p*Ry/num_p_rows,(p+1)*Ry/num_p_rows-1), Range::all())
				= ArSlab_temp(Range(p*maxrx,(p+1)*maxrx-1),Range::all(),Range::all());
		}
	}
}

void SpectralTransform::To_pencil(Array<Real,3> ArSlab, Array<Real,3> ArPencil)
{
	Array<Real,3> ArSlab_temp(num_p_rows*maxrx,Ry/num_p_rows,Rz);

	if (my_row_id==0) {
		for (int p=0; p<num_p_rows; p++) {
			ArSlab_temp(Range(p*maxrx,(p+1)*maxrx-1),Range::all(),Range::all())
				= ArSlab(Range::all(), Range(p*Ry/num_p_rows,(p+1)*Ry/num_p_rows-1), Range::all());
		}
	}

	MPI_Scatter(ArSlab_temp.data(), ArPencil.size(), MPI_Real,
		ArPencil.data(), ArPencil.size(), MPI_Real, 0,
        MPI_COMM_COL);
}
MPI_Comm SpectralTransform::Get_communicator(string which) {
	if (which=="ROW")
		return MPI_COMM_ROW;
	else if (which=="COL")
		return MPI_COMM_COL;
	else {
		if (my_id==0)
			cerr << "Invalid communicator: " << which << endl;
		MPI_Finalize();
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
}


TinyVector<int,3> SpectralTransform::Get_FA_shape()
{
	return FA_shape;
}
TinyVector<int,3> SpectralTransform::Get_IA_shape()
{
	return IA_shape;
}
TinyVector<int,3> SpectralTransform::Get_RA_shape()
{
	return RA_shape;
}

TinyVector<int,3> SpectralTransform::Get_FA_start()
{
	return TinyVector<int, 3>(fx_start,fy_start,fz_start);
}
TinyVector<int,3> SpectralTransform::Get_IA_start()
{
	return TinyVector<int, 3>(ix_start,iy_start,iz_start);
}
TinyVector<int,3> SpectralTransform::Get_RA_start()
{
	return TinyVector<int, 3>(rx_start,ry_start,rz_start);
}

int SpectralTransform::Get_row_id()
{
	return my_row_id;
}
int SpectralTransform::Get_col_id()
{
	return my_col_id;
}

int SpectralTransform::Get_num_p_rows()
{
	return num_p_rows;
}
int SpectralTransform::Get_num_p_cols()
{
	return num_p_cols;
}

void SpectralTransform::Finalize() {
	delete[] buffer_Y;
	delete[] buffer_Z;

	MPI_Type_free(&MPI_Vector_ZY_planes);
	MPI_Type_free(&MPI_Vector_resized_ZY_planes);

	FFTW_DESTROY_PLAN(fft_plan_Fourier_forward_x);
	FFTW_DESTROY_PLAN(fft_plan_Fourier_inverse_x);
	FFTW_DESTROY_PLAN(fft_plan_Sine_forward_x);
	FFTW_DESTROY_PLAN(fft_plan_Sine_inverse_x);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_forward_x);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_inverse_x);
	FFTW_DESTROY_PLAN(fft_plan_Fourier_forward_y);
	FFTW_DESTROY_PLAN(fft_plan_Fourier_inverse_y);
	FFTW_DESTROY_PLAN(fft_plan_Sine_forward_y);
	FFTW_DESTROY_PLAN(fft_plan_Sine_inverse_y);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_forward_y);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_inverse_y);
	FFTW_DESTROY_PLAN(fft_plan_Fourier_forward_z);
	FFTW_DESTROY_PLAN(fft_plan_Fourier_inverse_z);
	FFTW_DESTROY_PLAN(fft_plan_Sine_forward_z);
	FFTW_DESTROY_PLAN(fft_plan_Sine_inverse_z);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_forward_z);
	FFTW_DESTROY_PLAN(fft_plan_Cosine_inverse_z);
}
