#include <mpi.h>
#include <fftw3.h>

#include "spectral_transform.h"


void test(int Nx, int Ny, int Nz, int iter);
Real f(string basis, string basis_option, int rx, int ry, int rz);
void Print_large_Fourier_elements(Array<Complex,3> A, string array_name);

int my_id,numprocs, Nx,Ny,Nz;

int main(int argc, char** argv)
{

  	MPI_Init(&argc, &argv);

  	int Nx=atoi(argv[1]);
  	int Ny=atoi(argv[2]);
  	int Nz=atoi(argv[3]);

  	int iter=atoi(argv[4]);

  	test(Nx, Ny, Nz, iter);

	MPI_Finalize();
	return 0;
}


void test(int Nx, int Ny, int Nz, int iter)
{
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	string basis="FFF";
	string basis_option="FFF";

	//Nx=Ny=N;
	//Nz=N;

	int rows=1;

  	if (Nx<numprocs)		//For pencil
  		rows=numprocs/Nx;

	SpectralTransform spectralTransform;
	spectralTransform.Init(basis,Nx,Ny,Nz,rows);

	Array<Complex,3> A(spectralTransform.Get_FA_shape());
	Array<Real,3> Ar(spectralTransform.Get_RA_shape());
	Array<Real, 3> Ar_init(spectralTransform.Get_RA_shape());

	Ar_init=0;

	int maxrx=Ar.extent(0);
	int maxry=Ar.extent(1);
	int maxrz=Ar.extent(2);


	TinyVector<int,3> RA_start = spectralTransform.Get_RA_start();
	int rx_start=RA_start(0);
	int ry_start=RA_start(1);
	int rz_start=RA_start(2);

	for (int rx=0; rx<maxrx; rx++)
		for (int ry=0; ry<maxry; ry++)
			for (int rz=0; rz<maxrz; rz++){
				Ar_init(rx, ry, rz) = f(basis, basis_option, rx_start + rx, ry_start + ry, rz);
				// Ar_init(rx, ry, rz) =  my_id*Ar_init.size() + (rx*Ar_init.extent(1)*Ar_init.extent(2)) + ry*Ar_init.extent(2) + rz;
			}

	//Zero pad
	if (basis[2] != 'S')
		spectralTransform.Zero_pad_last_plane(Ar_init);

/*	if (my_id==1)
		sleep(1);
	cout << Ar_init << endl;	*/	


	// double local_energy = sum(sqr(Ar_init));
	// double total_energy;
	// MPI_Reduce(&local_energy, &total_energy, 1, MPI_Real, MPI_SUM, 0, MPI_COMM_WORLD);

/*	if (my_id == 0)
		cout << "Before transform energy = " << total_energy << endl;*/

    Real start_time;
    Real end_time;

	Ar = Ar_init;

    Real local_time=0, local_time_max;
    MPI_Barrier(MPI_COMM_WORLD);


	time_t start;
	time_t end;

	time (&start);

	for (int i=0; i<iter; i++) {
		spectralTransform.Forward_transform(basis_option, Ar, A);
		spectralTransform.Inverse_transform(basis_option, A, Ar);
	}
	time(&end);

	Real dif = difftime(end, start);
	if (my_id == 0) {
		cout << "Time taken = " << dif << endl;
		cout << "sum(abs(Ar-Ar_init)) = " << sum(abs(Ar-Ar_init)) << endl;

	}
}

void Print_large_Fourier_elements(Array<Complex,3> A, string array_name)
{

	for (int lx=0; lx<A.extent(0); lx++)
		for (int ly=0; ly<A.extent(1); ly++)
			for (int lz=0; lz<A.extent(2); lz++)
				if (abs(A(lx, ly, lz)) > 1E-6) {
					cout << "my_id = " << my_id <<  " vect(k) = (" << lx << "," << ly << "," << lz <<");  " << array_name << "(k) = " << A(lx, ly, lz) << '\n';
				}
	
	 if (my_id==0)
		 cout << endl;
}

Real f(string basis, string basis_option, int rx, int ry, int rz)
{
	Real k0 = 1;
	Real x,y,z;
	
	Real Lx=2*M_PI;
	Real Ly=2*M_PI;
	Real Lz=2*M_PI;
	if (basis == "FFF") {
		x = rx*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
		return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
	}

	else if (basis == "FFFW") {
		x = rx*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
		return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis == "SFF") {
		Lx=M_PI;
		x = (rx+0.5)*Lx/Nx;
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;

		if (basis_option=="SFF")
			return 8*sin(k0*x)*cos(k0*y)*cos(k0*z);
		else if (basis_option=="CFF")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis == "SSF") {
		x = (rx+0.5)*Lx/Nx;
		y = (ry+0.5)*Ly/Ny;
		z = rz*Lz/Nz;

		if (basis_option=="SSF")
			return 8*sin(k0*x)*sin(k0*y)*cos(k0*z);
		else if (basis_option=="CCF")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis == "SSS") {
		x = (rx+0.5)*Lx/Nx;
		y = (ry+0.5)*Ly/Ny;
		z = (rz+0.5)*Lz/Nz;

		if (basis_option=="SSS")
			return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
		else if (basis_option=="CCC")
			return 8*cos(k0*x)*cos(k0*y)*cos(k0*z);
	}
	
	else if (basis == "ChFF") {
		x = rx*M_PI/(Nx-1);
		y = ry*Ly/Ny;
		z = rz*Lz/Nz;
	}

	
	
//	return 4*sin(k0*x)*cos(k0*z);
//	return 4*(2*cos(k0*x)*cos(k0*x)-1)*cos(k0*y)*cos(k0*z); // Cheb
//	return 4*(cos(k0*x))*cos(k0*y)*cos(k0*z);  // Cheb
//	return 2*(cos(k0*x))*cos(k0*z);  // Cheb (101)
	// return 8*sin(k0*x)*sin(k0*y)*sin(k0*z);
	// return 4*sin(k0*x)*sin(k0*z);
}
