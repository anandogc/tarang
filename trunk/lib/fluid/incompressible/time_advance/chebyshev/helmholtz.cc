//#include "basis_basicfn_inline.h"
//#include <gsl/gsl_inline.h>
// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_errno.h>

// Tridogonal solved using http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
#include "Time_advance.h"

#ifndef _C_K_CHEBY
#define _C_K_CHEBY
inline int c_k(int k) { return (k==0) ? 2 : 1; }
inline int beta_k(int k) {return (k >= Nx-2) ? 0 : 1;}
#endif


//******************************************************************************

void Time_advance_incompress::Helmholtz_real(Array<Real,1> f, Array<Real,1> x, Real lambda, Real value_plus1, Real value_minus1)
{
    static Array<Real,1> even_x(Nx/2);
    static Array<Real,1> odd_x(Nx/2);
    static Array<Real,1> half_f(Nx/2);

    // Even part
    for (int i=0; i<Nx/2; i++)
        half_f(i) = f(2*i);
    
    Tridiagonal_solver(half_f, even_x, lambda, (value_plus1+value_minus1)/2, 0);
    
    // Odd part
    for (int i=0; i<Nx/2; i++)
        half_f(i) = f(2*i+1);
    
    Tridiagonal_solver(half_f, odd_x, lambda, (value_plus1-value_minus1)/2, 1);
    
    for (int i=0; i<Nx/2; i++) {
        x(2*i) = even_x(i);
        x(2*i+1) = odd_x(i);
    }
}

//******************************************************************************

// Solves A u = f with the top entry
// The matrix entries: [1,1,1..; a1,b1,c1, ..;  ] = [d0, d1, d2..]
void Time_advance_incompress::Helmholtz_complex(Array<Complex,1> f, Array<Complex, 1> u, Real lambda, Complex value_plus1, Complex value_minus1)
{
    static Array<Real,1> f_part(Nx);
    static Array<Real,1> u_part(Nx);
    
    f_part = real(f);
    Helmholtz_real(f_part, u_part, lambda, real(value_plus1), real(value_minus1));
    real(u) = u_part;
    
    f_part = imag(f);
    Helmholtz_real(f_part, u_part, lambda, imag(value_plus1), imag(value_minus1));
    imag(u) = u_part;
}

//******************************************************************************

// lambda = k_perp^2 + lambda_supplement
void Time_advance_incompress::Helmholtz_real_full_array(Array<Real,3> R, Array<Real,3> u, bool pressure_switch, Real lambda_supplement,  Real value_plus1, Real value_minus1)
{
    Real lambda, Kperpsqr;
    int ky,kz;
    
	for(int ly=0; ly<local_Ny; ly++)
		for (int lz=0; lz<local_Nz; lz++) {
			global.temp_array.in_helm_real = R(ly,lz,Range::all());
			
			ky = universal->Get_ky(ly);
			kz = universal->Get_kz(lz);
			
			Kperpsqr = pow2(ky*kfactor[2]) + pow2(kz*kfactor[3]);
        
            // k_perp = 0 for pressure yields del^2 p = 0 whose soln is 
            if ((Kperpsqr < MYEPS) && (pressure_switch)) {
                u(ly,lz,0) = (value_plus1+value_minus1)/2;
                u(ly,lz,1) = (value_plus1-value_minus1)/2;
                u(ly,lz,Range(2,Nx-1)) = 0.0;
            }
                
            else {
                lambda = Kperpsqr+lambda_supplement;
                
                Helmholtz_real(global.temp_array.in_helm_real, global.temp_array.out_helm_real, lambda, value_plus1, value_minus1);
                
                u(ly,lz,Range::all()) = global.temp_array.out_helm_real;
            }
        }
}

//******************************************************************************

// lambda = k_perp^2 + lambda_supplement
void Time_advance_incompress::Helmholtz_complex_full_array(Array<Complex,3> R, Array<Complex,3> u, Real lambda_supplement, Complex value_plus1, Complex value_minus1)
{
    Real lambda, Kperpsqr;
    int ky,kz;
    
	if ((basis_type == "ChFF") || (basis_type=="ChSF")) {
		for(int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				global.temp_array.in_helm_complex = R(ly,lz,Range::all());
				
				ky = universal->Get_ky(ly);
				kz = universal->Get_kz(lz);
				
				Kperpsqr = pow2(ky*kfactor[2]) + pow2(kz*kfactor[3]);
				
				lambda = Kperpsqr+lambda_supplement;
				
				Helmholtz_complex(global.temp_array.in_helm_complex, global.temp_array.out_helm_complex, lambda, value_plus1, value_minus1);	
				
				u(ly,lz,Range::all()) = global.temp_array.out_helm_complex;
			}
	}
	
	else if (basis_type == "ChSS") {
		for(int ly=0; ly<local_Ny; ly++)
			for (int lz=0; lz<local_Nz; lz++) {
				ky = universal->Get_ky(ly);
				kz = universal->Get_kz(lz);
				
				global.temp_array.in_helm_real = real(R(ly,lz,Range::all()));
				
				Kperpsqr = pow2(ky*kfactor[2]) + pow2(kz*kfactor[3]);
				
				lambda = Kperpsqr+lambda_supplement;
				
				Helmholtz_real(global.temp_array.in_helm_real, global.temp_array.out_helm_real, lambda, real(value_plus1), real(value_minus1));
				
				real(u(ly,lz,Range::all())) = global.temp_array.out_helm_real;
				
				// for imag
				global.temp_array.in_helm_real = imag(R(ly,lz,Range::all()));
				
				Kperpsqr = pow2(ky*kfactor[2]) + pow2((kz+1)*kfactor[3]);
				
				lambda = Kperpsqr+lambda_supplement;
				
				Helmholtz_real(global.temp_array.in_helm_real, global.temp_array.out_helm_real, lambda, real(value_plus1), real(value_minus1));
				
				imag(u(ly,lz,Range::all())) = global.temp_array.out_helm_real;
				
			}
	}
		
}


//******************************************************************************

// d=rhs(Nx/2)
void Time_advance_incompress::Tridiagonal_solver(Array<Real,1> f, Array<Real,1> y, Real lambda, Real d0, int odd_switch)
{
	int k;
    
    // The matrix entries: [1,1,1..; a1,b1,c1, ..;  ] = [d0, d1, d2..]
    static Array<Real,1> a(Nx/2);  // m=Nx/2;
    static Array<Real,1> b(Nx/2); 
    static Array<Real,1> c(Nx/2); 
    static Array<Real,1> d(Nx/2);
    
    
    a(0)=0; // not defined really
    b(0)=1; c(0)=1;
    d(0)=d0;
    
		// Nx = 2*Z, u(0:Nx-1), even k: 0, 2, .. 2*(Z-1); odd k: 1, 2, ...2*Z-1
		// if odd_switch=0: a,b,c,d(k) for k=2,4, ...N
		// if odd_switch=1: a,b,c,d(k) for k=1,4, ...N
    for (int i=1; i<Nx/2; i++) {
        k = 2*i+odd_switch;
        a(i) = c_k(k-2)*lambda/(4*k*(k-1));
        b(i) = -(1 + lambda*beta_k(k)/(2*(k*k-1)));
        c(i) = lambda*beta_k(k+2)/(4*k*(k+1));
        d(i) = c_k(k-2)*f(i-1)/(4*k*(k-1)) -beta_k(k)*f(i)/(2*(k*k-1)) +beta_k(k+2)*f(i+1)/(4*k*(k+1));
    }

    
    // Thual algorithm follows (Peret, Appendix B)
    //To do: Apply static to various arrays

    int n = Nx/2-1;

    Array<Real,1> X(n);
    Array<Real,1> Y(n);
   

    //STEP 1
    
    X(n-1) = - a(n)/b(n);
    Y(n-1) =  d(n)/b(n);

    for (int i=n-1; i>0; i--){
        X(i-1) = -a(i)/(b(i) + c(i) * X(i));
        Y(i-1) = ( d(i) - (c(i) * Y(i)) ) / ( b(i) + c(i) * X(i) );
    }

    //STEP 2
    Array<Real,1> theta(n+1);
    Array<Real,1> lamda(n+1);
    
    theta(0) = 1;
    lamda(0) = 0;

    for (int i=1; i<=n; i++){
        theta(i) = X(i-1) * theta(i-1);
        lamda(i) = ( X(i-1) * lamda(i-1) ) + Y(i-1);
    }

    //STEP 3
    Real THETA  = 0;
    Real LAMDA = 0;
    for (int i=0; i<=n; i++){
        THETA += theta(i);
        LAMDA += lamda(i);
    }

    //#STEP 4
    y(0) = (d(0) - LAMDA) / THETA;
    for (int i=0; i<n; i++)
        y(i+1) = X(i) * y(i) + Y(i);
      
}




