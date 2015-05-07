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

/*! \file  Global.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "Global.h"

//overload << to print std::vector nicely
template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(os, ", "));
    os << v.back();
    os << "]";
    return os;
}

template<class T>
string Print_linear_blitz_array(const Array<T,1> A)
{
    string str = "[";
    if (A.size()>0)
    {
        for (int i=0; i<A.size()-1; i++)
            str+=To_string(A(i))+", ";
        str+=To_string(A(A.size()-1));
    }
    str+="]";
    return str;
}

template<class T>
string Print_vector(const vector<T> vec, int start_index=0)
{
    string str = "[";
    if (vec.size()>start_index)
    {
        for (int i=start_index; i<vec.size()-1; i++)
            str+=To_string(vec[i])+", ";
        str+=To_string(vec[vec.size()-1]);
    }
    str+="]";
    return str;
}

string Indent(int level)
{
    int Indent_per_level=4;
    return string(Indent_per_level*level,' ');
}

void Global::Print()
{
	if (mpi.master){
		program.Print(0);
		PHYSICS.Print(0);
		field.Print(0);
		time.Print(0);
		force.Print(0);
		io.Print(0);
		energy_transfer.Print(0);
		mpi.Print(0);
		spectrum.Print(0);
		cout<<endl;
	}
}

void Global::program::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"program:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"version: " << version
        <<'\n'<<Indent(my_level+1)<<"kind: " << kind
        <<'\n'<<Indent(my_level+1)<<"iter_or_diag: " << iter_or_diag
		<<'\n'<<Indent(my_level+1)<<"decomposition: " << decomposition
		<<'\n'<<Indent(my_level+1)<<"basis_type: " << basis_type
        <<'\n'<<Indent(my_level+1)<<"alias_option: " << alias_option
        <<'\n'<<Indent(my_level+1)<<"integration_scheme: " << integration_scheme
        <<'\n'<<Indent(my_level+1)<<"LES_switch: " << LES_switch
        <<'\n'<<Indent(my_level+1)<<"T_exists: " << T_exists
        <<'\n'<<Indent(my_level+1)<<"W_exists: " << W_exists
        <<'\n'<<Indent(my_level+1)<<"apply_strong_realitycond_alltime_switch: " << apply_strong_realitycond_alltime_switch
        <<'\n'<<Indent(my_level+1)<<"apply_weak_realitycond_alltime_switch: " << apply_weak_realitycond_alltime_switch
        <<'\n'<<Indent(my_level+1)<<"low_dimensional_switch: " << low_dimensional_switch
        <<'\n'<<Indent(my_level+1)<<"two_and_half_dimension: " << two_and_half_dimension
        <<'\n'<<Indent(my_level+1)<<"two_dimension: " << two_dimension
        <<'\n'<<Indent(my_level+1)<<"dt_option: " << dt_option
        <<'\n'<<Indent(my_level+1)<<"helicity_switch: " << helicity_switch
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch: " << sincostr_switch
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_Vx: " << sincostr_switch_Vx
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_Vy: " << sincostr_switch_Vy
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_Vz: " << sincostr_switch_Vz
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_F: " << sincostr_switch_F
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_Visqr: " << sincostr_switch_Visqr
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_VxVy: " << sincostr_switch_VxVy
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_VxVz: " << sincostr_switch_VxVz
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_VyVz: " << sincostr_switch_VyVz
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_FVx: " << sincostr_switch_FVx
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_FVy: " << sincostr_switch_FVy
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_FVz: " << sincostr_switch_FVz
		<<'\n'<<Indent(my_level+1)<<"sincostr_switch_divergence " << sincostr_switch_divergence
        <<endl;

}

void Global::field::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"field:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"N: " << "["<< Nx <<", "<< Ny <<", "<< Nz <<"]"
        <<'\n'<<Indent(my_level+1)<<"kfactor: " << "["<<kfactor[1]<<", "<<kfactor[2]<<", "<<kfactor[3]<<"]"
        <<'\n'<<Indent(my_level+1)<<"L: " << "["<<L[1]<<", "<<L[2]<<", "<<L[3]<<"]"

        <<'\n'<<Indent(my_level+1)<<"diss_coefficients: "<<diss_coefficients
        <<'\n'<<Indent(my_level+1)<<"hyper_diss_coefficients: "<<hyper_diss_coefficients
        <<'\n'<<Indent(my_level+1)<<"hyper_diss_exponents: "<<hyper_diss_exponents
        <<endl;
}

void Global::time::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"time:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"init: "<<init
        <<'\n'<<Indent(my_level+1)<<"final: "<<final
        <<'\n'<<Indent(my_level+1)<<"dt_fixed: "<<dt_fixed
		<<'\n'<<Indent(my_level+1)<<"Courant_no: "<<Courant_no
		<<'\n'<<Indent(my_level+1)<<"job_time: "<< job_time
        <<endl;
}

void Global::force::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"force:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"U_switch: "<<U_switch
        <<'\n'<<Indent(my_level+1)<<"W_switch: "<<W_switch
        <<'\n'<<Indent(my_level+1)<<"T_switch: "<<T_switch
		<<'\n'<<Indent(my_level+1)<<"C_switch: "<<C_switch
        <<'\n'<<Indent(my_level+1)<<"field_procedure: "<<field_procedure
        <<'\n'<<Indent(my_level+1)<<"int_para : "<<Print_linear_blitz_array<int>(int_para)
        <<'\n'<<Indent(my_level+1)<<"double_para : "<<Print_linear_blitz_array<Real>(double_para)
        <<'\n'<<Indent(my_level+1)<<"string_para : "<<Print_linear_blitz_array<string>(string_para)
        <<endl;
}


//io

void Global::io::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"io:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"data_dir: "<<data_dir
        <<'\n'<<Indent(my_level+1)<<"input_field_procedure: "<<input_field_procedure
        <<'\n'<<Indent(my_level+1)<<"input_vx_vy_switch: "<<input_vx_vy_switch
        <<'\n'<<Indent(my_level+1)<<"output_vx_vy_switch: "<<output_vx_vy_switch
        <<'\n'<<Indent(my_level+1)<<"output_precision: "<<output_precision
	    <<'\n'<<Indent(my_level+1)<<"diagnostic_procedures: "<<Print_vector(diagnostic_procedures)
		<<'\n'<<Indent(my_level+1)<<"N_in_reduced: "<<Print_vector(N_in_reduced)
        <<'\n'<<Indent(my_level+1)<<"N_out_reduced: "<<Print_vector(N_out_reduced)
        <<'\n'<<Indent(my_level+1)<<"int_para : "<<Print_linear_blitz_array<int>(int_para)
        <<'\n'<<Indent(my_level+1)<<"double_para : "<<Print_linear_blitz_array<Real>(double_para)
        <<'\n'<<Indent(my_level+1)<<"string_para : "<<Print_linear_blitz_array<string>(string_para);

        probes.Print(my_level+1);
        init_cond_modes.Print(my_level+1);
        time.Print(my_level+1);
        cout<<endl;
}

void Global::io::probes::spectral_space::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"spectral_space:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"number: "<<number
        <<'\n'<<Indent(my_level+1)<<"coords: ";

    for (int i=0; i<number; i++)
        cout<<'\n'<<Indent(my_level+2)<<"- coord: "<<Print_linear_blitz_array<int>(coords(i, Range(1,3)));
    cout<<endl;
}

void Global::io::probes::real_space::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"real_space:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"number: "<<number
        <<'\n'<<Indent(my_level+1)<<"coords: ";

    for (int i=0; i<number; i++)
        cout<<'\n'<<Indent(my_level+2)<<"- coord: "<<Print_linear_blitz_array<int>(coords(i, Range(1,3)));
     cout<<endl;
}


void Global::io::probes::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"probes:";
    
    spectral_space.Print(my_level+1);
    real_space.Print(my_level+1);
    cout<<endl;
}

void Global::io::init_cond_modes::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"init_cond_modes:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"number: "<<number
        <<'\n'<<Indent(my_level+1)<<"number_components: "<<number_components
        <<'\n'<<Indent(my_level+1)<<"modes:";

    for (int i=0; i<number; i++){
        cout <<'\n'<<Indent(my_level+2)<<"- coord: "<< Print_linear_blitz_array<int>(coords(i, Range(1,3)));
        if (field_array_real.size()>0)
			 cout<<'\n'<<Indent(my_level+2)<<"- mode: "<< Print_linear_blitz_array<Real>(field_array_real(i, Range::all()));
		else if (field_array_complex.size()>0)
			 cout<<'\n'<<Indent(my_level+2)<<"- mode: "<< Print_linear_blitz_array<Complex>(field_array_complex(i, Range::all()));
	}
    cout<<endl;
}

void Global::io::time::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"time:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"global_save_next: "<<global_save_next
        <<'\n'<<Indent(my_level+1)<<"complex_field_save_next: "<<complex_field_save_next
        <<'\n'<<Indent(my_level+1)<<"field_frequent_save_next: "<<field_frequent_save_next
        <<'\n'<<Indent(my_level+1)<<"field_reduced_save_next: "<<field_reduced_save_next
        <<'\n'<<Indent(my_level+1)<<"real_field_save_next: "<<real_field_save_next
        <<'\n'<<Indent(my_level+1)<<"field_k_save_next: "<<field_k_save_next
        <<'\n'<<Indent(my_level+1)<<"field_r_save_next: "<<field_r_save_next
        <<'\n'<<Indent(my_level+1)<<"spectrum_save_next: "<<spectrum_save_next
        <<'\n'<<Indent(my_level+1)<<"pressure_save_next: "<<pressure_save_next
        <<'\n'<<Indent(my_level+1)<<"pressure_spectrum_save_next: "<<pressure_spectrum_save_next
        <<'\n'<<Indent(my_level+1)<<"flux_save_next: "<<flux_save_next
        <<'\n'<<Indent(my_level+1)<<"shell_to_shell_save_next: "<<shell_to_shell_save_next
        <<'\n'<<Indent(my_level+1)<<"ring_spectrum_save_next: "<<ring_spectrum_save_next
        <<'\n'<<Indent(my_level+1)<<"ring_to_ring_save_next: "<<ring_to_ring_save_next
        <<'\n'<<Indent(my_level+1)<<"cylindrical_ring_spectrum_save_next: "<<cylindrical_ring_spectrum_save_next
        <<'\n'<<Indent(my_level+1)<<"cylindrical_ring_to_ring_save_next: "<<cylindrical_ring_to_ring_save_next
        <<'\n'<<Indent(my_level+1)<<"structure_fn_save_next: "<<structure_fn_save_next
        <<'\n'<<Indent(my_level+1)<<"cout_save_next: "<<cout_save_next
        <<endl
        <<'\n'<<Indent(my_level+1)<<"global_save_interval: "<<global_save_interval
        <<'\n'<<Indent(my_level+1)<<"complex_field_save_interval: "<<complex_field_save_interval
        <<'\n'<<Indent(my_level+1)<<"field_frequent_save_interval: "<<field_frequent_save_interval
        <<'\n'<<Indent(my_level+1)<<"field_reduced_save_interval: "<<field_reduced_save_interval
        <<'\n'<<Indent(my_level+1)<<"real_field_save_interval: "<<real_field_save_interval
        <<'\n'<<Indent(my_level+1)<<"field_k_save_interval: "<<field_k_save_interval
        <<'\n'<<Indent(my_level+1)<<"field_r_save_interval: "<<field_r_save_interval
        <<'\n'<<Indent(my_level+1)<<"spectrum_save_interval: "<<spectrum_save_interval
        <<'\n'<<Indent(my_level+1)<<"pressure_save_interval: "<<pressure_save_interval
        <<'\n'<<Indent(my_level+1)<<"pressure_spectrum_save_interval: "<<pressure_spectrum_save_interval
        <<'\n'<<Indent(my_level+1)<<"flux_save_interval: "<<flux_save_interval
        <<'\n'<<Indent(my_level+1)<<"shell_to_shell_save_interval: "<<shell_to_shell_save_interval
        <<'\n'<<Indent(my_level+1)<<"ring_spectrum_save_interval: "<<ring_spectrum_save_interval
        <<'\n'<<Indent(my_level+1)<<"ring_to_ring_save_interval: "<<ring_to_ring_save_interval
        <<'\n'<<Indent(my_level+1)<<"cylindrical_ring_spectrum_save_interval: "<<cylindrical_ring_spectrum_save_interval
        <<'\n'<<Indent(my_level+1)<<"cylindrical_ring_to_ring_save_interval: "<<cylindrical_ring_to_ring_save_interval
        <<'\n'<<Indent(my_level+1)<<"structure_fn_save_interval: "<<structure_fn_save_interval
        <<'\n'<<Indent(my_level+1)<<"cout_save_interval: "<<cout_save_interval
        <<endl
		<<'\n'<<Indent(my_level+1)<<"global_save_last: "<<global_save_last
		<<'\n'<<Indent(my_level+1)<<"complex_field_save_last: "<<complex_field_save_last
		<<'\n'<<Indent(my_level+1)<<"field_frequent_save_last: "<<field_frequent_save_last
		<<'\n'<<Indent(my_level+1)<<"field_reduced_save_last: "<<field_reduced_save_last
		<<'\n'<<Indent(my_level+1)<<"real_field_save_last: "<<real_field_save_last
		<<'\n'<<Indent(my_level+1)<<"field_k_save_last: "<<field_k_save_last
		<<'\n'<<Indent(my_level+1)<<"field_r_save_last: "<<field_r_save_last
		<<'\n'<<Indent(my_level+1)<<"spectrum_save_last: "<<spectrum_save_last
		<<'\n'<<Indent(my_level+1)<<"pressure_save_last: "<<pressure_save_last
		<<'\n'<<Indent(my_level+1)<<"pressure_spectrum_save_last: "<<pressure_spectrum_save_last
		<<'\n'<<Indent(my_level+1)<<"flux_save_last: "<<flux_save_last
		<<'\n'<<Indent(my_level+1)<<"shell_to_shell_save_last: "<<shell_to_shell_save_last
		<<'\n'<<Indent(my_level+1)<<"ring_spectrum_save_last: "<<ring_spectrum_save_last
		<<'\n'<<Indent(my_level+1)<<"ring_to_ring_save_last: "<<ring_to_ring_save_last
		<<'\n'<<Indent(my_level+1)<<"cylindrical_ring_spectrum_save_last: "<<cylindrical_ring_spectrum_save_last
		<<'\n'<<Indent(my_level+1)<<"cylindrical_ring_to_ring_save_last: "<<cylindrical_ring_to_ring_save_last
		<<'\n'<<Indent(my_level+1)<<"structure_fn_save_last: "<<structure_fn_save_last
		<<'\n'<<Indent(my_level+1)<<"cout_save_last: "<<cout_save_last
		<<endl;
}

//energy_transfer
void Global::energy_transfer::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"energy_transfer:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"helicity_flux_switch: "<<helicity_flux_switch
        <<'\n'<<Indent(my_level+1)<<"helicity_shell_to_shell_switch: "<<helicity_shell_to_shell_switch
        <<'\n'<<Indent(my_level+1)<<"Elsasser: "<<Elsasser
        <<endl
        <<'\n'<<Indent(my_level+1)<<"Vpll_switch: "<<Vpll_switch <<endl;

    flux.Print(my_level+1);
    shell_to_shell.Print(my_level+1);
    ring_to_ring.Print(my_level+1);
    cylindrical_ring_to_ring.Print(my_level+1);
    cout<<endl;
}


void Global::energy_transfer::flux::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"flux:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_spheres: "<<no_spheres
        <<'\n'<<Indent(my_level+1)<<"radii: "<<Print_linear_blitz_array(radii)
        <<endl;
}

void Global::energy_transfer::shell_to_shell::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"shell_to_shell:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<'\n'<<Indent(my_level+1)<<"radii: "<<Print_linear_blitz_array(radii)
        <<endl;
}

void Global::energy_transfer::ring_to_ring::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"ring_to_ring:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<'\n'<<Indent(my_level+1)<<"no_sectors: "<<no_sectors
        <<'\n'<<Indent(my_level+1)<<"radii: "<<Print_linear_blitz_array<Real>(radii)
		<<'\n'<<Indent(my_level+1)<<"sector_option: "<<sector_option
        <<'\n'<<Indent(my_level+1)<<"sector_angles: "<<Print_linear_blitz_array<Real>(sector_angles)
        <<endl;
}

void Global::energy_transfer::cylindrical_ring_to_ring::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"cylindrical_ring_to_ring:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<'\n'<<Indent(my_level+1)<<"no_slabs: "<<no_slabs
        <<'\n'<<Indent(my_level+1)<<"radii: "<<Print_linear_blitz_array<Real>(radii)
        <<'\n'<<Indent(my_level+1)<<"kpll_option: "<<kpll_option
        <<'\n'<<Indent(my_level+1)<<"kpll_array: "<<Print_linear_blitz_array<Real>(kpll_array)
        <<endl;
}

void Global::mpi::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"mpi:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"numprocs: "<<numprocs
        <<'\n'<<Indent(my_level+1)<<"num_p_col: "<<num_p_cols
        <<'\n'<<Indent(my_level+1)<<"num_p_row: "<<num_p_rows
        <<endl; 
}

void Global::PHYSICS::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"PHYSICS:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"Pr_option: "<<Pr_option
        <<'\n'<<Indent(my_level+1)<<"Uscaling: "<<Uscaling
        <<'\n'<<Indent(my_level+1)<<"Rayleigh: "<<Rayleigh
        <<'\n'<<Indent(my_level+1)<<"Prandtl: "<<Prandtl
        <<'\n'<<Indent(my_level+1)<<"temperature_grad: "<<temperature_grad
		<<'\n'<<Indent(my_level+1)<<"Chandrasekhar: "<<Chandrasekhar
		<<'\n'<<Indent(my_level+1)<<"Prandtl_mag: "<<Prandtl_mag
		<<'\n'<<Indent(my_level+1)<<"Prandtl_c: "<<Prandtl_c
		<<'\n'<<Indent(my_level+1)<<"Reynolds: "<<Reynolds
		<<'\n'<<Indent(my_level+1)<<"Reynolds_mag: "<<Reynolds_mag
		<<'\n'<<Indent(my_level+1)<<"Peclet: "<< Peclet
		<<'\n'<<Indent(my_level+1)<<"Peclet_c: "<< Peclet_c
        <<endl; 
}

void Global::spectrum::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"spectrum:";
    shell.Print(my_level+1);
    ring.Print(my_level+1);
    cylindrical_ring.Print(my_level+1);
    cout<<endl;
}

void Global::spectrum::shell::Print(int my_level)
{
    cout<<'\n'<<Indent(my_level)<<"shell:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<endl;
}


void Global::spectrum::ring::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"ring:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<'\n'<<Indent(my_level+1)<<"no_sectors: "<<no_sectors
        <<'\n'<<Indent(my_level+1)<<"sector_option: "<<sector_option
        <<'\n'<<Indent(my_level+1)<<"sector_angles: "<<Print_linear_blitz_array<Real>(sector_angles)
        <<endl;
}

void Global::spectrum::cylindrical_ring::Print(int my_level)
{
    if (!turnon)
        return;

    cout<<'\n'<<Indent(my_level)<<"cylindrical_ring:";
    cout<<boolalpha
        <<'\n'<<Indent(my_level+1)<<"turnon: "<<turnon
        <<'\n'<<Indent(my_level+1)<<"no_shells: "<<no_shells
        <<'\n'<<Indent(my_level+1)<<"no_slabs: "<<no_slabs
        <<'\n'<<Indent(my_level+1)<<"kpll_option: "<<kpll_option
        <<'\n'<<Indent(my_level+1)<<"kpll_array: "<<Print_linear_blitz_array<Real>(kpll_array)
        <<endl;
}
