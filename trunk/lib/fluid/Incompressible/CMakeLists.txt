
 # Tarang-2
 # 
 # Copyright (C) 2008, 2009  Mahendra K. Verma
 #
 # Mahendra K. Verma
 # Indian Institute of Technology, Kanpur-208016
 # UP, India
 #
 # mkv@iitk.ac.in
 #
 # This file is part of Tarang-2 .
 #
 # Tarang-2 is free software; you can redistribute it and/or
 # modify it under the terms of the GNU General Public License
 # as published by the Free Software Foundation; either version 2
 # of the License, or (at your option) any later version.
 # Tarang-2 is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU General Public License
 # along with Tarang-2; if not, write to the Free Software
 # Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 #
#SET(CMAKE_CXX_COMPILER mpicxx)

# ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)

# ADD_DEFINITIONS(-DWAVENOACTUAL -DANISDIRN1)

# INCLUDE_DIRECTORIES(/Users/mkv/local/include )

# LINK_DIRECTORIES(/Users/mkv/local/lib )

ADD_LIBRARY(fluid_incompressible
	Pressure
	IncIO/IncIO
	IncIO/output/output_pressure.cc
	IncIO/output/output_main.cc
	IncIO/output/output_ET.cc
	nlin/compute_diag_offdiag.cc
	nlin/compute_nlin_ET_VF.cc
	nlin/compute_nlin_VF_MHD.cc
	nlin/compute_nlin_VF_2SF.cc
	nlin/Nlin.h
	nlin/compute_nlin.cc
	nlin/compute_nlin_GP.cc
	nlin/compute_nlin_SF.cc
	nlin/compute_nlin_2SF.cc
	nlin/compute_nlin_VF_SF.cc
	time_advance/RK.cc
	time_advance/RK_SF.cc
        time_advance/RK_GP.cc
	time_advance/RK_2SF.cc
	time_advance/RK_VF_MHD.cc
	time_advance/RK_VF_SF.cc
        time_advance/RK_VF_2SF.cc
	time_advance/Tadvance_functions.cc
	time_advance/Time_advance.cc
	time_advance/Time_advance.h
        time_advance/Chebyshev/Adam_Bashforth.cc
       time_advance/Chebyshev/Adam_Bashforth_SF.cc
	time_advance/Chebyshev/BDF1.cc
	time_advance/Chebyshev/BDF1_SF.cc
        time_advance/Chebyshev/helmholtz.cc
	time_advance/Chebyshev/influence_matrix.cc
	EnergyTr/compute_cylin_ring_tr.cc
       EnergyTr/compute_energy_tr_Vpll.cc
	EnergyTr/compute_flux.cc
	EnergyTr/compute_ring_tr.cc
	EnergyTr/compute_shell_tr.cc
	EnergyTr/energy_tr_B0.cc
	EnergyTr/EnergyTr.cc
	EnergyTr/EnergyTr.h
	EnergyTr/fill_space.cc
	EnergyTr/force_V.cc
	EnergyTr/prod_nlinV.cc
)
