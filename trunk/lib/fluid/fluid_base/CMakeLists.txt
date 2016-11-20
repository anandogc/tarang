
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
ADD_LIBRARY(fluid_base
	Correlation
	FluidSF
	FluidVF
	PlainFluidVF
	MHD

    io/FluidIO.cc
    io/input/ic_amp_phase.cc
    io/input/ic_model.cc
    io/input/init_cond_ABC.cc
    io/input/init_cond_dynamo_6mode.cc
    io/input/init_cond_dynamo_velocity_field.cc
    io/input/init_cond_energy.cc
    io/input/init_cond_field.cc
    io/input/init_cond_main.cc
    io/input/init_cond_modes.cc
    io/input/init_cond_RayleighTaylor.cc
    io/input/init_cond_TG.cc
    io/input/init_cond_channel.cc
    io/input/init_cond_vortex.cc
    io/output/output_energy.cc
    io/output/output_energy_real.cc
    io/output/output_field.cc
    #IO/output_str_fn.cc

    force/compute_force_ABC.cc
    force/compute_force_astro.cc
    force/compute_force_Coriolis.cc
    force/compute_force_decay.cc
    force/compute_force_dynamo_6mode.cc
    force/compute_force_ek_hk.cc
    force/compute_force_LM.cc
    force/compute_force_main.cc
    force/compute_force_modes.cc
    force/compute_force_MRBC.cc
    force/compute_force_random.cc
    force/compute_force_RBC.cc
    force/compute_force_TG.cc
    force/force_alpha_beta.cc
    force/force_amp_phase.cc
    force/force_model.cc
    force/FORCE.h
)
