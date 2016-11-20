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

/*! \file  compute_force_main.cc
 * 
 * @brief  Compute force and put it in F.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "Time_advance.h"


//*********************************************************************************************

void Time_advance_incompress::Time_advance_step(FluidVF& U, Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U, P, Force);
		
	else if (global.program.integration_scheme == "RK2")
		RK2(U, P, Force);
		
	else if (global.program.integration_scheme == "RK4")
		RK4(U, P, Force);	
	
	else if (global.program.integration_scheme == "BDF1")
		BDF1(U, P, Force);
	
	else if (global.program.integration_scheme == "Adam_Bashforth")
		BDF1(U, P, Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true)
        U.Satisfy_strong_reality_condition_field();
   
    if (global.program.apply_weak_realitycond_alltime_switch == true)
        U.Satisfy_weak_reality_condition_field();
}


void Time_advance_incompress::Time_advance_step(FluidVF& U, FluidSF& T, Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U, T, P, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(U, T, P, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(U, T, P,Force);	
	
	else if (global.program.integration_scheme == "BDF1")
		BDF1(U, T, P,Force);
	
	else if (global.program.integration_scheme == "Adam_Bashforth")
		BDF1(U, T, P,Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true) {
        U.Satisfy_strong_reality_condition_field();
        T.Satisfy_strong_reality_condition_field();
    }
	
    if (global.program.apply_weak_realitycond_alltime_switch == true) {
        U.Satisfy_weak_reality_condition_field();
        T.Satisfy_weak_reality_condition_field();
    }
}

void Time_advance_incompress::Time_advance_step(FluidVF& U, FluidSF& T1, FluidSF& T2, Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U, T1, T2, P, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(U, T1, T2, P, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(U, T1, T2, P,Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true) {
        U.Satisfy_strong_reality_condition_field();
        T1.Satisfy_strong_reality_condition_field();
        T2.Satisfy_strong_reality_condition_field();
    }
    
    if (global.program.apply_weak_realitycond_alltime_switch == true) {
        U.Satisfy_weak_reality_condition_field();
        T1.Satisfy_weak_reality_condition_field();
        T2.Satisfy_weak_reality_condition_field();
    }
	
}


void Time_advance_incompress::Time_advance_step(FluidVF& U, FluidVF& W, Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U, W, P, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(U, W, P, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(U, W, P, Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true) {
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
    }
    
    if (global.program.apply_weak_realitycond_alltime_switch == true) {
        U.Satisfy_weak_reality_condition_field();
        W.Satisfy_weak_reality_condition_field();
    }
				
}


void Time_advance_incompress::Time_advance_step(FluidVF& U, FluidVF& W, FluidSF& T, Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U,  W, T, P, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(U,  W, T, P, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(U,  W, T, P, Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true) {
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
        T.Satisfy_strong_reality_condition_field();
    }
    
    if (global.program.apply_weak_realitycond_alltime_switch == true) {
        U.Satisfy_weak_reality_condition_field();
        W.Satisfy_weak_reality_condition_field();
        T.Satisfy_weak_reality_condition_field();
    }
}

void Time_advance_incompress::Time_advance_step(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C,Pressure& P, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(U,  W, T, C, P, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(U,  W, T, C, P, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(U,  W, T, C, P, Force);
    
    if (global.program.apply_strong_realitycond_alltime_switch == true) {
        U.Satisfy_strong_reality_condition_field();
        W.Satisfy_strong_reality_condition_field();
        T.Satisfy_strong_reality_condition_field();
		C.Satisfy_strong_reality_condition_field();
    }
    
    if (global.program.apply_weak_realitycond_alltime_switch == true) {
        U.Satisfy_weak_reality_condition_field();
        W.Satisfy_weak_reality_condition_field();
        T.Satisfy_weak_reality_condition_field();
		C.Satisfy_weak_reality_condition_field();
    }
}



void Time_advance_incompress::Time_advance_step(FluidSF& T, FORCE& Force)
{
	if (global.program.integration_scheme == "EULER")
		Euler(T, Force);
	
	else if (global.program.integration_scheme == "RK2")
		RK2(T, Force);
	
	else if (global.program.integration_scheme == "RK4")
		RK4(T, Force);
    
    // No reality condition for GP since wavefn is complex.
  /*  if (global.program.apply_strong_realitycond_alltime_switch == true)
        T.Satisfy_strong_reality_condition_field();
    
    if (global.program.apply_weak_realitycond_alltime_switch == true) 
        T.Satisfy_weak_reality_condition_field(); */
}

//*******************************  End of compute_force_main.cc *******************************



