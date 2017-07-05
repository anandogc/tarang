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


/*! \file  init_cond_TG.cc
 * 
 * @brief Initial conditions as Taylor Green flow (TG).
 *
 * @note The parameters are read from parameter file. 
 *
 *		Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x)sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "FluidIO.h"


//********************************************************************************************

void FluidIO::Init_cond_user_defined1(FluidVF& U)
{
}

void FluidIO::Init_cond_user_defined2(FluidVF& U)
{
}

//******************************************************************************

void FluidIO::Init_cond_user_defined1(FluidVF& U, FluidSF& T)
{
}

void FluidIO::Init_cond_user_defined2(FluidVF& U, FluidSF& T)
{
}

//******************************************************************************

void FluidIO::Init_cond_user_defined1(FluidVF& U, FluidVF& W)
{
}

void FluidIO::Init_cond_user_defined2(FluidVF& U, FluidVF& W)
{
}

//******************************************************************************







