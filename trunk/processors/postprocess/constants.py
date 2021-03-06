"""
 *Tarang Post processor
 *
 * Copyright (C) 2014  Mahendra K. Verma (mkv@iitk.ac.in)
 *
 * @Author: Vijay Jain, Abhishek Kumar
 *
 * 
 *
 * Tarang Post processor is free software; you can redistribute it and/or
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
 *
"""
import numpy
import linecache

directory = ''
data = ''
spectrum_timelist = []
shell_to_shell_timelist = []
ring_spectrum_timelist = []
flux_timelist = []
shell_data = []
K = []
KF = []
