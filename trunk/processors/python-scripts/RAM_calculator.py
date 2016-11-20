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

 # \file  CMakeLists.txt
 # @author  Abhishek Kumar, A. G. Chatterjee, M. K. Verma
 # @date 2 April 2013 

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np


array = 18 #Fluid = 15, RBC = 18, MHD = 27
pre = 8 #double = 8, float = 4
n = 512
integration_schme = 4 #Euler = 1,  RK2 = 2, RK4 = 4

ram_2D = integration_schme*array*pre*2*(n**2)/(1024*1024*1024)
ram_3D = integration_schme*array*pre*(n**3)/(1024*1024*1024)

print "Required RAM for 3D in (GB) for Euler = ", ram_3D
print "Required RAM for 2D in (GB) for Euler = ", ram_2D



