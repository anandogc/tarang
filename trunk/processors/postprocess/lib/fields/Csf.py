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
import constants

class Csf(object):
	def __init__(self, flag, kind):
		if kind == "RBC":
			################ data of glob.d ########################
			self.total_energy = constants.data[:,2]
			self.total_dissipation = constants.data[:,4]			
			self.total_k2energy = None
			self.entropy = constants.data[:,11]
			self.kmax_eta = constants.data[:,13]
			self.diffusion_coefficient = constants.data[:,15]
			################ data of glob.d ########################
			################ data of spectrum.d ####################
			self.Tek_shell_ek = []
			self.FT_T_shell_ek = []
			################ data of spectrum.d ####################
