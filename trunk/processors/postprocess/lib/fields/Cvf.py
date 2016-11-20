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

class Cvf(object):
	def __init__(self, flag, kind):
		if kind == "FLUID_INCOMPRESS":
			############## Data of glob.d file ################
			self.total_energy = constants.data[:,1]
			self.total_dissipation = constants.data[:,2]
			self.total_helicity1 = constants.data[:,3]
			self.total_helicity2 = constants.data[:,4]
			self.total_k2H1 = constants.data[:,5]
			self.total_k2H2 = constants.data[:,6]
			self.entropy = constants.data[:,7]
			self.kmax_eta = constants.data[:,8]
			self.dissipation_coefficient = constants.data[:,9]
			self.total_E1 = constants.data[:,12]
			self.total_E2 = constants.data[:,13]
			self.total_E3 = constants.data[:,14]
			self.total_k2energy = None
			self.total_k2Hc = None
			############## Data of glob.d file ################

			############## Data of spectrum.d file ############
			# each variable is a list here which contains elements according to time given
			# these variables are filled when user called read function for spectrum
			self.Uek_shell_ek1 = []
			self.Uek_shell_ek2 = []
			self.Uek_shell_ek3 = []
			self.UDk_shell_dissk1 = []
			self.UDk_shell_dissk2 = []
			self.UDk_shell_dissk3 = []
			self.Fv_v_shell_ek1 = []
			self.Fv_v_shell_ek2 = []
			self.Fv_v_shell_ek3 = []
			############## Data of spectrum.d file ############

			############## Data of field_k_out.d file #########
			self.Modek = []
			############## Data of field_k_out.d file #########

			############## Data of field_r_out.d file #########
			self.Moder = []
			############## Data of field_r_out.d file #########

		elif kind == "RBC":
			############## Data of glob.d file ################
			self.total_energy = constants.data[:,1]
			self.total_dissipation = constants.data[:,3]
			self.total_helicity1 = constants.data[:,6]
			self.total_helicity2 = constants.data[:,7]
			self.total_k2H1 = constants.data[:,8]
			self.total_k2H2 = constants.data[:,9]
			self.entropy = constants.data[:,10]
			self.kmax_eta = constants.data[:,12]
			self.dissipation_coefficient = constants.data[:,14]
			self.total_E1 = constants.data[:,18]
			self.total_E2 = constants.data[:,19]
			self.total_E3 = constants.data[:,20]
			self.total_k2energy = None
			self.total_k2Hc = None
			############## Data of glob.d file ################

			############## Data of spectrum.d file ############
			# each variable is a list here which contains elements according to time given
			# these variables are filled when user called read function for spectrum
			self.Uek_shell_ek1 = []
			self.Uek_shell_ek2 = []
			self.Uek_shell_ek3 = []
			self.UDk_shell_dissk1 = []
			self.UDk_shell_dissk2 = []
			self.UDk_shell_dissk3 = []
			self.Fv_v_shell_ek1 = []
			self.Fv_v_shell_ek2 = []
			self.Fv_v_shell_ek3 = []
			############## Data of spectrum.d file ############

			############## Data of field_k_out.d file #########
			self.Modek = []
			############## Data of field_k_out.d file #########

			############## Data of field_r_out.d file #########
			self.Moder = []
			############## Data of field_r_out.d file #########
		elif kind == "MHD_INCOMPRESS":
			############## Data of glob.d file ################
			if flag == "U":
				self.total_energy = constants.data[:,1]
				self.total_dissipation = constants.data[:,3]
				self.total_helicity1 = constants.data[:,7]
				self.total_helicity2 = constants.data[:,8]
				self.total_k2H1 = constants.data[:,11]
				self.total_k2H2 = constants.data[:,12]
				self.entropy = constants.data[:,15]
				self.kmax_eta = constants.data[:,17]
				self.dissipation_coefficient = constants.data[:,19]
				self.total_E1 = constants.data[:,23]
				self.total_E2 = constants.data[:,24]
				self.total_E3 = constants.data[:,25]
				self.total_k2energy = None
				self.total_k2Hc = None
			elif flag == "B":
				self.total_energy = constants.data[:,2]
				self.total_dissipation = constants.data[:,4]
				self.total_helicity1 = constants.data[:,9]
				self.total_helicity2 = constants.data[:,10]
				self.total_k2H1 = constants.data[:,13]
				self.total_k2H2 = constants.data[:,14]
				self.entropy = constants.data[:,16]
				self.kmax_eta = constants.data[:,18]
				self.dissipation_coefficient = constants.data[:,20]
				self.total_E1 = constants.data[:,26]
				self.total_E2 = constants.data[:,27]
				self.total_E3 = constants.data[:,28]
				self.total_k2energy = None
				self.total_k2Hc = None
			############## Data of glob.d file ################

			############## Data of spectrum.d file ############
			# each variable is a list here which contains elements according to time given
			# these variables are filled when user called read function for spectrum
			self.Uek_shell_ek1 = []
			self.Uek_shell_ek2 = []
			self.Uek_shell_ek3 = []
			self.UDk_shell_dissk1 = []
			self.UDk_shell_dissk2 = []
			self.UDk_shell_dissk3 = []
			self.Fv_v_shell_ek1 = []
			self.Fv_v_shell_ek2 = []
			self.Fv_v_shell_ek3 = []
			############## Data of spectrum.d file ############

			############## Data of field_k_out.d file #########
			self.Modek = []
			############## Data of field_k_out.d file #########

			############## Data of field_r_out.d file #########
			self.Moder = []
			############## Data of field_r_out.d file #########


