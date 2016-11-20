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
class Field_k(object):
	def __init__(self, m1, m2, m3):
		self.mode1 = m1
		self.mode2 = m2
		self.mode3 = m3
		self.Ux_real = []
		self.Ux_imag = []
		self.Uy_real = []
		self.Uy_imag = []
		self.Uz_real = []
		self.Uz_imag = []
		self.T_real = []
		self.T_imag = []
		self.Bx_real = []
		self.Bx_imag = []
		self.By_real = []
		self.By_imag = []
		self.Bz_real = []
		self.Bz_imag = []
		self.Time = []

	def check_mode(self, m1, m2, m3):
		if(self.mode1 == m1 and self.mode2 == m2 and self.mode3 == m3):
			return True
		else:
			return False

	def insert_Fluid(self, Ux_real, Ux_imag, Uy_real, Uy_imag, Uz_real, Uz_imag, Time):
		self.Ux_real.append(Ux_real)
		self.Ux_imag.append(Ux_imag)
		self.Uy_real.append(Uy_real)
		self.Uy_imag.append(Uy_imag)
		self.Uz_real.append(Uz_real)
		self.Uz_imag.append(Uz_imag)
		self.Time.append(Time)
		
	def insert_RBC(self, Ux_real, Ux_imag, Uy_real, Uy_imag, Uz_real, Uz_imag, T_real, T_imag, Time):
		self.Ux_real.append(Ux_real)
		self.Ux_imag.append(Ux_imag)
		self.Uy_real.append(Uy_real)
		self.Uy_imag.append(Uy_imag)
		self.Uz_real.append(Uz_real)
		self.Uz_imag.append(Uz_imag)
		self.T_real.append(T_real)
		self.T_imag.append(T_imag)
		self.Time.append(Time)
		
	def insert_MHD(self, Ux_real, Ux_imag, Uy_real, Uy_imag, Uz_real, Uz_imag, Bx_real, Bx_imag, By_real, By_imag, Bz_real, Bz_imag, Time):
		self.Ux_real.append(Ux_real)
		self.Ux_imag.append(Ux_imag)
		self.Uy_real.append(Uy_real)
		self.Uy_imag.append(Uy_imag)
		self.Uz_real.append(Uz_real)
		self.Uz_imag.append(Uz_imag)
		self.Bx_real.append(Bx_real)
		self.Bx_imag.append(Bx_imag)
		self.By_real.append(By_real)
		self.By_imag.append(By_imag)		
		self.Bz_real.append(Bz_real)
		self.Bz_imag.append(Bz_imag)		
		self.Time.append(Time)
