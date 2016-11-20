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
class Field_r(object):
	def __init__(self, mode1, mode2, mode3):
		self.mode1 = mode1
		self.mode2 = mode2
		self.mode3 = mode3
		self.Ux = []
		self.Uy = []
		self.Uz = []
		self.Bx = []
		self.By = []
		self.Bz = []
		self.T = []
		self.Time = []

	def check_mode(self, m1, m2, m3):
		if(self.mode1 == m1 and self.mode2 == m2 and self.mode3 == m3):
			return True
		else:
			return False

	def insert_Fluid(self, Ux, Uy, Uz, Time):
		self.Ux.append(Ux)
		self.Uy.append(Uy)
		self.Uz.append(Uz)
		self.Time.append(Time)
	
	def insert_RBC(self, Ux, Uy, Uz, T, Time):
		self.Ux.append(Ux)
		self.Uy.append(Uy)
		self.Uz.append(Uz)
		self.T.append(T)
		self.Time.append(Time)
	
	def insert_MHD(self, Ux, Uy, Uz, Bx, By, Bz, Time):
		self.Ux.append(Ux)
		self.Uy.append(Uy)
		self.Uz.append(Uz)
		self.Bx.append(Bx)
		self.By.append(By)
		self.Bz.append(Bz)
		self.Time.append(Time)
