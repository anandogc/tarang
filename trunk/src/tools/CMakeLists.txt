
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


########################################
#tarang_field_version_convert
########################################
ADD_EXECUTABLE(tarang_field_version_convert
	tarang_field_version_convert
	)
	
TARGET_LINK_LIBRARIES(tarang_field_version_convert basis_basicfn global ${SYSTEM_LIBRARIES})

ADD_CUSTOM_COMMAND(TARGET tarang_field_version_convert
                   POST_BUILD
                   COMMAND mv -f tarang_field_version_convert ${PROJECT_BINARY_DIR}
                   COMMENT "moving tarang_field_version_convert executable to PROJECT_BINARY_DIR" )
########################################

########################################
#TARANG ASCII to HDF5
########################################
ADD_EXECUTABLE(tarang_ascii_to_hdf5
	tarang_ascii_to_hdf5
	)
	
TARGET_LINK_LIBRARIES(tarang_ascii_to_hdf5 basis_basicfn global ${SYSTEM_LIBRARIES})

ADD_CUSTOM_COMMAND(TARGET tarang_ascii_to_hdf5
                   POST_BUILD
                   COMMAND mv -f tarang_ascii_to_hdf5 ${PROJECT_BINARY_DIR}
                   COMMENT "moving tarang_ascii_to_hdf5 executable to PROJECT_BINARY_DIR" )
########################################
