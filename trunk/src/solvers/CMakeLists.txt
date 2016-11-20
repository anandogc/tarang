
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

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/src/solvers
)

ADD_EXECUTABLE(tarang
	main.cc
	GP/GP_main
	ifluid/Ifluid_main
	imhd/IMHD_main
	imhdastro/IMHDastro_main
	imhdscalar/IMHDscalar_main
	iscalar/Iscalar_main
	mrbc/MRBC_main
	mem_usage
	)

#TARGET_LINK_LIBRARIES(tarang fluid_incompressible fluid_base fields fff_pencil fff_slab  sff_pencil sff_slab ssf_pencil ssf_slab sss_pencil sss_slab chff_slab basis_basicfn global ${SYSTEM_LIBRARIES})

# To do fffw_slab cfff_slab
TARGET_LINK_LIBRARIES(tarang fluid_incompressible fluid_base fields sss_pencil ssf_pencil sff_pencil fff_pencil  basis_basicfn global ${SYSTEM_LIBRARIES})

ADD_CUSTOM_TARGET(
    tarang_source_dir_message
    COMMAND echo 
    COMMAND echo 'Compiling  from source directory: \"${PROJECT_SOURCE_DIR}\"'
    COMMAND echo 
)
ADD_DEPENDENCIES(tarang tarang_source_dir_message)


ADD_CUSTOM_COMMAND(TARGET tarang
                   POST_BUILD
                   COMMAND mv -f tarang ${PROJECT_BINARY_DIR}
                   COMMENT "moving ${EXECUTABLE_NAME} executable to PROJECT_BINARY_DIR" )



