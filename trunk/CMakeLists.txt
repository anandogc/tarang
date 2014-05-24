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
 # @author  M. K. Verma, A. G. Chatterjee


PROJECT(TARANG CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)

######################
#Environment Variables
######################
#export PATH=$HOME/local/bin:$PATH
#export PKG_CONFIG_DISABLE_UNINSTALLED=true
#export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH
#export HDF5_ROOT=$HOME/local
#export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
######################

#CMake usage:
#CXX=mpicxx cmake [-DDP=FLOAT] [-DGROSSMANN_LOHSE=ON] [-DFIND_LIBRARIES=OFF] [-DPOSTFIX=".mic"] [-DCMAKE_BUILD_TYPE=DEBUG] [-DMMIC] [-DFFTW_PLAN=MEASURE] /path/to/CMakeLists.txt

#Example:
#CC=mpicc CXX=mpicxx cmake ../tarang/trunk
############################################

#Default options
SET (DEFAULT_DP "DOUBLE")                        #FLOAT, DOUBLE
SET (DEFAULT_GROSSMANN_LOHSE OFF)                #ON, OFF
SET (DEFAULT_FIND_LIBRARIES ON)                  #ON, OFF - Try to locate packages using pkg-config when ON
SET (DEFAULT_CMAKE_BUILD_TYPE "DISTRIBUTION")    #DISTRIBUTION, DEBUG
SET (DEFAULT_FFTW_PLAN "MEASURE")                #ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE

SET (INCLUDE_DIRS
       #/path/to/include2
    )
SET (LIBRARY_DIRS
       #/path/to/lib2
    )

SET (SPECTRAL_TRANSFORM_VERSION v2)

#Required libraries, One of FFTW_FLOAT_LIBS or FFTW_DOUBLE_LIBS is selected depending on DP
SET (FFTW_FLOAT_LIBS fftw3f_mpi fftw3f)
SET (FFTW_DOUBLE_LIBS fftw3_mpi fftw3)
SET (OTHER_LIBS blitz m hdf5 yaml-cpp dl z)


#################################################
#Users need not modify any thing beyond this line
#################################################

#Set version
SET (VERSION "2.4.0")

#Set default value for DP
IF (NOT DEFINED DP)
    SET(DP ${DEFAULT_DP})
ENDIF()

#Set default value for GROSSMANN_LOHSE
IF (NOT DEFINED GROSSMANN_LOHSE)
    SET(GROSSMANN_LOHSE ${DEFAULT_GROSSMANN_LOHSE})
ENDIF()

#Set default value for FIND_LIBRARIES
IF(NOT DEFINED FIND_LIBRARIES)
    SET (FIND_LIBRARIES ${DEFAULT_FIND_LIBRARIES})
ENDIF()

#Set default value for CMAKE_BUILD_TYPE
IF(NOT CMAKE_BUILD_TYPE)
   SET (CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE})
ENDIF()

#Set default value for FFTW_PLAN
IF (NOT DEFINED FFTW_PLAN)
    SET(FFTW_PLAN ${DEFAULT_FFTW_PLAN})
ENDIF()

#Check VERSION is not empty
IF (NOT VERSION)
    MESSAGE (FATAL_ERROR "VERSION can not be empty")
ENDIF()

SET (DP ${DP} CACHE STRING "DP")
SET (FIND_LIBRARIES ${FIND_LIBRARIES} CACHE BOOL "FIND_LIBRARIES")
SET (FFTW_PLAN ${FFTW_PLAN} CACHE STRING "FFTW_PLAN")


#Check CC compiler supports MPI
#IF (NOT CC_SUPPORTS_MPI)        #CC_SUPPORTS_MPI is undefined when run for the first time.
#    MESSAGE("")
#    MESSAGE(STATUS "MPI support in C compiler: ${CMAKE_C_COMPILER}")
#    FILE(WRITE ${PROJECT_BINARY_DIR}/CMakeFiles/C_has_MPI.c 
#               "#include <mpi.h>
#                int main(int argc, char** argv)
#                { MPI_Init(&argc, &argv);
#                  MPI_Finalize();} ")

#    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} ${PROJECT_BINARY_DIR}/CMakeFiles/C_has_MPI.c
#                    ERROR_VARIABLE err
#                    RESULT_VARIABLE CC_HAS_MPI)

#    IF (CC_HAS_MPI MATCHES 0)
#        MESSAGE(STATUS "MPI support in C compiler: ${CMAKE_C_COMPILER} -- supported")
#    ELSE()
#        MESSAGE(FATAL_ERROR "${CMAKE_C_COMPILER} does not support MPI")
#    ENDIF()

#    SET (CC_SUPPORTS_MPI TRUE CACHE BOOL "CC_SUPPORTS_MPI")
#ENDIF()

#Check CXX compiler supports MPI
IF (NOT CXX_SUPPORTS_MPI)        #CXX_SUPPORTS_MPI is undefined when run for the first time.
    MESSAGE(STATUS "MPI support in CXX compiler: ${CMAKE_CXX_COMPILER}")
    FILE(WRITE ${PROJECT_BINARY_DIR}/CMakeFiles/CXX_has_MPI.cc 
               "#include <mpi.h>
                int main(int argc, char** argv)
                { MPI_Init(&argc, &argv);
                  MPI_Finalize();} ")

    EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} ${PROJECT_BINARY_DIR}/CMakeFiles/CXX_has_MPI.cc
                    ERROR_VARIABLE err
                    RESULT_VARIABLE CXX_HAS_MPI)

    IF (CXX_HAS_MPI MATCHES 0)
        MESSAGE(STATUS "MPI support in CXX compiler: ${CMAKE_CXX_COMPILER} -- supported")
    ELSE()
        MESSAGE(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not support MPI")
    ENDIF()
    SET (CXX_SUPPORTS_MPI TRUE CACHE BOOL "CXX_SUPPORTS_MPI")
ENDIF()

#Configure for FFTW library
MESSAGE("")
IF (DP MATCHES "FLOAT")
    ADD_DEFINITIONS (-DFLOAT_DP)
ELSEIF(DP MATCHES "DOUBLE")
    ADD_DEFINITIONS (-DDOUBLE_DP)
ELSE()
    MESSAGE(FATAL_ERROR "DP must be FLOAT or DOUBLE")
ENDIF()

IF (NOT DP_DISPLAYED)
	MESSAGE(STATUS "DP: ${DP}")
    SET(DP_DISPLAYED TRUE CACHE BOOL "DP_DISPLAYED")
ENDIF()

#Find the required libraries
IF (NOT ALL_LIBRARIES_FOUND AND FIND_LIBRARIES)
    FUNCTION(DISPLAY_FOUND_LIBRARY DISPLAY_NAME LIB_NAME INCLUDE_PATH ERROR_MESSAGE)
        IF (${LIB_NAME}_FOUND)
            EXECUTE_PROCESS(COMMAND echo ${INCLUDE_PATH}
                            COMMAND rev
                            COMMAND cut -d/ -f2-
                            COMMAND rev
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            OUTPUT_VARIABLE ${LIB_NAME}_ROOT_DIR)

            MESSAGE(STATUS "Found ${DISPLAY_NAME}: ${${LIB_NAME}_ROOT_DIR} (version \"${${LIB_NAME}_VERSION}\")")
        ELSE()
            MESSAGE(FATAL_ERROR "${DISPLAY_NAME} not found. ${ERROR_MESSAGE}")
        ENDIF()
    ENDFUNCTION()

#Check whether pkg-config is installed
    EXECUTE_PROCESS(COMMAND which pkg-config
                    RESULT_VARIABLE PKG_CONFIG_FOUND
                    OUTPUT_VARIABLE PKG_CONFIG_LOCATION)

    IF (PKG_CONFIG_FOUND MATCHES 0)
	    FIND_PACKAGE(PkgConfig REQUIRED)

        PKG_CHECK_MODULES(YAML_CPP QUIET yaml-cpp)
        DISPLAY_FOUND_LIBRARY(YAML-CPP YAML_CPP "${YAML_CPP_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig") 

        PKG_CHECK_MODULES(BLITZ QUIET blitz)
        DISPLAY_FOUND_LIBRARY(Blitz++ BLITZ "${BLITZ_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig") 

        IF (DP MATCHES "FLOAT")
            PKG_CHECK_MODULES(FFTW3F QUIET fftw3f)
            DISPLAY_FOUND_LIBRARY(FFTW3F FFTW3F "${FFTW3F_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig") 
        ELSEIF(DP MATCHES "DOUBLE")
            PKG_CHECK_MODULES(FFTW3 QUIET fftw3)
            DISPLAY_FOUND_LIBRARY(FFTW3 FFTW3 "${FFTW3_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig")
        ELSE()
            MESSAGE(FATAL_ERROR "DP must be FLOAT or DOUBLE")
        ENDIF()

        FIND_PACKAGE(HDF5 QUIET)
        IF (HDF5_FOUND)  
            #Find installed version of HDF5
            EXECUTE_PROCESS(COMMAND ${HDF5_DIFF_EXECUTABLE} --version
                            COMMAND cut -d\  -f3
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            OUTPUT_VARIABLE HDF5_VERSION
                            ERROR_VARIABLE HDF5_VERSION_ERROR)
        ENDIF()
        DISPLAY_FOUND_LIBRARY(HDF5 HDF5 "${HDF5_INCLUDE_DIRS}" "Set HDF5_ROOT, e.g. /usr/local")
    #IF (HDF5_FOUND) #Check whether HDF5 has parallen installation or not
	#		IF (HDF5_IS_PARALLEL)
	#			MESSAGE(STATUS "HDF5 supports parallel IO: yes")
	#		ELSE()
	#			MESSAGE(FATAL_ERROR "HDF5 supports parallel IO: no")
	#		ENDIF()
	#ENDIF()
	ELSE()
		MESSAGE(STATUS "Warning: pkg-config not found, libraries will not be searched for. Supply -DFIND_LIBRARIES=OFF to CMake to disable library search.")
	ENDIF()
       SET(ALL_LIBRARIES_FOUND TRUE CACHE BOOL "ALL_LIBRARIES_FOUND")
ENDIF()


#Set FFTW_PLAN
IF (FFTW_PLAN MATCHES "ESTIMATE")
    ADD_DEFINITIONS(-DESTIMATE)
ELSEIF (FFTW_PLAN MATCHES "MEASURE")
    ADD_DEFINITIONS(-DMEASURE)
ELSEIF (FFTW_PLAN MATCHES "PATIENT")
    ADD_DEFINITIONS(-DPATIENT)
ELSEIF (FFTW_PLAN MATCHES "EXHAUSTIVE")
    ADD_DEFINITIONS(-DPATIENT)
ELSE()
    MESSAGE(FATAL_ERROR "FFTW_PLAN can be one of ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE.")
ENDIF()


IF(GROSSMANN_LOHSE)
    ADD_DEFINITIONS(-DGROSSMANN_LOHSE)
    MESSAGE(STATUS "GROSSMANN_LOHSE enabled")
ENDIF()

IF (NOT FFTW_PLAN_DISPLAYED)
    MESSAGE(STATUS "FFTW_PLAN: ${FFTW_PLAN}")
    SET (FFTW_PLAN_DISPLAYED TRUE CACHE BOOL "FFTW_PLAN_DISPLAYED")
    MESSAGE("")
ENDIF()



#Folders where header files will be searched for
INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/global
${CMAKE_SOURCE_DIR}/lib/basis
${CMAKE_SOURCE_DIR}/lib/basis/spectral_transform/${SPECTRAL_TRANSFORM_VERSION}
${CMAKE_SOURCE_DIR}/lib/basis/basicfn
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/fluid
${CMAKE_SOURCE_DIR}/lib/fluid/fluid_base
${CMAKE_SOURCE_DIR}/lib/fluid/fluid_base/force
${CMAKE_SOURCE_DIR}/lib/fluid/fluid_base/io
${CMAKE_SOURCE_DIR}/lib/fluid/incompressible
${CMAKE_SOURCE_DIR}/lib/fluid/incompressible/energytr
${CMAKE_SOURCE_DIR}/lib/fluid/incompressible/incio
${CMAKE_SOURCE_DIR}/lib/fluid/incompressible/nlin
${CMAKE_SOURCE_DIR}/lib/fluid/incompressible/time_advance
${YAML_CPP_INCLUDE_DIRS}
${BLITZ_INCLUDE_DIRS}
${FFTW3F_INCLUDE_DIRS}
${FFTW3_INCLUDE_DIRS}
${HDF5_INCLUDE_DIRS}
${INCLUDE_DIRS}
)

#Folders where libraries will be searched for
LINK_DIRECTORIES(

# ${CMAKE_SOURCE_DIR}/lib/global
# ${CMAKE_SOURCE_DIR}/lib/basis/basicfn
# ${CMAKE_SOURCE_DIR}/lib/basis/fff/pencil
# ${CMAKE_SOURCE_DIR}/lib/basis/fff/slab
# ${CMAKE_SOURCE_DIR}/lib/basis/SFF/pencil
# ${CMAKE_SOURCE_DIR}/lib/basis/SFF/slab
# ${CMAKE_SOURCE_DIR}/lib/basis/SSF/pencil
# ${CMAKE_SOURCE_DIR}/lib/basis/SSF/slab
# ${CMAKE_SOURCE_DIR}/lib/basis/SSS/pencil
# ${CMAKE_SOURCE_DIR}/lib/basis/SSS/slab
# ${CMAKE_SOURCE_DIR}/lib/fields
# ${CMAKE_SOURCE_DIR}/lib/fluid/fluid_base
# ${CMAKE_SOURCE_DIR}/lib/fluid/Incompressible

${YAML_CPP_LIBRARY_DIRS}
${BLITZ_LIBRARY_DIRS}
${FFTW3F_LIBRARY_DIRS}
${FFTW3_LIBRARY_DIRS}
${HDF5_LIBRARY_DIRS}
${LIBRARY_DIRS}
)

ADD_DEFINITIONS(-DVERSION=${VERSION})

IF (DP MATCHES "FLOAT")
    SET(SYSTEM_LIBRARIES ${FFTW_FLOAT_LIBS} ${OTHER_LIBS})
ELSEIF(DP MATCHES "DOUBLE")
    SET(SYSTEM_LIBRARIES ${FFTW_DOUBLE_LIBS} ${OTHER_LIBS})
ENDIF()

#Build type
IF (NOT CMAKE_BUILD_TYPE_DISPLAYED)
    IF ( (NOT ${CMAKE_BUILD_TYPE} MATCHES "DISTRIBUTION") AND (NOT ${CMAKE_BUILD_TYPE} MATCHES "DEBUG") )
        MESSAGE(FATAL_ERROR "Build type can be one of DISTRIBUTION, DEBUG")
    ENDIF()

    MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    SET (CMAKE_BUILD_TYPE_DISPLAYED TRUE CACHE BOOL "CMAKE_BUILD_TYPE_DISPLAYED")
ENDIF ()

MESSAGE(STATUS "Source Folder: ${TARANG_SOURCE_DIR}")
MESSAGE("")

SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")

SET(CMAKE_CXX_FLAGS_DEBUG "-g")
#SET(CMAKE_C_FLAGS_DEBUG "-g")

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(src)



