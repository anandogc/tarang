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
#CXX=mpicxx cmake [-DREAL=FLOAT] [-DGROSSMANN_LOHSE=ON] [-DFIND_LIBRARIES=OFF] [-DCMAKE_BUILD_TYPE=DEBUG] [-DFFTW_PLAN=MEASURE] /path/to/CMakeLists.txt

#Example:
#CC=mpicc CXX=mpicxx cmake ../tarang/trunk
############################################

#Default options
SET (DEFAULT_REAL "DOUBLE")                      #FLOAT, DOUBLE
SET (DEFAULT_GROSSMANN_LOHSE OFF)                #ON, OFF
SET (DEFAULT_FIND_LIBRARIES ON)                  #ON, OFF - Try to locate packages using pkg-config when ON
SET (DEFAULT_CMAKE_BUILD_TYPE "DISTRIBUTION")    #DISTRIBUTION, DEBUG
SET (DEFAULT_FFTW_PLAN "MEASURE")                #ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE

SET (INCLUDE_DIRS
        #/path/to/include
    )
SET (LIBRARY_DIRS
       #/path/to/lib
    )


#Required libraries, One of FFTW_FLOAT_LIBS or FFTW_DOUBLE_LIBS is selected depending on REAL
SET (FFTW_FLOAT_LIBS fftkf fftw3f_mpi fftw3f)
SET (FFTW_DOUBLE_LIBS fftk fftw3_mpi fftw3)
SET (OTHER_LIBS blitz m hdf5 h5si dl yaml-cpp)


#################################################
#Users need not modify any thing beyond this line
#################################################

#Set version
SET (VERSION "2.5.2-dev")

SET (SPECTRAL_TRANSFORM_VERSION v2)

#Set default value for REAL
IF (NOT DEFINED REAL)
    SET(REAL ${DEFAULT_REAL})
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

SET (REAL ${REAL} CACHE STRING "REAL")
SET (FIND_LIBRARIES ${FIND_LIBRARIES} CACHE BOOL "FIND_LIBRARIES")
SET (FFTW_PLAN ${FFTW_PLAN} CACHE STRING "FFTW_PLAN")

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
IF (REAL MATCHES "FLOAT")
    ADD_DEFINITIONS (-DREAL_FLOAT)
ELSEIF(REAL MATCHES "DOUBLE")
    ADD_DEFINITIONS (-DREAL_DOUBLE)
ELSE()
    MESSAGE(FATAL_ERROR "REAL must be FLOAT or DOUBLE")
ENDIF()

IF (NOT REAL_DISPLAYED)
	MESSAGE(STATUS "REAL: ${REAL}")
    SET(REAL_DISPLAYED TRUE CACHE BOOL "REAL_DISPLAYED")
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
        
        IF (REAL MATCHES "FLOAT")
            PKG_CHECK_MODULES(FFTW3F QUIET fftw3f)
            DISPLAY_FOUND_LIBRARY(FFTW3F FFTW3F "${FFTW3F_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig") 
        ELSEIF(REAL MATCHES "DOUBLE")
            PKG_CHECK_MODULES(FFTW3 QUIET fftw3)
            DISPLAY_FOUND_LIBRARY(FFTW3 FFTW3 "${FFTW3_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig")
        ELSE()
            MESSAGE(FATAL_ERROR "REAL must be FLOAT or DOUBLE")
        ENDIF()
        
        FIND_PACKAGE(HDF5 QUIET)
        IF (HDF5_FOUND)  
            #Find installed version of HDF5
            EXECUTE_PROCESS(COMMAND ${HDF5_DIFF_EXECUTABLE} --version
                            COMMAND cut -d\  -f3
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            OUTPUT_VARIABLE HDF5_VERSION
                            ERROR_VARIABLE HDF5_VERSION_ERROR)
            DISPLAY_FOUND_LIBRARY(HDF5 HDF5 "${HDF5_INCLUDE_DIRS}" "Set HDF5_ROOT, e.g. /usr/local")
            IF (HDF5_IS_PARALLEL) #Check whether HDF5 has parallen installation or not
                    MESSAGE(STATUS "HDF5 supports parallel IO: yes")
            ELSE()
                    MESSAGE(FATAL_ERROR "HDF5 supports parallel IO: no")
            ENDIF()
        ENDIF()
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
${CMAKE_SOURCE_DIR}/lib/basis/spectral_transform/
${CMAKE_SOURCE_DIR}/lib/basis/basicfn
${CMAKE_SOURCE_DIR}/lib/basis/basicfn/h5si
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
${CACHED_INCLUDE_DIRS}
)

#Folders where libraries will be searched for
LINK_DIRECTORIES(
${YAML_CPP_LIBRARY_DIRS}
${BLITZ_LIBRARY_DIRS}
${FFTW3F_LIBRARY_DIRS}
${FFTW3_LIBRARY_DIRS}
${HDF5_LIBRARY_DIRS}
${LIBRARY_DIRS}
${CACHED_LIBRARY_DIRS}
)

ADD_DEFINITIONS(-DVERSION=${VERSION} ${HDF5_DEFINITIONS})

IF (REAL MATCHES "FLOAT")
    SET(SYSTEM_LIBRARIES ${FFTW_FLOAT_LIBS} ${FFTW3F_LIBRARIES})
ELSEIF(REAL MATCHES "DOUBLE")
    SET(SYSTEM_LIBRARIES ${FFTW_DOUBLE_LIBS} ${FFTW3_LIBRARIES})
ENDIF()

SET(SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OTHER_LIBS} ${HDF5_LIBRARIES} ${FFTW3_LIBRARIES} ${BLITZ_LIBRARIES} ${YAML_CPP_LIBRARIES} ${CACHED_LIBRARIES} )

# MESSAGE(STATUS "SYSTEM_LIBRARIES = ${SYSTEM_LIBRARIES}")
# SET (SYSTEM_LIBRARIES CACHE BOOL "CMAKE_BUILD_TYPE_DISPLAYED")

#  IF (NOT CACHED_INCLUDE_DIRS)
#      SET (CACHED_INCLUDE_DIRS ${YAML_CPP_INCLUDE_DIRS} ${BLITZ_INCLUDE_DIRS} ${FFTW3F_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${INCLUDE_DIRS} CACHE STRING "CACHED_INCLUDE_DIRS")
#  ENDIF ()

#  IF (NOT CACHED_LIBRARY_DIRS)
#      SET (CACHED_LIBRARY_DIRS ${YAML_CPP_LIBRARY_DIRS} ${BLITZ_LIBRARY_DIRS} ${FFTW3F_LIBRARY_DIRS} ${FFTW3_LIBRARY_DIRS} ${HDF5_LIBRARY_DIRS} ${LIBRARY_DIRS} CACHE STRING "CACHED_LIBRARY_DIRS")
#  ENDIF ()

#  IF (NOT CACHED_LIBRARIES)
#      SET (CACHED_LIBRARIES ${SYSTEM_LIBRARIES} CACHE STRING "CACHED_LIBRARIES")
#  ENDIF ()


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

SET (CONFIGURE_OPTIONS '"CXX=${CMAKE_CXX_COMPILER} cmake -DREAL=${REAL} -DGROSSMANN_LOHSE=${GROSSMANN_LOHSE} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DFIND_LIBRARIES=${FIND_LIBRARIES} -DFFTW_PLAN=${FFTW_PLAN} ${TARANG_SOURCE_DIR}"')
SET(CMAKE_CXX_FLAGS "-DCONFIGURE_OPTIONS=${CONFIGURE_OPTIONS}")


SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -Wno-deprecated -DH5SI_ENABLE_MPI")

SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DH5SI_ENABLE_MPI")

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(src)



