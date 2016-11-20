/* Tarang-2
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-2 .
 *
 * Tarang-2 is free software; you can redistribute it and/or
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
 */

/*! \file  Global.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author  A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

//*********************************************************************************************

#include "Global.h"
#include <limits>

#define XSTR(x) #x
#define STR(x) XSTR(x)


Global::mpi::mpi():master_id(0){}

#ifdef FLOAT_Real
Global::io::io():output_precision(6){}
#else
Global::io::io():output_precision(12){}
#endif	

Global::program::program():
    version(STR(VERSION)){}

Global::myconstant::myconstant():
    I(Complex(0,1.0)),
    minusI(Complex(0,-1.0)),
    minus2I(Complex(0,-2.0)),
    MYEPS(1E-15),
    MYEPS2(1E-5),
	MY_MAX_INT(numeric_limits<int>::infinity()),
    INF_RADIUS(numeric_limits<Real>::infinity()),
    INF_TIME(numeric_limits<Real>::infinity()),
    MAX_NO_GLOB_BUFFER_PACKETS(1),
    MAX_NO_PROBE_PACKETS(5) {} 

//*********************************************************************************************

bool Global::Input_provided(const YAML::Node& node, const string parameter)
{
    if (!node.FindValue(parameter))
        return false;
        
    int type = node[parameter].Type();
    switch (type)
    {
        case YAML::NodeType::Scalar:
            return true;
            break;
        case YAML::NodeType::Sequence:
            return (node[parameter].size() > 0);
            break;
        case YAML::NodeType::Map:
            return true;
            break;
        case YAML::NodeType::Null:
            return false;
            break;
    }

    return false;
}


void Global::Show_error(string message)
{
	if (mpi.master){
		cerr << "ERROR: "+message << endl;
	}
	exit(1);
}

Global::Global()
{
	//FFF
	basis_table["FFF"]["FFF"]["Vx"]  ="FFF";
	basis_table["FFF"]["FFF"]["Vy"]  ="FFF";
	basis_table["FFF"]["FFF"]["Vz"]  ="FFF";
	basis_table["FFF"]["FFF"]["div"] ="FFF";
	basis_table["FFF"]["FFF"]["Vi2"] ="FFF";
	basis_table["FFF"]["FFF"]["VxVy"]="FFF";
	basis_table["FFF"]["FFF"]["VxVz"]="FFF";
	basis_table["FFF"]["FFF"]["VyVz"]="FFF";
	
	// ChFF: Chebyshev basis
	basis_table["ChFF"]["ChFF"]["Vx"]  ="ChFF";
	basis_table["ChFF"]["ChFF"]["Vy"]  ="ChFF";
	basis_table["ChFF"]["ChFF"]["Vz"]  ="ChFF";
	basis_table["ChFF"]["ChFF"]["div"] ="ChFF";
	basis_table["ChFF"]["ChFF"]["Vi2"] ="ChFF";
	basis_table["ChFF"]["ChFF"]["VxVy"]="ChFF";
	basis_table["ChFF"]["ChFF"]["VxVz"]="ChFF";
	basis_table["ChFF"]["ChFF"]["VyVz"]="ChFF";
	
	// SFF
	basis_table["SFF"]["SFF"]["Vx"]  ="SFF";
	basis_table["SFF"]["SFF"]["Vy"]  ="CFF";
	basis_table["SFF"]["SFF"]["Vz"]  ="CFF";
	basis_table["SFF"]["SFF"]["div"] ="CFF";
	basis_table["SFF"]["SFF"]["Vi2"] ="CFF";
	basis_table["SFF"]["SFF"]["VxVy"]="SFF";
	basis_table["SFF"]["SFF"]["VxVz"]="SFF";
	basis_table["SFF"]["SFF"]["VyVz"]="CFF";
	
	basis_table["SFF"]["CFF"]["Vx"]  ="CFF";
	basis_table["SFF"]["CFF"]["Vy"]  ="SFF";
	basis_table["SFF"]["CFF"]["Vz"]  ="SFF";
	basis_table["SFF"]["CFF"]["div"] ="SFF";
	basis_table["SFF"]["CFF"]["Vi2"] ="CFF";
	basis_table["SFF"]["CFF"]["VxVy"]="SFF";
	basis_table["SFF"]["CFF"]["VxVz"]="SFF";
	basis_table["SFF"]["CFF"]["VyVz"]="CFF";
	
	basis_table["SFF"]["S0F"]["Vx"]  ="S0F";
	basis_table["SFF"]["S0F"]["Vy"]  ="C0F";
	basis_table["SFF"]["S0F"]["Vz"]  ="C0F";
	basis_table["SFF"]["S0F"]["div"] ="C0F";
	basis_table["SFF"]["S0F"]["Vi2"] ="C0F";
	basis_table["SFF"]["S0F"]["VxVy"]="S0F";
	basis_table["SFF"]["S0F"]["VxVz"]="S0F";
	basis_table["SFF"]["S0F"]["VyVz"]="C0F";
	
	basis_table["SFF"]["C0F"]["Vx"]  ="C0F";
	basis_table["SFF"]["C0F"]["Vy"]  ="S0F";
	basis_table["SFF"]["C0F"]["Vz"]  ="S0F";
	basis_table["SFF"]["C0F"]["div"] ="S0F";
	basis_table["SFF"]["C0F"]["Vi2"] ="C0F";
	basis_table["SFF"]["C0F"]["VxVy"]="S0F";
	basis_table["SFF"]["C0F"]["VxVz"]="S0F";
	basis_table["SFF"]["C0F"]["VyVz"]="C0F";
	
	basis_table["SFF"]["SCF"]["Vx"]  ="SCF";
	basis_table["SFF"]["SCF"]["Vy"]  ="CSF";
	basis_table["SFF"]["SCF"]["Vz"]  ="CCF";
	basis_table["SFF"]["SCF"]["div"] ="CCF";
	basis_table["SFF"]["SCF"]["Vi2"] ="CCF";
	basis_table["SFF"]["SCF"]["VxVy"]="SSF";
	basis_table["SFF"]["SCF"]["VxVz"]="SCF";
	basis_table["SFF"]["SCF"]["VyVz"]="CSF";
	
	//SSF
	basis_table["SSF"]["SCF"]["Vx"]  ="SCF";
	basis_table["SSF"]["SCF"]["Vy"]  ="CSF";
	basis_table["SSF"]["SCF"]["Vz"]  ="CCF";
	basis_table["SSF"]["SCF"]["div"] ="CCF";
	basis_table["SSF"]["SCF"]["Vi2"] ="CCF";
	basis_table["SSF"]["SCF"]["VxVy"]="SSF";
	basis_table["SSF"]["SCF"]["VxVz"]="SCF";
	basis_table["SSF"]["SCF"]["VyVz"]="CSF";
	
	basis_table["SSF"]["CSF"]["Vx"]  ="CSF";
	basis_table["SSF"]["CSF"]["Vy"]  ="SCF";
	basis_table["SSF"]["CSF"]["Vz"]  ="SSF";
	basis_table["SSF"]["CSF"]["div"] ="SSF";
	basis_table["SSF"]["CSF"]["Vi2"] ="CCF";
	basis_table["SSF"]["CSF"]["VxVy"]="SSF";
	basis_table["SSF"]["CSF"]["VxVz"]="SCF";
	basis_table["SSF"]["CSF"]["VyVz"]="CSF";
	
	basis_table["SSF"]["SSF"]["Vx"]  ="SSF";
	basis_table["SSF"]["SSF"]["Vy"]  ="CCF";
	basis_table["SSF"]["SSF"]["Vz"]  ="CSF";
	basis_table["SSF"]["SSF"]["div"] ="CSF";
	basis_table["SSF"]["SSF"]["Vi2"] ="CCF";
	basis_table["SSF"]["SSF"]["VxVy"]="SSF";
	basis_table["SSF"]["SSF"]["VxVz"]="SCF";
	basis_table["SSF"]["SSF"]["VyVz"]="CSF";
	
	basis_table["SSF"]["CCF"]["Vx"]  ="CCF";
	basis_table["SSF"]["CCF"]["Vy"]  ="SSF";
	basis_table["SSF"]["CCF"]["Vz"]  ="SCF";
	basis_table["SSF"]["CCF"]["div"] ="SCF";
	basis_table["SSF"]["CCF"]["Vi2"] ="CCF";
	basis_table["SSF"]["CCF"]["VxVy"]="SSF";
	basis_table["SSF"]["CCF"]["VxVz"]="SCF";
	basis_table["SSF"]["CCF"]["VyVz"]="CSF";
	
	//SSS
	basis_table["SSS"]["SCC"]["Vx"]  ="SCC";
	basis_table["SSS"]["SCC"]["Vy"]  ="CSC";
	basis_table["SSS"]["SCC"]["Vz"]  ="CCS";
	basis_table["SSS"]["SCC"]["div"] ="CCC";
	basis_table["SSS"]["SCC"]["Vi2"] ="CCC";
	basis_table["SSS"]["SCC"]["VxVy"]="SSC";
	basis_table["SSS"]["SCC"]["VxVz"]="SCS";
	basis_table["SSS"]["SCC"]["VyVz"]="CSS";
	
	basis_table["SSS"]["CSS"]["Vx"]  ="CSS";
	basis_table["SSS"]["CSS"]["Vy"]  ="SCS";
	basis_table["SSS"]["CSS"]["Vz"]  ="SSC";
	basis_table["SSS"]["CSS"]["div"] ="SSS";
	basis_table["SSS"]["CSS"]["Vi2"] ="CCC";
	basis_table["SSS"]["CSS"]["VxVy"]="SSC";
	basis_table["SSS"]["CSS"]["VxVz"]="SCS";
	basis_table["SSS"]["CSS"]["VyVz"]="CSS";
	
	basis_table["SSS"]["CCS"]["Vx"]  ="CCS";
	basis_table["SSS"]["CCS"]["Vy"]  ="SSS";
	basis_table["SSS"]["CCS"]["Vz"]  ="SCC";
	basis_table["SSS"]["CCS"]["div"] ="SCS";
	basis_table["SSS"]["CCS"]["Vi2"] ="CCC";
	basis_table["SSS"]["CCS"]["VxVy"]="SSC";
	basis_table["SSS"]["CCS"]["VxVz"]="SCS";
	basis_table["SSS"]["CCS"]["VyVz"]="CSS";
	
	basis_table["SSS"]["SSC"]["Vx"]  ="SSC";
	basis_table["SSS"]["SSC"]["Vy"]  ="CCC";
	basis_table["SSS"]["SSC"]["Vz"]  ="CSS";
	basis_table["SSS"]["SSC"]["div"] ="CSC";
	basis_table["SSS"]["SSC"]["Vi2"] ="CCC";
	basis_table["SSS"]["SSC"]["VxVy"]="SSC";
	basis_table["SSS"]["SSC"]["VxVz"]="SCS";
	basis_table["SSS"]["SSC"]["VyVz"]="CSS";
	
	basis_table["SSS"]["CSC"]["Vx"]  ="CSC";
	basis_table["SSS"]["CSC"]["Vy"]  ="SCC";
	basis_table["SSS"]["CSC"]["Vz"]  ="SSS";
	basis_table["SSS"]["CSC"]["div"] ="SSC";
	basis_table["SSS"]["CSC"]["Vi2"] ="CCC";
	basis_table["SSS"]["CSC"]["VxVy"]="SSC";
	basis_table["SSS"]["CSC"]["VxVz"]="SCS";
	basis_table["SSS"]["CSC"]["VyVz"]="CSS";
	
	basis_table["SSS"]["SCS"]["Vx"]  ="SCS";
	basis_table["SSS"]["SCS"]["Vy"]  ="CSS";
	basis_table["SSS"]["SCS"]["Vz"]  ="CCC";
	basis_table["SSS"]["SCS"]["div"] ="CCS";
	basis_table["SSS"]["SCS"]["Vi2"] ="CCC";
	basis_table["SSS"]["SCS"]["VxVy"]="SSC";
	basis_table["SSS"]["SCS"]["VxVz"]="SCS";
	basis_table["SSS"]["SCS"]["VyVz"]="CSS";
	
	basis_table["SSS"]["CCC"]["Vx"]  ="CCC";
	basis_table["SSS"]["CCC"]["Vy"]  ="SSC";
	basis_table["SSS"]["CCC"]["Vz"]  ="SCS";
	basis_table["SSS"]["CCC"]["div"] ="SCC";
	basis_table["SSS"]["CCC"]["Vi2"] ="CCC";
	basis_table["SSS"]["CCC"]["VxVy"]="SSC";
	basis_table["SSS"]["CCC"]["VxVz"]="SCS";
	basis_table["SSS"]["CCC"]["VyVz"]="CSS";
	
	basis_table["SSS"]["SSS"]["Vx"]  ="SSS";
	basis_table["SSS"]["SSS"]["Vy"]  ="CCS";
	basis_table["SSS"]["SSS"]["Vz"]  ="CSC";
	basis_table["SSS"]["SSS"]["div"] ="CSS";
	basis_table["SSS"]["SSS"]["Vi2"] ="CCC";
	basis_table["SSS"]["SSS"]["VxVy"]="SSC";
	basis_table["SSS"]["SSS"]["VxVz"]="SCS";
	basis_table["SSS"]["SSS"]["VyVz"]="CSS";
	
	basis_table["SSS"]["S0S"]["Vx"]  ="S0S";
	basis_table["SSS"]["S0S"]["Vy"];
	basis_table["SSS"]["S0S"]["Vz"]  ="C0C";
	basis_table["SSS"]["S0S"]["div"] ="C0S";
	basis_table["SSS"]["S0S"]["Vi2"] ="C0C";
	basis_table["SSS"]["S0S"]["VxVy"];
	basis_table["SSS"]["S0S"]["VxVz"]="S0S";
	basis_table["SSS"]["S0S"]["VyVz"];
	
	basis_table["SSS"]["C0C"]["Vx"]  ="C0C";
	basis_table["SSS"]["C0C"]["Vy"];
	basis_table["SSS"]["C0C"]["Vz"]  ="S0S";
	basis_table["SSS"]["C0C"]["div"] ="S0C";
	basis_table["SSS"]["C0C"]["Vi2"] ="C0C";
	basis_table["SSS"]["C0C"]["VxVy"];
	basis_table["SSS"]["C0C"]["VxVz"]="S0S";
	basis_table["SSS"]["C0C"]["VyVz"];
	
	basis_table["SSS"]["S0C"]["Vx"]  ="S0C";
	basis_table["SSS"]["S0C"]["Vy"]  ="C0C";
	basis_table["SSS"]["S0C"]["Vz"]  ="C0S";
	basis_table["SSS"]["S0C"]["div"] ="C0C";
	basis_table["SSS"]["S0C"]["Vi2"] ="C0C";
	basis_table["SSS"]["S0C"]["VxVy"];
	basis_table["SSS"]["S0C"]["VxVz"]="S0S";
	basis_table["SSS"]["S0C"]["VyVz"];
	
	basis_table["SSS"]["C0S"]["Vx"]  ="C0S";
	basis_table["SSS"]["C0S"]["Vy"];
	basis_table["SSS"]["C0S"]["Vz"]  ="S0C";
	basis_table["SSS"]["C0S"]["div"] ="S0S";
	basis_table["SSS"]["C0S"]["Vi2"] ="C0C";
	basis_table["SSS"]["C0S"]["VxVy"];
	basis_table["SSS"]["C0S"]["VxVz"]="S0S";
	basis_table["SSS"]["C0S"]["VyVz"];
	
	
	// no of components for reading modes in initial condition and modes in force
	//no_table[SOLVER_NAME] 
    no_components_table["FLUID_INCOMPRESS"]=2;
    no_components_table["FLUID_COMPRESS"]=3;
    
    no_components_table["RBC"]=3;
    no_components_table["STRATIFIED"]=3;
	no_components_table["MRBC"]=4;
    
    no_components_table["SCALAR_INCOMPRESS"]=3;
    no_components_table["SCALAR_COMPRESS"]=4;
    
    no_components_table["MHD_INCOMPRESS"]=4;
    no_components_table["MHD_COMPRESS"]=6;
    
    no_components_table["KEPLERIAN"]=4;
    no_components_table["KEPLERIAN_COMPRESS"]=6;
    
    no_components_table["MHD_SCALAR_INCOMPRESS"]=5;
    no_components_table["RBMC"];
	
	no_components_table["MHD_ASTRO_INCOMPRESS"]=6;
    
    //global data packet size
    //global_data_packet_size_table[SOLVER_NAME]
    global_data_packet_size_table["FLUID_INCOMPRESS"]=15;
    global_data_packet_size_table["RBC"]=21;
    global_data_packet_size_table["STRATIFIED"]=21;
	global_data_packet_size_table["MRBC"]=24;
    global_data_packet_size_table["SCALAR_INCOMPRESS"]=20;
    global_data_packet_size_table["MHD_INCOMPRESS"]=29;
    global_data_packet_size_table["KEPLERIAN"]=29;
    global_data_packet_size_table["MHD_SCALAR_INCOMPRESS"]=34;
	global_data_packet_size_table["MHD_ASTRO_INCOMPRESS"]=38;
    global_data_packet_size_table["RBMC"];
    
    //spectral probe packet size
    //spectral_probe_packet_size_table[SOLVER_NAME][BASIS_TYPE]
	// (global.io.probes.spectral_space.number * 6) +1 + Tk for each field; 
    spectral_probe_packet_size_table["FLUID_INCOMPRESS"]["SSS"]=8;
    spectral_probe_packet_size_table["FLUID_INCOMPRESS"]["SSF"]=11;
    spectral_probe_packet_size_table["FLUID_INCOMPRESS"]["SFF"]=11;
    spectral_probe_packet_size_table["FLUID_INCOMPRESS"]["FFF"]=11;
    spectral_probe_packet_size_table["FLUID_INCOMPRESS"]["FFFW"]=11;
    
    spectral_probe_packet_size_table["RBC"]["SSS"]=10;
    spectral_probe_packet_size_table["RBC"]["SSF"]=14;
    spectral_probe_packet_size_table["RBC"]["SFF"]=14;
    spectral_probe_packet_size_table["RBC"]["FFF"]=14;
    spectral_probe_packet_size_table["RBC"]["FFFW"]=14;
    
    spectral_probe_packet_size_table["STRATIFIED"]["SSS"]=10;
    spectral_probe_packet_size_table["STRATIFIED"]["SSF"]=14;
    spectral_probe_packet_size_table["STRATIFIED"]["SFF"]=14;
    spectral_probe_packet_size_table["STRATIFIED"]["FFF"]=14;
    spectral_probe_packet_size_table["STRATIFIED"]["FFFW"]=14;
    
    spectral_probe_packet_size_table["SCALAR_INCOMPRESS"]["SSS"]=10;
    spectral_probe_packet_size_table["SCALAR_INCOMPRESS"]["SSF"]=14;
    spectral_probe_packet_size_table["SCALAR_INCOMPRESS"]["SFF"]=14;
    spectral_probe_packet_size_table["SCALAR_INCOMPRESS"]["FFF"]=14;
    spectral_probe_packet_size_table["SCALAR_INCOMPRESS"]["FFFW"]=14;
	
	spectral_probe_packet_size_table["MRBC"]["SSS"]=12;
    spectral_probe_packet_size_table["MRBC"]["SSF"]=17;
    spectral_probe_packet_size_table["MRBC"]["SFF"]=17;
    spectral_probe_packet_size_table["MRBC"]["FFF"]=17;
    spectral_probe_packet_size_table["MRBC"]["FFFW"]=17;

    spectral_probe_packet_size_table["MHD_INCOMPRESS"]["SSS"]=12;
    spectral_probe_packet_size_table["MHD_INCOMPRESS"]["SSF"]=18;
    spectral_probe_packet_size_table["MHD_INCOMPRESS"]["SFF"]=18;
    spectral_probe_packet_size_table["MHD_INCOMPRESS"]["FFF"]=18;
    spectral_probe_packet_size_table["MHD_INCOMPRESS"]["FFFW"]=18;
    
    spectral_probe_packet_size_table["KEPLERIAN"]["FFF"]=18;
    
    spectral_probe_packet_size_table["MHD_SCALAR_INCOMPRESS"];
    spectral_probe_packet_size_table["RBMC"];

	//real probe packet size
    //real_probe_packet_size_table[SOLVER_NAME]
    real_probe_packet_size_table["FLUID_INCOMPRESS"]=7;
    real_probe_packet_size_table["RBC"]=8;
    real_probe_packet_size_table["STRATIFIED"]=8;
    real_probe_packet_size_table["SCALAR_INCOMPRESS"]=8;
	real_probe_packet_size_table["MRBC"]=9;
    real_probe_packet_size_table["MHD_INCOMPRESS"]=10;
    real_probe_packet_size_table["KEPLERIAN"]=10;
}

void Global::io::Dump_buffer(ofstream& file, Array<Real,1> buffer, unsigned long& buffer_index, unsigned long packet_size){
	unsigned long num_packets=buffer_index/packet_size;
	unsigned long i=0;
	
	for (i=0; i<num_packets; i++){
		copy(buffer.data()+i*packet_size, buffer.data()+(i+1)*packet_size, ostream_iterator<Real>(file, "\t"));
		file << '\n';
	}
	
	file.flush();
	
	buffer_index = 0;
}

int Global::io::init_cond_modes::In_table(int kx, int ky, int kz)
{
    for (int i=0; i<number; i++)
        if ( (kx == coords(i,1)) && (ky == coords(i,2)) && (kz == coords(i,3)) ) {
            return i;
        }

    return -1;
}
