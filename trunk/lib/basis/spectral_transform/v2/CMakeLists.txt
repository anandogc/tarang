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

ADD_LIBRARY(spectral_transform
spectral_transform.cc
spectral_plan.cc


plans/FFF_slab/fff_slab_isend_recv_overlap_isend_forward_3d.cc
plans/FFF_slab/fff_slab_isend_recv_overlap_isend_both_3d.cc
plans/FFF_slab/fff_slab_alltoall_3d.cc

plans/SFF_slab/sff_slab_isend_recv_3d.cc
plans/SFF_slab/sff_slab_alltoall_3d.cc

)
