# -*- coding: utf-8 -*-
# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

"""
Created on Thu Sep 29 17:12:39 2016

@author: py12rw
"""

import sys

try:
    from ffeatools.modules.FFEA_trajectory import FFEA_trajectory # python package
except:
    try:
        from ffeatools.ffea_trajectory import ffea_trajectory as FFEA_trajectory
    except ImportError:
        print("Failure to import FFEA_trajectory")
        sys.exit(1) # failure to import

try:
    test_load_traj = FFEA_trajectory("unit_test_traj.ftj")
    sys.exit(0)
except IOError:
    print("Couldn't find trajectory file that's supposed to be packed in with these tests. Probably a CMake issue!")
    sys.exit(1)
except Exception as e:
    print(e)
    sys.exit(1)
