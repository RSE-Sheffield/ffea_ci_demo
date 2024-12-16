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

set (vdwsiteposdir "${PROJECT_BINARY_DIR}/tests/rods/unit/rod_vdw_site_placement/")
file (COPY rod_vdw_site_placement.ffeatest DESTINATION ${vdwsiteposdir})
file (COPY params.yml DESTINATION ${vdwsiteposdir})
file (COPY create_rod_with_vdw_sites.py DESTINATION ${vdwsiteposdir})
file (COPY straight.rod DESTINATION ${vdwsiteposdir})
file (COPY straight.rodvdw DESTINATION ${vdwsiteposdir})
file (COPY straight_vdw_pos.csv DESTINATION ${vdwsiteposdir})
file (COPY curvy.rod DESTINATION ${vdwsiteposdir})
file (COPY curvy.rodvdw DESTINATION ${vdwsiteposdir})
file (COPY curvy_vdw_pos.csv DESTINATION ${vdwsiteposdir})
add_test(NAME rod_vdw_site_placement COMMAND "${PROJECT_BINARY_DIR}/src/ffea" rod_vdw_site_placement.ffeatest)
