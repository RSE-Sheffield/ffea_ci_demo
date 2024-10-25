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

import numpy as np
import sys

ffeatoolsFound = False
try:
    import ffeatools # python package
    ffeatoolsFound = True
    FFEA_script = ffeatools.modules.FFEA_script.FFEA_script
except:
    try:
        from ffeatools.ffea_script import ffea_script as FFEA_script
    except ImportError:
        print("Failure to import relevant FFEA modules")
        sys.exit(1) # failure to import

# Load trajectory
start = 5000
#start = 0
end = 10000
try:
	script = FFEA_script(sys.argv[1])
except:
	sys.exit("Script please!")
traj = script.load_trajectory(start=start, num_frames = end-start)
end = traj.num_frames + start
mat = script.load_material(0)
top = script.load_topology(0)

# Analyse trajctory in sets of 1ps and test against theoretical diffusion
r2 = [None for i in range(start, end-1)]

for f in traj.blob[0][0].frame:
	f.calc_centroid()

for i in range((end - start) - 1):
	r2[i] = traj.blob[0][0].frame[i + 1].get_centroid() - traj.blob[0][0].frame[i].get_centroid()
	r2[i] = np.dot(r2[i], r2[i])

r2mean = np.mean(r2, axis=0)
r2stdev = np.std(r2, axis=0)
r2err = r2stdev / np.sqrt(len(r2))

mass = traj.blob[0][0].frame[0].calc_mass(top, mat)
r2theory = 3 * (script.params.kT / mass) * ((script.params.dt * script.params.check)**2)
#print "Calculated Diffusion: <r^2> = %6.2e +/- %6.2e" % (r2mean, r2err)
#print "Theoretical Diffusion: <r^2> = %6.2e" % (r2theory)
#print "Percent Error: d<r^2>/<r^2> = %6.2f%%" % (100 * np.fabs(r2mean - r2theory) / r2theory)
err = 0
tolerance = 0.1
if np.fabs(r2mean - r2theory) / r2theory <= tolerance:
	print("Inertial Diffusion", ": correct ", np.fabs(r2mean - r2theory) / r2theory , " < ", tolerance)
else:
	print("Inertial Diffusion", ": failed ", np.fabs(r2mean - r2theory) / r2theory , " > ", tolerance)
	err = 1
sys.exit(err)
