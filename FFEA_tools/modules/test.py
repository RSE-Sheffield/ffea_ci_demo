import FFEA_trajectory
import sys
from Vectors import *

mytraj = FFEA_trajectory.FFEA_trajectory(sys.argv[1], 4000)
mytraj.blob[1].define_sub_blob_radially(0, vector3(0.0, 0.0, 0.0), 250e-10, 1000e-10)
mytraj.blob[1].calc_sub_blob_centroid_moment_about_x0(2, 0)
