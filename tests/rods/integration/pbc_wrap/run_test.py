import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import ffeatools.ffea_rod as ffea_rod

import matplotlib as mpl
mpl.use('Agg')

"""
This test expects a single runarg which points to the location of the ffea binary to be tested
e.g. python run_test <path to ffea>
"""
def main():

    subprocess.call(["bash", "cleanup.sh"])
    py_return_status = subprocess.call([sys.executable, "create_straight_rod.py"])

    if py_return_status != 0:
        sys.exit("Error calling create_straight_rod.py")

    if len(glob.glob("*.ffea")) > 1:
        sys.exit("Error: more than one .ffea file - cannot determine test "
                 "name")
    else:
        ffea_script_name = glob.glob("*.ffea")[0].split(".")[0]

    # test function located in: src/ffea_test.cpp
    print("Calling FFEA from Python script")
    subprocess.call([sys.argv[1], ffea_script_name + ".ffea"])

    # load both rod trajectories
    print("Analysing output of FFEA simulation" + str("\n"))
    rod1 = ffea_rod.ffea_rod(filename="1.rodtraj")

    centre_node_r = rod1.current_r[:, 2, :]*1e9

    box_dim_nm = 1176.47*1.7e-10*1e9
    steps = np.arange(0,centre_node_r.shape[0])

    plt.plot(steps, centre_node_r[:,0],'b-', label="x")
    plt.plot(steps, centre_node_r[:,1],'r-', label="y")
    plt.plot(steps, centre_node_r[:,2],'g-', label="z")
    plt.plot(steps, np.zeros(steps.size), 'k--')
    plt.plot(steps, np.ones(steps.size)*box_dim_nm, 'k--')
    plt.xlabel("t (s)")
    plt.ylabel("y (nm)")
    plt.legend()
    plt.savefig("xyzt_2d.png",dpi=400)

    print(f"box_dim     : {box_dim_nm:e} nm")
    print(f"y(step=95)  : {centre_node_r[95, 1]:e} nm")
    print(f"y(step=105) : {centre_node_r[105, 1]:e} nm")
    # a few steps either side of the crossing point (+ve y wall)
    if (box_dim_nm/2 < centre_node_r[95, 1] < box_dim_nm) and (0 < centre_node_r[105, 1] < box_dim_nm/2):
        return 0
    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
