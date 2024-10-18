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
    rod1a = ffea_rod.ffea_rod(filename="1A.rodtraj")
    rod1b = ffea_rod.ffea_rod(filename="1B.rodtraj")
    rod2a = ffea_rod.ffea_rod(filename="2A.rodtraj")
    rod2b = ffea_rod.ffea_rod(filename="2B.rodtraj")

    r = {
        "1a" : rod1a.current_r,
        "1b" : rod1b.current_r,
        "2a" : rod2a.current_r,
        "2b" : rod2b.current_r
    }

    energy = {
        "1a" : rod1a.steric_energy,
        "1b" : rod1b.steric_energy,
        "2a" : rod2a.steric_energy,
        "2b" : rod2b.steric_energy
    }

    ax = plt.figure().add_subplot(projection='3d')
    ax.legend()
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_zlabel('z (nm)')

    for rod in r.keys():
        ax.plot(r[rod][1, :, 0]*1e9, r[rod][1, :, 1]*1e9, r[rod][1, :, 2]*1e9, '-', label=rod)

    ax.view_init(elev=20., azim=-35, roll=0)
    plt.legend()
    plt.savefig("position_preview.png", dpi=400)
    plt.close()

    # Check that all rods have a finite interaction energy on at least one node
    no_zero_energy = energy["1a"][-1].any() and energy["1b"][-1].any() and energy["2a"][-1].any() and energy["2b"][-1].any()

    # Surface-surface energy 1A and 2A should be identical, as should 1B and 2B
    diff_energy_ab = (energy["1a"][-1] - energy["2a"][-1]) + (energy["1b"][-1] - energy["2b"][-1])

    print("ENERGIES (t=tmax):")
    print("1A (centre)\n", energy["1a"][-1])
    print("1B (centre)\n", energy["1b"][-1])
    print("2A (edge)\n", energy["2a"][-1])
    print("2B (edge)\n", energy["2b"][-1])
    print("diff\n", diff_energy_ab)
    print(f"<diff>: {np.mean(diff_energy_ab):e}")
    print(f"no zero energy rods: {no_zero_energy}\n")

    # Using absolute node positions, the central rods should move further apart,
    # but the edge rods should move closer together, since the latter will repel
    # from their periodic images.
    r_diff_cent = [r["1b"][0] - r["1a"][0], r["1b"][-1] - r["1a"][-1]]
    r_diff_edge = [r["2b"][0] - r["2a"][0], r["2b"][-1] - r["2a"][-1]]

    # only need to consider displacements in the y-axis
    diff_y_cent = r_diff_cent[1][:, 1] - r_diff_cent[0][:, 1]
    diff_y_edge = r_diff_edge[1][:, 1] - r_diff_edge[0][:, 1]

    print("NODE-NODE DISPLACEMENTS:")
    print("Central (t=0)\n", r_diff_cent[0])
    print("Central (t=tmax)\n", r_diff_cent[1])
    print(f"<diff_y>: {np.mean(diff_y_cent):e}")
    print("Edge (t=0)\n", r_diff_edge[0])
    print("Edge (t=tmax)\n", r_diff_edge[1])
    print(f"<diff_y>: {np.mean(diff_y_edge):e}")

    plt.title("Node 3 displacement from initial position")
    nsteps = r["1a"].shape[0]
    plt.plot(np.arange(nsteps), np.linalg.norm((r["1a"][:, 2, :] - r["1a"][0, 2, :])*1e9, axis=1), '-', label="Rod 1a")
    plt.plot(np.arange(nsteps), np.linalg.norm((r["1b"][:, 2, :] - r["1b"][0, 2, :])*1e9, axis=1), '-', label="Rod 1b")
    plt.plot(np.arange(nsteps), np.linalg.norm((r["2a"][:, 2, :] - r["2a"][0, 2, :])*1e9, axis=1), '-', label="Rod 2a")
    plt.plot(np.arange(nsteps), np.linalg.norm((r["2b"][:, 2, :] - r["2b"][0, 2, :])*1e9, axis=1), '-', label="Rod 2b")
    plt.ylabel("Displacement (nm)")
    plt.xlabel("Step")
    plt.legend()
    plt.savefig("Node3Pos.png", dpi=300)

    print("Pass conditions:")
    print("  Central rods repel towards edge: ", np.mean(diff_y_cent) > 0)
    print("  Edge rods repel towards centre:  ", np.mean(diff_y_edge) < 0)
    print("  Same energy on all rods:         ", abs(np.mean(diff_energy_ab)) < 1e-12)
    print("  No zero energy rods:             ", no_zero_energy)

    if np.mean(diff_y_cent) > 0 and np.mean(diff_y_edge) < 0 and abs(np.mean(diff_energy_ab)) < 1e-12 and no_zero_energy:
        return 0

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
