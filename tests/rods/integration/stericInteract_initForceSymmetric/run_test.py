# Test the calculation and interpolation of steric repulsion forces
# Criteria for pass:
#   Intra-element forces are equal
#   Force acts only on central elements
#   Inter-element forces are equal and opposite
#   Forces match expected values
import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import re

import ffeatools.ffea_rod as ffea_rod

def steric_potential(k, x):
    """Hard-coded! Square repulsive potential"""
    return k * x * x

def steric_force(k, x):
    """Hard-coded! Force of square repulsive potential"""
    return 2 * k * x

def steric_constant(U_max, x_max):
    """Hard-coded! Constant defining the strength of repulsive potential"""
    return U_max / (x_max*x_max)


"""
This test expects a single runarg which points to the location of the ffea binary to be tested
e.g. python run_test <path to ffea>
"""
def main():

    shutil.rmtree("outputs", ignore_errors=True)
    os.makedirs("outputs")

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
    rod0 = ffea_rod.ffea_rod(filename="outputs/1.rodtraj")
    rod1 = ffea_rod.ffea_rod(filename="outputs/2.rodtraj")

    # some hard-coded bits
    force_dim = 2.4364705882352941e-11
    length_dim = 1.7e-10
    energy_dim = 4.141945559999999e-21

    print(f"FFEA force unit: {force_dim:e}")
    print(f"FFEA length unit: {length_dim:e}")
    print(f"FFEA energy unit: {energy_dim:e}")

    # no force is applied at step 0, which is for initialisation
    force0 = rod0.steric_force[1, :, :]
    force1 = rod1.steric_force[1, :, :]

    # all elements have the same radius
    radius0 = rod0.material_params[1, 0, 2]
    radius1 = rod1.material_params[1, 0, 2]

    err_lim = 1e-12
    print(f"Error threshold: {err_lim:e}\n")

    # force on nodes of central element
    force_02 = force0[2, :]
    force_03 = force0[3, :]
    force_12 = force1[2, :]
    force_13 = force1[3, :]

    # force occurs halfway along element
    print("# === INTRA-ELEMENT FORCES EQUAL? === #")
    if np.linalg.norm(force_02 - force_03) < err_lim:
        print("Rod 0, nodes 2 and 3: force magnitudes are equal (PASS)")
        intra_element = True
    else:
        print("Rod 0, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
        intra_element = False

    print("node 2: ", force_02)
    print("node 3: ", force_03)
    print("delta: ", np.linalg.norm(force_02 - force_03))
    print("")

    if np.linalg.norm(force_12 - force_13) < err_lim:
        print("Rod 1, nodes 2 and 3: force magnitudes are equal (PASS)")
        intra_element = True
    else:
        print("Rod 1, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
        intra_element = False

    print("node 2: ", force_12)
    print("node 3: ", force_13)
    print("delta: ", np.linalg.norm(force_12 - force_13))
    print("")

    print("# === INTER-ELEMENT FORCES SYMMETRIC? === #")
    # force_02 + force_03 = -(force_12 + force_13)
    if np.linalg.norm(force_02 + force_03 + force_12 + force_13) < err_lim:
        print("Inter-element forces are symmetric (PASS)")
        inter_element = True
    else:
        print("Inter-element forces are NOT symmetric (FAIL)")
        inter_element = False

    print("rod 0, node 2: ", force_02)
    print("rod 0, node 3: ", force_03)
    print("rod 1, node 2: ", force_12)
    print("rod 1, node 3: ", force_13)
    print("delta: ", np.linalg.norm(force_02 + force_03 + force_12 + force_13))
    print("")

    print("# === FORCE ONLY ON CENTRAL ELEMENTS? === #")
    force0_outer = np.delete(force0, [6, 7, 8, 9, 10, 11])
    force1_outer = np.delete(force1, [6, 7, 8, 9, 10, 11])
    if np.any(force0_outer != 0):
        print("Rod 0, steric force is on NON-central elements (FAIL)")
        print("offending nodes: ", list(np.where(force0_outer != 0)[0] // 3))
        central_element = False
    else:
        print("Rod 0, steric force is on central element only (PASS)")
        central_element = True

    print("Rod 0 forces:")
    print(force0)
    print("")

    if np.any(force1_outer != 0):
        print("Rod 1, steric force is on NON-central elements (FAIL)")
        print("offending nodes: ", list(np.where(force1_outer != 0)[0] // 3))
        central_element = False
    else:
        print("Rod 1, steric force is on central element only (PASS)")
        central_element = True

    print("Rod 1 forces:")
    print(force1)
    print("")

    print("Rod 0 num neighbours:")
    print(rod0.get_num_nbrs(1)["steric"])
    print("")

    print("Rod 1 num neighbours:")
    print(rod1.get_num_nbrs(1)["steric"])
    print("")

    # Get max steric energy directly from rod code
    try:
        ffea_src = os.environ.get("FFEA_SRC", os.environ.get("FFEA_HOME"))
        with open(f"{ffea_src}/include/rod_structure.h", "r") as f:
            for line in f:
                if "max_steric_energy" in line:
                    U_max = float(line.split("=")[-1].split(";")[0])
                    break
    except (KeyError, FileNotFoundError):
        U_max = 50

    U_max *= energy_dim

    # rods are initially separated in x only
    try:
        ss_dist = abs(rod0.current_r[1, 0, 0] - rod1.current_r[1, 0, 0]) - (radius0+radius1)
    except:
        ss_dist = -2.5e-9

    k_steric = steric_constant(U_max, x_max = radius0 + radius1)
    boe_energy = steric_potential(k_steric, ss_dist)
    boe_force = steric_force(k_steric, ss_dist)
    print(f"Surf-surf distance in x:    {ss_dist:.3e} m")
    print(f"Maximum overlap:")
    print(f"  Max steric energy:        {U_max:.3e} J")
    print(f"  Radius 0:                 {radius0:.3e} m")
    print(f"  Radius 1:                 {radius1:.3e} m")
    print(f"  Steric force constant:    {k_steric:.3e} N/m")
    print("Expected values:")
    print(f"  Expected energy:          {boe_energy:.3e} J")
    print(f"  Expected element force:   {boe_force:.3e} N")
    print(f"  Expected node force:      {boe_force/2:.3e} N")

    if (force_02[0] - boe_force/2 < err_lim) and (force_03[0] - boe_force/2) and (force_12[0] - boe_force/2) and (force_13[0] - boe_force/2):
        match_boe = True
        print("Node forces match expected values (PASS)")
    else:
        print("Node forces do not match expected values (FAIL)")
        print("deltas:")
        print(force_02[0] - boe_force/2)
        print(force_03[0] - boe_force/2)
        print(force_12[0] - boe_force/2)
        print(force_13[0] - boe_force/2)
        match_boe = False

    if intra_element and inter_element and central_element and match_boe:
        return 0
    else:
        print("intra-element forces equal: ", intra_element)
        print("inter-element forces symmetric: ", inter_element)
        print("forces only on central elements: ", central_element)
        print("forces match BoE calculation: ", match_boe)

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
