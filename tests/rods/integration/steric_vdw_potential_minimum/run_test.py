# Test the calculation of the van der Waals forces
# Criteria for pass:
#   Inter-element forces are equal and opposite (first frame)
#   Forces are equal on both nodes of a given element (first frame)
#   Forces match expected values (first frame)
#   After some time, elements settle into infection points:
#       Steric:  rods touching (zero surface-surface distance, or sigma)
#       Interp:  potential well (r_min)
#       LJ 6_12: potential well (r_min)
import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import pprint

from omegaconf import OmegaConf

import ffeatools.ffea_rod as ffea_rod

params = OmegaConf.load("params.yml")


# A bunch of hard-coded potentials, to match stuff from the C++ code. BAD!!
def steric_constant(steric_energy_max, r_max):
    """Hard-coded! Constant defining the strength of repulsive potential"""
    return steric_energy_max / (r_max * r_max)


def lj_6_12(r, r_eq, eps):
    """
        Typical Lennard-Jones potential

        r:     surface-surface distance
        sigma: distance where U = 0
        r_eq:  distance where U = -eps
        eps:   potential well depth

        r > 0 when called

    return U, F
    """
    sigma = r_eq / np.power(2, 1 / 6)
    r_6 = np.power(sigma / r, 6)
    r_12 = np.power(sigma / r, 12)
    return (4 * eps * (r_12 - r_6), -(24 * eps / r) * (2 * r_12 - r_6))


def lj_interp(r, r_eq, eps):
    """
        Interpolative potential, ben Hanson thesis p74-77

        r > 0 when called

    return U, F
    """
    r_2 = np.power(r / r_eq, 2)
    r_3 = np.power(r / r_eq, 3)
    return (eps * (2 * r_3 - 3 * r_2), (6 * eps / r_eq) * (r_2 - r / r_eq))


def steric(r, k):
    """
        Quadratic repulsive potential.

        r < 0 when called

    return U, F
    """
    return (k * r * r, 2 * k * r)


def get_energy_force(r, k_steric, r_eq, eps, steric_func, int_func, lj_func):
    """r = surface-surface distance between rod elements"""

    if r < 0:
        return steric_func(r, k_steric)
    elif r < r_eq and r >= 0:
        return int_func(r, r_eq, eps)
    elif r >= r_eq:
        return lj_func(r, r_eq, eps)
    else:
        raise Exception(f"Out of range, r = {r / r_eq:e} r_eq = {r * 1e10:.3f} A")


def get_max_steric_energy():
    # Get max steric energy directly from rod code
    try:
        ffea_src = os.environ.get("FFEA_SRC", os.environ.get("FFEA_HOME"))
        with open(f"{ffea_src}/include/rod_structure.h", "r") as f:
            for line in f:
                if "max_steric_energy" in line:
                    return float(line.split("=")[-1].split(";")[0])

    except (KeyError, FileNotFoundError):
        return 50


"""
This test expects a single runarg which points to the location of the ffea binary to be tested
e.g. python run_test <path to ffea>
"""
def main():
    # =========================== PREAMBLE =========================== #

    shutil.rmtree("outputs", ignore_errors=True)
    os.makedirs("outputs")

    py_return_status = subprocess.call([sys.executable, "create_straight_rod_vdw.py"])

    if py_return_status != 0:
        sys.exit("Error calling create_straight_rod_vdw.py")

    ffea_scripts = glob.glob("*.ffea")

    if len(ffea_scripts) == 0:
        raise Exception("No FFEA input scripts found!")

    # distance measured between rod element centrelines.
    # also the x translation of the rods in the .ffea files.
    #
    # steric: rods begin overlapped, so they should strongly repel.
    # interp: rods begin in interpolative potential, on the 'left'
    #         side of the LJ potential minimum. They should slightly
    #         repel.
    # LJ:     rods begin in attractive region of LJ potential, on
    #         'right' side of minimum. They should attract.
    centreline_dist = {"steric": 2.5e-9, "interp": 5.5e-9, "lj": 6.5e-9}

    # some hard-coded bits
    force_dim = params["ffea_force"]
    length_dim = params["ffea_length"]
    energy_dim = params["ffea_energy"]

    steric_energy_max = get_max_steric_energy() * energy_dim

    # single float precision
    err_lim = 1e-10
    print(f"Error threshold: {err_lim:e}\n")

    pass_count = 0
    passes = {}

    # =========================== FFEA SCRIPTS =========================== #
    for i, ffea_path in enumerate(ffea_scripts):
        ffea_script_name = ffea_path.split(".")[0]
        regime_name = ffea_script_name.split("_")[0]

        # test function located in: src/ffea_test.cpp
        print("Calling FFEA from Python script")
        ffea_proc = subprocess.run(
            [sys.argv[1], ffea_script_name + ".ffea"], capture_output=True, text=True
        )

        with open("outputs/" + ffea_script_name + "_stdout.txt", "w") as f:
            f.write(ffea_proc.stdout)

        with open("outputs/" + ffea_script_name + "_stderr.txt", "w") as f:
            f.write(ffea_proc.stderr)

        if ffea_proc.returncode != 0:
            sys.exit(f"Error calling FFEA ({ffea_proc.returncode:d})")

        # load both rod trajectories
        print(f"Analysing output of FFEA simulation: {ffea_script_name:s}\n")
        rod0 = ffea_rod.ffea_rod(filename=f"outputs/{regime_name:s}1.rodtraj")
        rod1 = ffea_rod.ffea_rod(filename=f"outputs/{regime_name:s}2.rodtraj")

        # frame index
        fi = 1

        # no force is applied at step 0, which is for initialisation (steric, vdw)
        if regime_name == "steric":
            force0 = rod0.steric_force[fi, :, :]
            force1 = rod1.steric_force[fi, :, :]
            inflection_dist = 0
        elif regime_name == "interp" or regime_name == "lj":
            force0 = rod0.vdw_force[fi, :, :]
            force1 = rod1.vdw_force[fi, :, :]
            inflection_dist = params["r_min"]
        else:
            raise Exception(f"Invalid regime name: {regime_name:s}.")

        # all elements have the same radius
        radius0 = rod0.material_params[1, 0, 2]
        radius1 = rod1.material_params[1, 0, 2]

        # force on nodes of central element
        force_02 = force0[2, :]
        force_03 = force0[3, :]
        force_12 = force1[2, :]
        force_13 = force1[3, :]

        # force occurs halfway along element
        print("# === INTRA-ELEMENT FORCES EQUAL? === #")
        if np.linalg.norm(force_02 - force_03) < err_lim:
            print(" Rod 0, nodes 2 and 3: force magnitudes are equal (PASS)")
            intra_element = True
        else:
            print(" Rod 0, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
            intra_element = False

            print(" node 2: ", force_02)
            print(" node 3: ", force_03)
            print(" delta: ", np.linalg.norm(force_02 - force_03))

        if np.linalg.norm(force_12 - force_13) < err_lim:
            print(" Rod 1, nodes 2 and 3: force magnitudes are equal (PASS)")
            intra_element = True
        else:
            print(" Rod 1, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
            intra_element = False

            print(" node 2: ", force_12)
            print(" node 3: ", force_13)
            print(" delta: ", np.linalg.norm(force_12 - force_13))
        print("")

        print("# === INTER-ELEMENT FORCES EQUAL AND OPPOSITE? === #")
        # force_02 + force_03 = -(force_12 + force_13) sum to zero
        if np.linalg.norm(force_02 + force_03 + force_12 + force_13) < err_lim:
            print(" Inter-element forces are symmetric (PASS)")
            inter_element = True
        else:
            print(" Inter-element forces are NOT symmetric (FAIL)")
            inter_element = False

            print(" rod 0, node 2: ", force_02)
            print(" rod 0, node 3: ", force_03)
            print(" rod 1, node 2: ", force_12)
            print(" rod 1, node 3: ", force_13)
            print(" delta: ", np.linalg.norm(force_02 + force_03 + force_12 + force_13))
        print("")

        print("# === FORCE ONLY ON CENTRAL ELEMENTS? === #")
        force0_outer = np.delete(force0, [6, 7, 8, 9, 10, 11])
        force1_outer = np.delete(force1, [6, 7, 8, 9, 10, 11])
        if np.any(force0_outer != 0):
            print(" Rod 0, force is on NON-central elements (FAIL)")
            print(" offending nodes: ", list(np.where(force0_outer != 0)[0] // 3))
            central_element = False
            print(" Rod 0 forces:")
            print(force0)
        else:
            print(" Rod 0, force is on central element only (PASS)")
            central_element = True

        if np.any(force1_outer != 0):
            print(" Rod 1, force is on NON-central elements (FAIL)")
            print(" offending nodes: ", list(np.where(force1_outer != 0)[0] // 3))
            central_element = False
            print(" Rod 1 forces:")
            print(force1)
        else:
            print(" Rod 1, force is on central element only (PASS)")
            central_element = True
        print("")

        if regime_name == "steric":
            nbrs0 = rod0.get_num_nbrs(1)["steric"]
            nbrs1 = rod1.get_num_nbrs(1)["steric"]
        elif regime_name == "interp" or regime_name == "lj":
            nbrs0 = rod0.get_num_nbrs(1)["vdw"]
            nbrs1 = rod1.get_num_nbrs(1)["vdw"]
        else:
            raise Exception(f"Invalid regime name: {regime_name:s}.")

        print("Rod 0 num neighbours:")
        print(nbrs0)
        print("")

        print("Rod 1 num neighbours:")
        print(nbrs1)
        print("")

        # rods are separated in x only
        radius_sum = radius0 + radius1
        expect_dist = centreline_dist[regime_name] - radius_sum

        # expected energy and force are computed from one of three functions
        expect_energy, expect_force = get_energy_force(
            r=expect_dist,
            k_steric=steric_constant(steric_energy_max, r_max=radius_sum),
            r_eq=params["r_min"],
            eps=params["eps_kT"] * params["ffea_energy"],
            steric_func=steric,
            int_func=lj_interp,
            lj_func=lj_6_12,
        )

        print("# === EXPECTED FORCE ON NODES? === #")
        print(f" Expected surf-surf distance:    {expect_dist:.3e} m")
        print(f" Expected energy:                {expect_energy:.3e} J")
        print(f" Expected element force:         {expect_force:.3e} N")
        print(f" Expected node force:            {expect_force/2:.3e} N\n")

        print(" node forces:")
        print("     02", force_02[0])
        print("     03", force_03[0])
        print("     12", force_12[0])
        print("     13", force_13[0])

        if (
            #sign
            (abs(force_02[0] - expect_force) / 2 < err_lim)
            and (abs(force_03[0] - expect_force) / 2 < err_lim)
            # magnitude
            and (abs(force_02[0]) - abs(expect_force) / 2 < err_lim)
            and (abs(force_03[0]) - abs(expect_force) / 2 < err_lim)
            and (abs(force_12[0]) - abs(expect_force) / 2 < err_lim)
            and (abs(force_13[0]) - abs(expect_force) / 2 < err_lim)
            # non-zero
            and (force_02[0] != 0)
            and (force_03[0] != 0)
            and (force_12[0] != 0)
            and (force_13[0] != 0)
        ):
            match_boe = True
            print(" Node forces match expected values (PASS)")
        else:
            print(" Node forces do not match expected values (FAIL)")
            match_boe = False

        print("")

        print("# === ENERGY MINIMUM? === #")
        # after enough time, the elements should have become separated by r_min, the distance at which the LJ potential is at a minimum
        # NOTE: will this work? won't the internal rod elasticity provide resistance?
        p0 = rod0.get_p_i(rod0.current_r)
        p1 = rod1.get_p_i(rod1.current_r)
        mid0 = rod0.current_r[-1, 2, :] + 0.5 * p0[-1, 2]
        mid1 = rod1.current_r[-1, 2, :] + 0.5 * p1[-1, 2]

        ss_dist = np.linalg.norm(mid0 - mid1) - radius_sum
        delta = ss_dist - inflection_dist

        print("Centroid rod 0:            ", mid0)
        print("Centroid rod 1:            ", mid1)
        print(f"Surface-surface distance:  {ss_dist:.3e}")
        print(f"inflection point:          {inflection_dist:.3e}")
        print(f"delta:                     {delta:.3e}")
        if abs(delta) < err_lim:
            print("Rod elements are separated by r_min (PASS)")
            energy_minimum = True
        else:
            print("Rod elements are not separated by r_min (FAIL)")
            energy_minimum = False
        print("")

        if (
            intra_element
            and inter_element
            and central_element
            and match_boe
            and energy_minimum
        ):
            pass_count += 1
            passes[ffea_script_name] = True
        else:
            print("intra-element forces equal:              ", intra_element)
            print("inter-element forces equal and opposite: ", inter_element)
            print("forces only on central elements:         ", central_element)
            print("forces match expected values:            ", match_boe)
            print("potential energy minimum reached:        ", energy_minimum)
            passes[ffea_script_name] = False
        print("")

    print(f"Pass count: {pass_count:d}")
    pprint.pprint(passes)
    if pass_count == len(ffea_scripts):
        return 0
    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
