
try:
    # ffeatools 1 (python 2.7)
    import ffeatools.modules.FFEA_rod as ffea_rod
    from ffeatools.modules.FFEA_rod import rod_creator as rc
except ModuleNotFoundError:
    # ffeatools 2 (python 3)
    import ffeatools.ffea_rod as ffea_rod
    from ffeatools.rod.rod_creator import rod_creator as rc

import pprint as pp
import numpy as np
from omegaconf import OmegaConf
import glob

# Default rod parameters
params = OmegaConf.load("params.yml")
num_nodes = params["num_nodes"]
stretch = 3.5e-11  # N
twist = 5e-29  # N.m^2
bend = 3.5e-29  # m^4.Pa
radius = params["radius"]  # m
diameter = 2 * params["radius"]
length = params["length"]

# Straight line in z-axis
def x_func(t):
    return 0

def y_func(t):
    return 0

def z_func(t):
    return t

def main():

    rod_name = "z-axis"

    # only one VDW site: in the middle of the rod
    sites = [(0, 0.5)]

    # Blank rod
    my_rod = ffea_rod.ffea_rod(num_elements=num_nodes, num_vdw_sites=len(sites))

    my_nodes = rc.create_rod_parametric(x_func, y_func, z_func, 0, length, num_nodes)

    my_rod.current_r[0] = my_nodes
    my_rod.equil_r[0] = my_nodes

    my_rod.current_m[0], my_rod.equil_m[0] = rc.create_material_frame(my_rod)

    rc.set_params(rod=my_rod, stretch_constant=stretch,
                torsion_constant=twist, radius=radius,
                bending_modulus=bend)

    # LJ parameter matrix (SI units)
    sig = float(params["r_min"]) / np.power(2, 1.0/6)
    eps = params["eps_kT"] * params["ffea_energy"]
    lj_params = [(0, 0, eps, sig)]

    # each .ffea script can use the same .rodlj file
    rc.write_lj_matrix("example", lj_params)

    # Types and positions of interaction sites
    rc.write_vdw_sites(rod=my_rod, vdw_fname=f"{rod_name:s}.rodvdw", vdw_sites=sites)

    my_rod.write_rod(f"{rod_name:s}.rod")

    print(f"Number of VDW sites: {len(sites):d}")
    print("Contents of .rodvdw file. List of (face, pos):")
    pp.pprint(sites, width=30)
    print("Contents of .rodlj file. List of (face_i, face_j, epsilon_ij, sigma_ij)")
    pp.pprint(lj_params, width=60)

if __name__ == "__main__":
    main()