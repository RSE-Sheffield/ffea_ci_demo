# Create a FFEA / KOBRA rod that is straight in the z-axis.
from omegaconf import OmegaConf
import ffeatools.ffea_rod as ffea_rod
from ffeatools.rod.rod_creator import rod_creator as rc
import pprint as pp
import sys
import numpy as np

# Default rod parameters
params = OmegaConf.load("params.yml")
stretch = 3.5e-11  # N
twist = 5e-29  # N.m^2
bend = 3.5e-29  # m^4.Pa
diameter = 2 * params["radius"]

# Straight line in z-axis
def x_line(t):
    return 0

def y_line(t):
    return 0

def z_line(t):
    return t

def x_curve(t):
    return np.sin(t)

def y_curve(t):
    return np.cos(t)

def main():

    ffea_input_name = params["ffea_input_name"]

    # Blank rods
    rod_line = ffea_rod.ffea_rod(num_elements=params["num_nodes"], num_rods=params["num_rods"], num_vdw_sites=params["num_vdw_sites"])
    nodes_line = rc.create_rod_parametric(
        x_line, y_line, z_line, 0, params["length"], params["num_nodes"]
    )

    rod_curve = ffea_rod.ffea_rod(num_elements=params["num_nodes"], num_rods=params["num_rods"], num_vdw_sites=params["num_vdw_sites"])
    nodes_curve = rc.create_rod_parametric(
        x_curve, y_curve, z_line, 0, params["length"], params["num_nodes"]
    )

    for my_rod, my_nodes, rod_name in zip([rod_line, rod_curve], [nodes_line, nodes_curve], ["straight", "curvy"]):

        my_rod.current_r[0] = my_nodes
        my_rod.equil_r[0] = my_nodes

        my_rod.current_m[0], my_rod.equil_m[0] = rc.create_material_frame(my_rod)

        rc.set_params(
            rod=my_rod,
            stretch_constant=stretch,
            torsion_constant=twist,
            radius=params["radius"],
            bending_modulus=bend,
        )

        # LJ parameter matrix (SI units)
        sig = float(params["r_min"]) / np.power(2, 1.0/6)
        lj_params = [(0, 0, 10 * float(params["kT"]), sig)]
        rc.write_lj_matrix(ffea_input_name, lj_params)

        sites = []
        for x in np.linspace(0, 1, num=params["num_vdw_sites"]):
            sites.append((0, x))

        print(f"Number of VDW sites: {params['num_vdw_sites']:d}")
        print("Contents of .rodvdw file. List of (face, pos):")
        pp.pprint(sites, width=30)
        print("Contents of .rodlj file. List of (face_i, face_j, epsilon_ij, sigma_ij)")
        pp.pprint(lj_params, width=60)

        # Types and positions of interaction sites
        for i in range(params["num_rods"]):
            rc.write_vdw_sites(rod=my_rod, vdw_fname=f"{rod_name:s}_{i+1:d}.rodvdw", vdw_sites=sites)

        my_rod.write_rod(f"{params['in_dir']}{rod_name:s}.rod")

if __name__ == "__main__":
    main()
