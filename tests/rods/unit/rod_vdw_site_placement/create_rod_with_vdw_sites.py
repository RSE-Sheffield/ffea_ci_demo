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
def x_zero(t):
    return 0

def y_zero(t):
    return 0

def z_line(t):
    return t

def spiral_rod(radius, t_min, t_max, theta_min, theta_max, num_nodes):
    """Specific version of create_rod_parametric() from rod_creator.py"""

    t_values = np.linspace(t_min, t_max, num=num_nodes)
    th_values = np.linspace(theta_min, theta_max, num=num_nodes)
    rod = np.empty([num_nodes, 3])
    for i in range(len(t_values)):
        rod[i][0] = radius * np.sin(th_values[i])
        rod[i][1] = radius * np.cos(th_values[i])
        rod[i][2] = t_values[i]

    return rod

def get_parent_element(r_frame, contour_length, norm_length_along_rod, num_nodes):
    """Determine the rod element that a VDW site belongs to."""
    r1 = np.zeros(3)
    r2 = np.zeros(3)
    p = np.zeros(3)
    psum = 0
    psum_prev = 0
    length_along_rod = norm_length_along_rod * contour_length

    # Traverse along the rod
    for node in range(num_nodes - 1):
        r1 = r_frame[node, :]
        r2 = r_frame[node + 1, :]

        p = r2 - r1

        pmag = np.linalg.norm(p)
        psum += pmag
        norm_length_along_elem = (length_along_rod - psum_prev) / pmag
        psum_prev = psum

        # Return if we overstep
        if psum > length_along_rod:
            site_pos = r1 + p * norm_length_along_elem
            return node, norm_length_along_elem, site_pos

    # Otherwise, set to the end of the rod
    return num_nodes - 1, 0, r1 + p


def main():

    # ffea_input_name = params["ffea_input_name"]

    # Blank rods
    rod_line = ffea_rod.ffea_rod(num_elements=params["num_nodes"], num_rods=params["num_rods"], num_vdw_sites=params["num_vdw_sites"])
    nodes_line = rc.create_rod_parametric(
        x_zero, y_zero, z_line, 0, params["length"], params["num_nodes"]
    )

    rod_curve = ffea_rod.ffea_rod(num_elements=params["num_nodes"], num_rods=params["num_rods"], num_vdw_sites=params["num_vdw_sites"])

    nodes_curve = spiral_rod(
        radius = params["length"]/5, t_min = 0, t_max = params["length"], theta_min = 0, theta_max = 2*np.pi, num_nodes = params["num_nodes"]
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
        # rc.write_lj_matrix(ffea_input_name, lj_params)

        sites = []
        for x in np.linspace(0, 1, num=params["num_vdw_sites"]):
            sites.append((0, x))

        print(f"Number of VDW sites: {params['num_vdw_sites']:d}")
        print("Contents of .rodvdw file. List of (face, pos):")
        pp.pprint(sites, width=30)
        print("Contents of .rodlj file. List of (face_i, face_j, epsilon_ij, sigma_ij)")
        pp.pprint(lj_params, width=60)

        # Types and positions of interaction sites
        rc.write_vdw_sites(rod=my_rod, vdw_fname=f"{rod_name:s}.rodvdw", vdw_sites=sites)

        my_rod.write_rod(f"{params['in_dir']}{rod_name:s}.rod")

        site_pos = np.zeros((params['num_vdw_sites'], 3))
        for i, (face, L0) in enumerate(sites):
            elem_id, L_elem, site_pos[i, :] = get_parent_element(my_rod.equil_r[0], my_rod.get_contour_length(), L0, params["num_nodes"])

        with open (f"{rod_name}_vdw_pos.csv", "w") as f:
            for pos in site_pos:
                f.write(f"{pos[0]:6e},{pos[1]:6e},{pos[2]:6e}\n")

        # rc.preview_rod(my_rod.current_r[0], vdw_site=site_pos, scale=1e9, save=True)

if __name__ == "__main__":
    main()
