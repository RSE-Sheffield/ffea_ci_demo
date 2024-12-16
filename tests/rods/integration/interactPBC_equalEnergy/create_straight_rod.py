
try:
    # ffeatools 2 (python 3)
    import ffeatools.ffea_rod as ffea_rod
    from ffeatools.rod.rod_creator import rod_creator as rc
except ModuleNotFoundError:
    # ffeatools 1 (python 2.7)
    import ffeatools.modules.FFEA_rod as ffea_rod
    from ffeatools.modules.FFEA_rod import rod_creator as rc


# Rod is 5 nm x 50 nm
num_nodes = 6
radius = 2.5e-9  # m
length = 5e-8

stretch = 3.5e-11  # N
twist = 5e-29  # N.m^2
bend = 3.5e-29  # m^4.Pa

# Straight line in z-axis
def x_func(t):
    return 0

def y_func(t):
    return 0

def z_func(t):
    return t

def main():

    # Blank rod
    try:
        my_rod = ffea_rod.ffea_rod(num_elements=num_nodes)
    except AttributeError:
        my_rod = ffea_rod.FFEA_rod(num_elements=num_nodes)

    my_nodes = rc.create_rod_parametric(x_func, y_func, z_func, 0, length, num_nodes)

    my_rod.current_r[0] = my_nodes
    my_rod.equil_r[0] = my_nodes

    my_rod.current_m[0], my_rod.equil_m[0] = rc.create_material_frame(my_rod)

    rc.set_params(rod=my_rod, stretch_constant=stretch,
                torsion_constant=twist, radius=radius,
                bending_modulus=bend)

    my_rod.write_rod("z-axis.rod")

if __name__ == "__main__":
    main()