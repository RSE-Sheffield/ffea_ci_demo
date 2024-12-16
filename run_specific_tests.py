# Run only the tests I have written
import os

ffea_build_dir = os.environ["FFEA_BUILD"]
tests = ["steric_energy_two_rod_elements",
        "rod_neighbour_list_construction",
        "line_connecting_rod_elements",
        "point_lies_within_rod_element"]

os.chdir(ffea_build_dir)
for t in tests:
    os.system("ctest --verbose -R " + str(t))
    print("")
