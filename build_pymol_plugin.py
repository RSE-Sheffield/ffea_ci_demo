import os
import tarfile

software_dir = os.environ['HOME'] + "/Software"
ffeatools_dir = os.environ['FFEA_SRC'] + "/ffeatools"
tar_dir = os.environ['FFEA_SRC']  + "/../FFEAplugin"
ffea_build_pymol = os.environ['FFEA_BUILD'] + "/share/ffea/plugins/pymol"
ffea_debug_pymol = os.environ['FFEA_BUILD'] + "/share/ffea/plugins/pymol"
pymol_files =  ["FFEA_analysis/pymol_plugin/Blob.py",
    "FFEA_analysis/pymol_plugin/__init__.py",
    "FFEA_analysis/pymol_plugin/mtTkinter.py",
    "FFEA_analysis/pymol_plugin/system_names_pokemon.txt",
    "FFEA_analysis/pymol_plugin/system_names_dbzcharacters.txt",
    "FFEA_analysis/pymol_plugin/system_names_greekletters.txt",
    "modules/FFEA_trajectory.py", 
    "modules/FFEA_turbotrajectory.py",
    "modules/FFEA_measurement.py", 
    "modules/FFEA_frame.py",
    "modules/FFEA_universe.py",
    "modules/FFEA_binding_sites.py", 
    "modules/FFEA_script.py", 
    "modules/FFEA_pin.py", 
    "modules/FFEA_springs.py", 
    "modules/FFEA_material.py", 
    "modules/FFEA_topology.py", 
    "modules/FFEA_surface.py",
    "modules/FFEA_node.py", 
    "modules/FFEA_stokes.py", 
    "modules/FFEA_vdw.py", 
    "modules/FFEA_lj.py", 
    "modules/FFEA_pdb.py", 
    "modules/FFEA_beads.py", 
    "modules/FFEA_ctforces.py", 
    "modules/FFEA_exceptions.py", 
    "modules/FFEA_kinetic_map.py",
    "modules/FFEA_kinetic_rates.py",
    "modules/FFEA_kinetic_states.py",
    "modules/FFEA_rod.py",
    "modules/FFEA_skeleton.py"]

try:
    os.mkdir(tar_dir)
except OSError:
    pass
for f in pymol_files:
    os.system("cp -v " + ffeatools_dir + "/" + f + " " + tar_dir + "/")

os.system("tar cfzv FFEAplugin.tar.gz FFEAplugin/")
os.system("cp FFEAplugin.tar.gz " + ffea_build_pymol + "/")
os.system("cp FFEAplugin.tar.gz " + ffea_debug_pymol + "/")
os.system("rm -r FFEAplugin.tar.gz FFEAplugin")
