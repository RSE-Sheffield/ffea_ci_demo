import sys
from math import *

if len(sys.argv) != 4:
	sys.exit("Usage: python transform_vol.py [INPUT VOL FILE] [INPUT PDB FILE] [OUTPUT .VOL FILE]")

inputvol = open(sys.argv[1], "r")
inputpdb = open(sys.argv[2], "r")
outputvol = open(sys.argv[3], "w")

#Calculating Centroid of pdb

pdbcentroid_x = 0.0
pdbcentroid_y = 0.0
pdbcentroid_z = 0.0
num_atoms = 0

pdblines = inputpdb.readlines()
for line in pdblines:
	if line.split()[0] == "ATOM":
		num_atoms = num_atoms + 1
		pdbcentroid_x = pdbcentroid_x + float(line[30:38])
		pdbcentroid_y = pdbcentroid_y + float(line[38:46])
		pdbcentroid_z = pdbcentroid_z + float(line[46:54])

print "Calculating Centoid for " + str(num_atoms) + " atoms\n"
pdbcentroid_x = pdbcentroid_x/num_atoms
pdbcentroid_y = pdbcentroid_y/num_atoms
pdbcentroid_z = pdbcentroid_z/num_atoms
print ".pdb Centroid is (" + str(pdbcentroid_x) + ", " + str(pdbcentroid_y) + ", " + str(pdbcentroid_z) + "\n"

#Calculating Centroid of vol

volcentroid_x = 0.0
volcentroid_y = 0.0
volcentroid_z = 0.0
num_nodes = 0

vollines = inputvol.readlines()
begin = vollines.index("points\n") + 2
for line in vollines:
	if vollines.index(line) >= begin:
		num_nodes = num_nodes + 1
		sline = line.split()
		volcentroid_x =	volcentroid_x + float(sline[0])	
		volcentroid_y =	volcentroid_y + float(sline[1])
		volcentroid_z =	volcentroid_z + float(sline[2])

if num_nodes == int(vollines[begin - 1]):
	print "Calculating Centroid for " + str(num_nodes) + " nodes\n"
	volcentroid_x =	volcentroid_x /num_nodes
	volcentroid_y =	volcentroid_y /num_nodes
	volcentroid_z =	volcentroid_z /num_nodes
	print ".vol Centroid is (" + str(volcentroid_x) + ", " + str(volcentroid_y) + ", " + str(volcentroid_z) + "\n"
else:
	sys.exit("Calculated Centroid for " + str(num_nodes) + "nodes. Actually we have " + vollines[begin - 1] + " nodes\n")

#Transform vector
translate_x = pdbcentroid_x - volcentroid_x	
translate_y = pdbcentroid_y - volcentroid_y
translate_z = pdbcentroid_z - volcentroid_z
print "Translating .vol by (" + str(translate_x) + ", " + str(translate_y) + ", " + str(translate_z) + "\n"

#Printing to output vol
for line in vollines:
	if vollines.index(line) < begin:
		outputvol.write(line)
	else:
		sline = line.split()
		new_line = str(float(sline[0]) + translate_x).rjust(13) + "  " + str(float(sline[1]) + translate_y).rjust(13) + "  " + str(float(sline[2]) + translate_z).rjust(13) + "\n" 
		outputvol.write(new_line)
		
