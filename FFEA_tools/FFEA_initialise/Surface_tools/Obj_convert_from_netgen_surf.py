import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT SURF FNAME] [OUTPUT OBJ FNAME]")

# Get surf files
surf = open(sys.argv[1], "r")
obj = open(sys.argv[2], "w")

# Parse it and write to out
surf.readline()
num_nodes = int(surf.readline())
for i in range(num_nodes):
	obj.write("v " + surf.readline().strip() + "\n")

num_faces = int(surf.readline())
for i in range(num_faces):
	obj.write("f " + surf.readline().strip() + "\n")

surf.close()
obj.close()

