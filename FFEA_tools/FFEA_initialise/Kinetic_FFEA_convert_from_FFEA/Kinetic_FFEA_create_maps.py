import sys, os

# Function for getting rid of all files created by this program
def delete_files():
	
	os.system("rm " + working_spring_fname + " " + working_script_fname + " " + "working.pin" + " " + "final0.node" + " " + "final1.node")

# Function for getting files from script
def get_structure_from_script(fname):
	
	# Open script and get lines
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()

	# Get structure files
	node_fname = []
	top_fname = []
	for line in lines:
		if "=" not in line:
			continue
		lvalue, rvalue = line.strip()[1:-1].split("=")
		if lvalue.strip() == "nodes":
			node_fname.append(rvalue.strip())
		elif lvalue.strip() == "topology":
			top_fname.append(rvalue.strip())

	# Quit if not enough files
	if len(node_fname) != 2 or len(top_fname) != 2:
		delete_files()
		sys.exit("Error. '" + fname + "' does not contain two blobs as required.")

	return node_fname[0], node_fname[1], top_fname[0], top_fname[1]

def get_parameters_from_script(fname, param_lvalue):
	
	# Open script and get lines
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get values
	parameters = []
	for line in lines:
		if "=" not in line:
			continue
		lvalue, rvalue = line.strip()[1:-1].split("=")
		if lvalue.strip() == param_lvalue:
			parameters.append(rvalue.strip())
	
	return parameters

# Main program begins here
if len(sys.argv) < 4:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .ffea script (must have 2 blobs)] [spring_constant] [tolerance length] [output check rate (num_frames)] [Pairs of nodes {[i,j]}]\n")

# Get initial scripts
script_fname = sys.argv[1]
in_script_base, in_script_ext = os.path.splitext(os.path.abspath(script_fname))

k = float(sys.argv[2])
l = float(sys.argv[3])
check = int(sys.argv[4])
node_pairs = []
node_pair_string = ""
for i in range(5, len(sys.argv), 1):
	try:
		node_pairs.append([int(sys.argv[i].strip()[1:-1].split(",")[0]), int(sys.argv[i].strip()[1:-1].split(",")[1])])
		node_pair_string += sys.argv[i].strip()[1:-1].split(",")[0] + " " + sys.argv[i].strip()[1:-1].split(",")[1] + " "
	except IndexError:
		sys.exit("Error. '" + sys.argv[i] + "' improperly formatted. Node pairs should look like this: [i,j]")
	except ValueError:
		sys.exit("Error. '" + sys.argv[i] + "' contains a non integer character. Or improperly formatted. Node pairs should look like this: [i,j]")

# Get node and topology fnames
orig_node_fname = ["",""]
orig_top_fname = ["",""]
orig_node_fname[0], orig_node_fname[1], orig_top_fname[0], orig_top_fname[1] = get_structure_from_script(script_fname)

# Make spring file
working_spring_fname = in_script_base + ".spring"
while os.path.isfile(working_spring_fname):
	working_spring_fname = working_spring_fname[0:-7] + "_a" + ".spring"

os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)

# Make working script
working_script_fname = in_script_base + "_working.ffea"
while os.path.isfile(working_script_fname):
	working_script_fname = working_script_fname[0:-5] + "_a" + ".ffea"
os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_script.py " + script_fname + " " + working_script_fname + " " + working_spring_fname + " " + "0" + " " + str(check))

# Run working script repeatadly until user is satisfied
run = 0
while True:

	# Make restart after first go
	if run > 0:
		os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_script.py " + script_fname + " " + working_script_fname + " " + working_spring_fname + " " + str(run) + " " + str(check))

	os.system("ffea " + working_script_fname)
	run += 1
	finished = raw_input("Are you satisfied with this level of overlap (y/n)? (check '" + os.path.basename(working_script_fname) + "' in the viewer):")
	if finished.lower() == "y":
		break
	
	extra_spring = raw_input("Would you like to add another spring (y/n)?:")
	if extra_spring.lower() == "y":
		extra_node_pair = []
		extra_node_pair.append(int(raw_input("Enter node of first blob to add a spring to:")))
		extra_node_pair.append(int(raw_input("Enter node of second blob to connect this spring to:")))
		node_pairs.append(extra_node_pair)
		node_pair_string += str(extra_node_pair[0]) + " " + str(extra_node_pair[1]) + " "
	
		os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)
		continue

	delete_spring = raw_input("Would you like to delete the last spring (y/n)?:")
	if delete_spring.lower() == "y":
		node_pairs.pop()
		node_pair_string = node_pair_string[0:-4]
		os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)
		continue

# Satisfied! Now make two node files and create two maps!
final_node_fname = ["final0.node", "final1.node"]
working_traj_fname = os.path.abspath(get_parameters_from_script(working_script_fname, "trajectory_out_fname")[0])
scale = [0.0,0.0]
scale[0], scale[1] = get_parameters_from_script(working_script_fname, "scale")
os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_analysis/Traj_tools/FFEA_traj_to_nodes.py " + working_traj_fname + " " + "1000" + " " + orig_node_fname[0] + " " + final_node_fname[0] + " " + scale[0] + " " + orig_node_fname[1] + " " + final_node_fname[1] + " " + scale[1])
os.system("/localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Node_tools/make_structure_map " + " " + final_node_fname[0] + " " + orig_top_fname[0] + " " + final_node_fname[1] + " " + orig_top_fname[1] + " " + "01.map")
os.system("/localhome/py09bh/Software/FFEA/FFEA_git/FFEA_tools/FFEA_initialise/Node_tools/make_structure_map " + " " + final_node_fname[1] + " " + orig_top_fname[1] + " " + final_node_fname[0] + " " + orig_top_fname[0] + " " + "10.map")

# Remove stuff that was made with this script
delete_files()
