import os, sys, math
from Vectors import vector3

class Force:
	def __init__(self):
		self.node = 0
		self.force_type = "no_force"
		self.constant_force = vector3(0.0, 0.0, 0.0)
		self.constant_force_magnitude = 0.0
		self.constant_torque = vector3(0.0, 0.0, 0.0)
		self.point_in_torque_axis = vector3(0.0, 0.0, 0.0)
		self.constant_torque_magnitude = 0.0
		self.variable_torque_base = 0
		self.variable_torque_tip = 0
		self.variable_torque_magnitude = 0.0
		
def add_constant_force(nodes, forces, forced_nodes):
	print "By volume, which nodes will we apply a constant force to:\n"
	x_min = input("Enter x_min:")
	x_max = input("Enter x_max:")
	y_min = input("Enter y_min:")
	y_max = input("Enter y_max:")
	z_min = input("Enter z_min:")
	z_max = input("Enter z_max:")

	print "Calculating nodes to force..."
	nodes_to_force = []
	for i in range(0, num_nodes):
		if nodes[i].x > x_min and nodes[i].x < x_max and nodes[i].y > y_min and nodes[i].y < y_max  and nodes[i].z > z_min and nodes[i].z < z_max:
			nodes_to_force.append(i) 
			forced_nodes.append(i)
			
			
	print "We are forcing " + str(len(nodes_to_force)) + " nodes in this range."
	if len(nodes_to_force) == 0:
		return
	
	print "Input Force Vector (doesn't have to be unit) and a magnitude:"
	temp = vector3(0.0, 0.0, 0.0)	
	temp.x = input("Enter x component:")
	temp.y = input("Enter y component:")
	temp.z = input("Enter z component:")
	magnitude = input("Enter magnitude:")
	
	for i in range(len(nodes_to_force)):
		force = Force()
		force.node = nodes_to_force[i]
		force.force_type = "constant"
		force.constant_force = temp
		force.constant_force_magnitude = magnitude
		forces.append(force)

def add_constant_torque(nodes, forces, forced_nodes):
	print "By volume, which nodes will we apply a constant torque to:\n"
	x_min = input("Enter x_min:")
	x_max = input("Enter x_max:")
	y_min = input("Enter y_min:")
	y_max = input("Enter y_max:")
	z_min = input("Enter z_min:")
	z_max = input("Enter z_max:")

	print "Calculating nodes to force..."
	nodes_to_force = []
	for i in range(0, num_nodes):
		if nodes[i].x > x_min and nodes[i].x < x_max and nodes[i].y > y_min and nodes[i].y < y_max  and nodes[i].z > z_min and nodes[i].z < z_max:
			nodes_to_force.append(i) 
			forced_nodes.append(i)
			
	print "We are forcing " + str(len(nodes_to_force)) + " nodes in this range."
	if len(nodes_to_force) == 0:
		return

	print "Defining torque:"
	print "Input vector defining a point in torque axis, vector defining direction of torque (doesn't have to be unit) and a magnitude:"
	temp_point = vector3(0.0, 0.0, 0.0)
	temp_torque = vector3(0.0, 0.0, 0.0)
	temp_point.x = input("Enter point in axis x component:")
	temp_point.y = input("Enter point in axis y component:")
	temp_point.z = input("Enter point in axis z component:")
	temp_torque.x = input("Enter torque direction x component:")
	temp_torque.y = input("Enter torque direction y component:")
	temp_torque.z = input("Enter torque direction z component:")
	magnitude = input("Enter torque magnitude:")

	for i in range(len(nodes_to_force)):
		force = Force()
		force.node = nodes_to_force[i]
		force.force_type = "constant_torque"
		force.point_in_torque_axis = temp_point
		force.constant_torque = temp_torque
		force.constant_torque_magnitude = magnitude
		forces.append(force)

def add_variable_torque(nodes, forces, forced_nodes):
	print "By volume, which nodes will we apply a constant torque to:\n"
	x_min = input("Enter x_min:")
	x_max = input("Enter x_max:")
	y_min = input("Enter y_min:")
	y_max = input("Enter y_max:")
	z_min = input("Enter z_min:")
	z_max = input("Enter z_max:")

	print "Calculating nodes to force..."
	nodes_to_force = []
	for i in range(0, num_nodes):
		if nodes[i].x > x_min and nodes[i].x < x_max and nodes[i].y > y_min and nodes[i].y < y_max  and nodes[i].z > z_min and nodes[i].z < z_max:
			nodes_to_force.append(i) 
			forced_nodes.append(i)
			
	print "We are forcing " + str(len(nodes_to_force)) + " nodes in this range."
	if len(nodes_to_force) == 0:
		return

	print "Defining torque:"
	print "Input two nodes, which will define the torque axis, in the direction you want the axis to point, followed by a magnitude:"
	temp_node1 = input("Enter first node:")
	temp_node2 = input("Enter second node:")
	magnitude = input("Enter torque magnitude:")

	for i in range(len(nodes_to_force)):
		force = Force()
		force.node = nodes_to_force[i]
		force.force_type = "variable_torque"
		force.variable_torque_base = temp_node1
		force.variable_torque_tip = temp_node2
		force.variable_torque_magnitude = magnitude
		forces.append(force)

def write_output(force_fname, forces, forced_nodes):
	force_out = open(force_fname, "w")
	force_out.write("ffea additional forces file\n")
	force_out.write("num_forced_nodes " + str(len(list(set(forced_nodes)))) + "\n")
	force_out.write("node|force type|stuff for type:\n")
	for force in forces:
		if force.force_type == "constant":
			force_out.write(str(force.node) + " " + force.force_type + " " + str(force.constant_force.x) + " " + str(force.constant_force.y) + " " + str(force.constant_force.z) + " " + str(force.constant_force_magnitude) + "\n")
		elif force.force_type == "constant_torque":
			force_out.write(str(force.node) + " " + force.force_type + " " + str(force.point_in_torque_axis.x) + " " + str(force.point_in_torque_axis.y) + " " + str(force.point_in_torque_axis.z) + " " + str(force.constant_torque.x) + " " + str(force.constant_torque.y) + " " + str(force.constant_torque.z) + " " + str(force.constant_torque_magnitude) + "\n")
		elif force.force_type == "variable_torque":
			force_out.write(str(force.node) + " " + force.force_type + " " + str(force.variable_torque_base) + " " + str(force.variable_torque_tip) + " " + force.variable_torque_magnitude + "\n")

		
if len(sys.argv) != 3:
	sys.exit("Usage: python make_force_file.py [INPUT .NODE FILE] [OUTPUT .FORCE FILE]")

inputnode = open(sys.argv[1], "r")
inputnode.readline()
num_nodes = int(inputnode.readline().split()[1])
num_surface_nodes = int(inputnode.readline().split()[1])
num_interior_nodes = int(inputnode.readline().split()[1])
inputnode.readline()
nodes = []

for i in range(0,num_surface_nodes):
	sline = inputnode.readline().split()
	nodes.append(vector3(float(sline[0]), float(sline[1]), float(sline[2])))

inputnode.readline()
for i in range(0,num_interior_nodes):
	sline = inputnode.readline().split()
	nodes.append(vector3(float(sline[0]), float(sline[1]), float(sline[2])))

inputnode.close()

done = 0
forces = []
forced_nodes = []
while(done == 0):
	force_to_add = raw_input("Which force would you like to add, Constant (C), Constant Torque (CT), Variable Torque (VT)? Or shall we finish (f)?:")
	if(force_to_add == "C" or force_to_add == "c"):
		add_constant_force(nodes, forces, forced_nodes)
	elif(force_to_add == "CT" or force_to_add == "ct"):
		add_constant_torque(nodes, forces, forced_nodes)
	elif(force_to_add == "VT" or force_to_add == "vt"):
		add_constant_torque(nodes, forces, forced_nodes)
	elif(force_to_add == "f"):
		write_output(sys.argv[2], forces, forced_nodes)
		done = 1
	else:
		print "Please enter either 'C', 'CT', 'VT' or 'f'\n"

