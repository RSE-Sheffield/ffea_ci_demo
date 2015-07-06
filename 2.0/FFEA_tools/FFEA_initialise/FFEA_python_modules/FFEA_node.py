import sys, os
import numpy as np

class FFEA_node:

	def __init__(self, node_fname):
		
		# Open file and check initial stuff
		try:
			fin = open(node_fname, "r")
		except(IOError):
			self.num_nodes = 0
			self.num_interior_nodes = 0
			self.num_surface_nodes = 0
			self.pos = []
			return
			
		if(fin.readline().strip() != "ffea node file"):
			sys.exit("Error. Expected but didn't find 'ffea node fname'\n")
		
		self.num_nodes = int(fin.readline().split()[1].strip())
		self.num_surface_nodes = int(fin.readline().split()[1].strip())
		self.num_interior_nodes = int(fin.readline().split()[1].strip())
		self.pos = np.array([[0.0, 0.0, 0.0] for i in range(self.num_nodes)])
		self.centroid = np.array([0.0,0.0,0.0])

		# Surface nodes
		fin.readline()
		for i in range(self.num_surface_nodes):
			sline = fin.readline().split()
			self.pos[i][0] = float(sline[0])
			self.pos[i][1] = float(sline[1])
			self.pos[i][2] = float(sline[2])
			self.centroid += self.pos[i]

		# Interior nodes
                fin.readline()
                for i in range(self.num_surface_nodes, self.num_nodes, 1):
                        sline = fin.readline().split()
                        self.pos[i][0] = float(sline[0])
                        self.pos[i][1] = float(sline[1])
                        self.pos[i][2] = float(sline[2])
			self.centroid += self.pos[i]

		self.centroid *= 1.0/self.num_nodes
		fin.close()

	def write_nodes_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea node file\n")
		fout.write("num_nodes %d\nnum_surface_nodes %d\nnum_interior_nodes %d\n" % (self.num_nodes, self.num_surface_nodes, self.num_interior_nodes))
		
		fout.write("surface nodes:\n")
		for i in range(self.num_nodes):
			if i == self.num_surface_nodes:
				fout.write("interior nodes:\n")
			fout.write("%6.3e %6.3e %6.3e\n" % (self.pos[i][0], self.pos[i][1], self.pos[i][2]))

	def set_centroid(self, centroid):

		# Get translation vector
		new_centroid = np.array(centroid)
		trans = new_centroid - self.centroid
		for i in range(self.num_nodes):
			self.pos[i] += trans
