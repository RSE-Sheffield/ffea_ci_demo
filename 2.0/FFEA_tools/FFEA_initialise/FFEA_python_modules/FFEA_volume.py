import sys
import numpy as np
import FFEA_tetrahedron

class FFEA_volume:

	def __init__(self, vol_fname):
		
		# Open file and check type
		fin = open(vol_fname, "r")
		if fin.readline().strip() == "mesh3d":
			fin.close()
			self.load_from_vol(vol_fname)
		else:
			sys.exit("Error. This is not a volume file recognized by FFEA.\n")
	
	def load_from_vol(self, fname):
	
		fin = open(fname, "r")
		fin.readline()
		
		# Get tetrahedra into array
		while fin.readline().strip() != "volumeelements":
			continue
		
		self.num_elements = int(fin.readline().strip())
		self.element = []
		for i in range(self.num_elements):
			sline = fin.readline().split()
			self.element.append(FFEA_tetrahedron.FFEA_tetrahedron(int(sline[2].strip()) - 1, int(sline[3].strip()) - 1, int(sline[4].strip()) - 1, int(sline[5].strip()) - 1))
		
		# Get nodes into a numpy array
		while fin.readline().strip() != "points":
			continue
		
		self.num_nodes = int(fin.readline().strip())
		temp_nodes = []
		for i in range(self.num_nodes):
			sline = fin.readline().split()
			temp_nodes.append(np.array([float(sline[0].strip()), float(sline[1].strip()), float(sline[2].strip())]))
		
		self.node = np.array(temp_nodes)
		del temp_nodes
		fin.close()

	def calc_volume(self):
	
		self.volume = 0.0
		for e in self.element:
			self.volume += e.calc_volume(self.node)
			
		return self.volume
