import sys, os
import numpy as np
import copy
from matplotlib import pyplot as plt

if len(sys.argv) < 6:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .evals file 1] [INPUT .evals file 2] [INPUT .evecs file 1] [INPUT .evecs file 2] [OUTPUT fname] OPTIONAL{[EIG SWITCH (x,y)]}")

inevalsfile = [sys.argv[1], sys.argv[2]]
infile = [sys.argv[3], sys.argv[4]]
out_fname = sys.argv[5]
evec_fname = [[],[]]
evecs = [[],[]]
evals = [[],[]]

eig_switch = []
if len(sys.argv) > 6:
	for arg in sys.argv[4:]:
		try:
			eig_switch.append([int(a) for a in arg.split(",")])
		except:
			sys.exit("EIG SWITCH input error. Pairs of integer values separated by columns please")

# Get motion eigenvalues (and convert to spatial ones for comparison)
print("Reading eigenvalues from inputs...")
for i in range(2):
	fin = open(inevalsfile[i])
	for line in fin.readlines():
		evals[i].append(float(line))
	fin.close()

	# This means spring constants (hopefully...)
	if evals[i][0] < 1.0:
		conversion = True
	else:
		conversion = False
	if conversion:
		for j in range(len(evals[i])):
			evals[i][j] = 0.411 / evals[i][j]	# Converted to angstroms	


print("done!")

# Get motion eigenvectors
print("Reading eigenvectors from inputs...")
for i in range(2):
	fin = open(infile[i])
	for line in fin.readlines():
		sline = line.split()
		evecs[i].append(np.array([float(s) for s in sline]))
	fin.close()
	#evecs[i] = np.array(evecs[i])

print("done!")

# In one of these, switch the eigenvectors of necessary
if eig_switch != []:
	for e in eig_switch:
		try:
			evec1 = copy.copy(evecs[0][e[0]])
			evec2 = copy.copy(evecs[0][e[1]])
			eval1 = copy.copy(evals[0][e[0]])
			eval2 = copy.copy(evals[0][e[1]])
		except:
			sys.exit("EIG SWITCH values greater than number of eigenvectors present. Please try again")

		evecs[0][e[0]] = evec2
		evecs[0][e[1]] = evec1
		evals[0][e[0]] = eval2
		evals[0][e[1]] = eval1

# Do some dot products!
print("Building the Eigenvectors Dot Product Matrix...")
num_modes = len(evecs[0])
dot_prod = np.array([[0.0 for i in range(num_modes)] for j in range(num_modes)])

for i in range(num_modes):
	for j in range(num_modes):
		dot_prod[i][j] = np.fabs(np.dot(evecs[0][i], evecs[1][j]))

print("done!")

# Print to file
print("Writing Matrix to file " + os.path.basename(out_fname) + "...")
outfile = open(out_fname, "w")
outfile.write("Compare pyPca modes\n\nEigen Set 1 x Eigen Set 2\n\n")
for i in range(num_modes):
	for j in range(num_modes):
		outfile.write("%6.3f " % (dot_prod[i][j]))

	outfile.write("\n")

outfile.close()
print("done!")

# Remove crap
for i in range(2):
	for fname in evec_fname[i]:
		os.system("rm " + fname)

#
# Plot the eigenvector heatmap thing
#

print("Visualising eigenvectors as a heatmap...")
#os.system(os.path.dirname("python " + os.path.abspath(sys.argv[0])) + "/plot_eigensystem_comparison.py " + out_fname)

# Get figure properties
fig = plt.figure(figsize=(20,10))
fig.suptitle("Eigensystem Comparison", fontsize=24)
ax = fig.add_subplot(1,2,1)
#fig, ax = plt.subplots(figsize=(13,10))

#column_labels = list('43210')
#row_labels = list('01234')
column_labels = []
row_labels = []
#xticks = []	# cols
#yticks = []	# rows

# We only want 11 labels if num_modes > 20!

if num_modes <= 20:
	for i in range(num_modes):
		row_labels.append(str(i))
		#column_labels.append(str((num_modes - 1) - i))
		column_labels.append(str(i))
else:
	count = 0
	for i in range(num_modes + 1):
		
		if i % int(num_modes / 10.0) == 0:
			row_labels.append(int((count / 10.0) * num_modes))
			#column_labels.append(int((1 - (count / 10.0)) * num_modes))
			column_labels.append(int((count / 10.0) * num_modes))
			count += 1
		else:
			pass


# Plot data
data = range(dot_prod.shape[0])
for i in range(num_modes):
	data[i] = dot_prod[(num_modes - 1) - i]

data = np.array(data)
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
heatmap.set_clim([0, 1])

# And 'plot' colorbar
cbar = plt.colorbar(heatmap, ticks=[i / 10.0 for i in range(11)])
	
# All tick
if num_modes <= 20:
	ax.set_xticks(np.arange(data.shape[0]) + 0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[1], 0, -1) - 0.5, minor=False)
else:
	ax.set_xticks(np.array([int((i / 10.0) * num_modes) for i in range(11)]), minor=False)
	ax.set_yticks(np.array([int((1 - (i / 10.0)) * num_modes) for i in range(11)]), minor=False)

# Labels
ax.set_yticklabels(row_labels, fontsize=18)
ax.set_xticklabels(column_labels, fontsize=18)

# And now titles and stuff
ax.set_xlim(0, num_modes)
ax.set_ylim(0, num_modes)

ax.set_title("Eigenvector Dot Product Array", fontsize=24)
ax.set_xlabel("Eigensystem A Modes", fontsize=18)
ax.set_ylabel("Eigensystem B Modes", fontsize=18)

cbar.ax.set_yticklabels(["%2.1f" % (i / 10.0) for i in range(11)], fontsize=18)
cbar.set_label('Normalised Dot Product Values', fontsize=18)

#
# Now eigenvalues
#
print("Visualising eigenvalues as a histogram...")
ax = fig.add_subplot(1,2,2)

evals[0] = np.array(evals[0])
evals[1] = np.array(evals[1])

ind = np.arange(evals[0].size)  # the x locations for the groups
width = 0.35
eva = ax.bar(ind, evals[0], width, color='r')
evb = ax.bar(ind + width, evals[1], width, color='g')
ax.set_title("Eigenvalue Bar Chart Comparison", fontsize=24)
ax.set_ylabel("Eigenvalue " + r"$(\AA ^2)$", fontsize=18)
ax.set_xlabel("Index", fontsize=18)
ax.legend([eva,evb], ['Eigensystem A', 'Eigensystem B'], loc = 1, fontsize=12)

base, ext = os.path.splitext(out_fname)
plt.savefig(base + ".png")
plt.show()
print("done!")
