import sys, os
import FFEA_script

try:
	from matplotlib import pyplot as plt
	import numpy as np
except(ImportError):
	sys.exit("ImportError recieved. Without the numpy and matplotlib libraries, this script cannot run. Please install these, or manually use gnuplot :)")

if len(sys.argv) != 2:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .ffea script]")

# Get args and build objects
script = FFEA_script.FFEA_script(sys.argv[1])

if not os.path.exists(script.params.measurement_out_fname):
	sys.exit("Error. Measurement file not found. Please supply a script with a completed measurement file.")

meas = script.load_measurement()
top = [script.load_topology(i) for i in range(script.params.num_blobs)]

# We need to plot a global measurement graph, and a graph for every blob
total_num_nodes = 0
total_num_mass_nodes = 0

for i in range(script.params.num_blobs):
	total_num_nodes += len(top[i].get_linear_nodes())
	if script.blob[i].solver != "CG_nomass":
		total_num_mass_nodes += len(top[i].get_linear_nodes())

kT = script.params.kT

#
# First, global graph (strain and kinetic)
#

plt.figure(0)
cmeas = meas.global_meas
num_steps = len(cmeas["Time"])

# Get x axis data (time)
x = cmeas["Time"] * 1e9	# ns

# And y axis data (all energies)
ys = []
yk = []
for i in range(num_steps):
	ys.append(np.mean(cmeas["StrainEnergy"][0:i + 1]) / kT)
	if cmeas["KineticEnergy"] != None:
		yk.append(np.mean(cmeas["KineticEnergy"][0:i + 1]) / kT)

print "\nGlobal System:\n"

ysEXP = [(3 * total_num_nodes - 6) / 2.0 for i in range(num_steps)]
ysh, = plt.plot(x, ys, label='ysh')
ysEXPh, = plt.plot(x, ysEXP, "-", label='ysEXPh')

yserr = (np.fabs(ys[-1] - ysEXP[-1]) / ysEXP[-1]) * 100.0
print "\tTheoretical strain energy = %f; Simulation strain energy = %f; Error is %f%%" % (ysEXP[-1], ys[-1], yserr) 
if cmeas["KineticEnergy"] != None:
	ykEXP = [(3 * total_num_mass_nodes) / 2.0 for i in range(num_steps)]
	ykh, = plt.plot(x, yk, label='ykh')
	ykEXPh, = plt.plot(x, ykEXP, "-", label='ykEXPh')

	ykerr = (np.fabs(yk[-1] - ykEXP[-1]) / ykEXP[-1]) * 100.0
	print "\tTheoretical kinetic energy = %f; Simulation kinetic energy = %f; Error is %f%%" % (ykEXP[-1], yk[-1], ykerr)

plt.xlabel("Time (ns)")
plt.ylabel("Energy (kT)")
plt.title("Global Energy - Running Average")

# Put details on the graph
if cmeas["KineticEnergy"] != None:
	plt.legend([ysh, ysEXPh, ykh, ykEXPh], ['Global Strain Energy - Sim', 'Global Strain Energy - Theory', 'Global Kinetic Energy - Sim', 'Global Kinetic Energy - Theory'], loc = 4)
else:
	plt.legend([ysh, ysEXPh], ['Global Strain Energy - Sim', 'Global Strain Energy - Theory'], loc = 4)
plt.show()

#
# Now, individual blobs (strain and kinetic), if we can
#
if meas.blob_meas == []:
	sys.exit()

fig = plt.figure(1)
for i in range(script.params.num_blobs):

	ax = fig.add_subplot(np.ceil(np.sqrt(script.params.num_blobs)), np.ceil(np.sqrt(script.params.num_blobs)), i + 1)

	cmeas = meas.blob_meas[i]
	num_nodes = len(top[i].get_linear_nodes())

	# Already got the x axis

	# And y axis data (all energies)
	ys = []
	yk = []
	for j in range(num_steps):
		ys.append(np.mean(cmeas["StrainEnergy"][0:j + 1]) / kT)
		if cmeas["KineticEnergy"] != None:
			yk.append(np.mean(cmeas["KineticEnergy"][0:j + 1]) / kT)

	print "\nBlob %d:\n" % (i)
	ysEXP = [(3 * num_nodes - 6) / 2.0 for j in range(num_steps)]
	ysh, = ax.plot(x, ys, label='ysh')
	ysEXPh, = ax.plot(x, ysEXP, "-", label='ysEXPh')

	yserr = (np.fabs(ys[-1] - ysEXP[-1]) / ysEXP[-1]) * 100.0
	print "\tTheoretical strain energy = %f; Simulation strain energy = %f; Error is %f%%" % (ysEXP[-1], ys[-1], yserr) 
	if cmeas["KineticEnergy"] != None:
		ykEXP = [(3 * num_nodes) / 2.0 for j in range(num_steps)]
		ykh, = ax.plot(x, yk, label='ykh')
		ykEXPh, = ax.plot(x, ykEXP, "-", label='ykEXPh')

		ykerr = (np.fabs(yk[-1] - ykEXP[-1]) / ykEXP[-1]) * 100.0
		print "\tTheoretical kinetic energy = %f; Simulation kinetic energy = %f; Error is %f%%" % (ykEXP[-1], yk[-1], ykerr)

	ax.set_xlabel("Time (ns)")
	ax.set_ylabel("Energy (kT)")
	ax.set_title("Blob %d Energy - Running Average" % (i))

	# Put details on the graph
	if cmeas["KineticEnergy"] != None:
		ax.legend([ysh, ysEXPh, ykh, ykEXPh], ['Blob %d Strain Energy - Sim' % (i), 'Blob %d Strain Energy - Theory' % (i), 'Blob %d Kinetic Energy - Sim' % (i), 'Blob %d Kinetic Energy - Theory' % (i)], loc = 4)
	else:
		ax.legend([ysh, ysEXPh], ['Blob %d Strain Energy - Sim' % (i), 'Blob %d Strain Energy - Theory' % (i)], loc = 4)
plt.show()
