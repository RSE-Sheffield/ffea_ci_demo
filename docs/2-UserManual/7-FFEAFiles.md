
Input and Output files {#ioFiles}
======================

Input files {#iFiles}
=====================

Main runner input file: .ffea   {#iffea}
------------------------------
 Corresponding to the main input file for the FFEA runner, this file is widely described
  [here](\ref ffea_iFile), and [here](\ref keywordReference)



Nodes file: .node   {#ifnode}
-----------------

Contains a list of white-space separated coordinates specifying the position
 of each node in the blob. The file starts with a header of four lines, and
 is followed by two separated lists with the interior and exterior nodes:


    ffea node file
    num_nodes 305
    num_surface_nodes 242
    num_interior_nodes 63
    surface nodes:
    x0 y0 z0
    x1 y1 z1
    ...
    x241 y241 z241
    interior nodes:
    x242 y242 z242
    x243 y243 z243
    ...
    x304 y304 z304



Topology file: .top   {#iftop}
-------------------

Contains a list a white-space separated indices specifying the node indices
contained within each element. There are ten indices per element as they are second
order elements (every tetrahedron has a node at each corner
 and at the midpoint of each edge). Separated into internal (no surface faces)
 and external (surface faces) elements for use in surface-surface interactions.


    ffea topology file
    num_elements 130
    num_surface elements 120
    num_interior elements 10
    surface elements:
    n00 n01 n02 n03 n04 n05 n06 n07 n08 n09
    n10 n11 n12 n13 n14 n15 n16 n17 n18 n19
    ...
    n110 n111 n112 n113 n114 n115 n116 n117 n118 n119
    interior elements:
    n120 n121 n122 n123 n124 n125 n126 n127 n128 n129



Surface file: .surf   {#ifsurf}
-------------------

Contains a list of white-space separated indices specifying the indices of the
nodes on each surface face and the index of the containing element.


    ffea surface file
    num_surface_faces 480
    faces:
    el0 n00 n01 n02
    el1 n10 n11 n12
    ...
    el480 n4800 n4801 n4802



Material file: .mat   {#ifmat}
-------------------

Contains a list of white-space separated material parameters for each element.


    ffea material params file
    num elements 120
    density0 shear_visc0 bulk_visc0 shear_modulus0 bulk_modulus0 dielectric_const0
    density1 shear_visc1 bulk_visc1 shear_modulus1 bulk_modulus1 dielectric_const1
    ...
    density119 shear_visc119 bulk_visc119 shear_modulus119 bulk_modulus119 dielectric_const119



Stokes radii file: .stokes   {#ifstokes}
--------------------------

Contains a list of effective hydrodynamic radii for each node in the system
for the calculation of external drag.

    ffea stokes radii file
    num_nodes 305
    radius0
    radius1
    ...
    radius 304


Van der Waals file: .vdw    {#ifvdw}
-------------------------
Contains a list of integers, ranging from -1 to 6, describing the type of vdw interaction
 this face will undergo, as described by the Lennard-Jones files, where:
   * -1 - stands for inactive faces.
   * 0 to 6 - specifies the type of an active face, interacting
              with parameters defined in .lj file.


    ffea vdw file
    num_faces 480
    vdw_params:
    type0
    type1
    ...
    type479


Lennard-Jones file: .lj    {#iflj}
-----------------------

Contains a matrix of parameter pairs for the Lennard-Jones interactions
 between the different face types specified in the .vdw file.
 A description of this file can be found in the [Lennard-Jones](\ref ljPotential) section.

<!-- In the following example:


    ffea vdw forcefield params file
    num vdw face types 7
    (1e+15, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)
    (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09) (1e+12, 1e-09)


face type0 interacts with face type0 through a Lennard-Jones potential
 that has a well depth of \f$\epsilon_{0,0} = 10^{15} J/m^2\f$ and an
 equilibrium distance \f$\sigma_{0,0} = 10^{-9} m\f$,
 whereas face type 0 interacts with face type 1
 (and in fact, all other interaction types) with a Lennard-Jones potential
 with parameters \f$\epsilon_{0,1} = 10^{12} J/m^2\f$ and
  \f$\sigma_{0,1} = 10^{-9} m\f$. Therefore, this matrix should be symmetric.
-->


Pre-computed potentials: .pot and .force  {#potfile}
----------------------------------------
These files contain values for the potential and force in a two column format:


    0 1606.41515966667
    0.002 1595.931897
    0.004 1580.257003
    0.006 1564.662109


where the first column give positions and the second column has energies and forces
 respectively. The runner input file will expect a conversion factor to \f$m\f$ and
 \f$J\f$ (and \f$Jm\f$) for them (see the [pre-computed potentials](\ref fm_inputfile)
 for more details).


Springs file: .springs   {#ifsprings}
----------------------

This file contains a set of Hookean springs joining pairs of nodes to allow the addition of arbitrary interactions, if necessary.
 After an introductory header, every line specifies the details of a single spring through:
 spring constant (in ` N/m `), equilibrium distance (in meters)
 blob i, blob j, conformation k, conformation l, node m, node n.
Specifically, a springs file with 3 springs could look like:


    ffea springs file
    num_springs 3
    springs:
    1e-2 5e-9 0 1 0 0 40 40
    1e-2 5e-9 0 1 0 0 57 57
    1e-2 3e-9 0 1 0 0 66 83


Pinned nodes file: .pin   {#ifpin}
----------------------

This file provides a list of nodes for the corresponding blob that will remain fixed
 in its initial position. The format for a valid ` .pin ` file with three nodes pinned would be:


    ffea pinned nodes file
    num_pinned_nodes 0
    pinned nodes:
    56
    59
    445



Constant forces file: .ctforces {#ifctforces}
---------------------------------------------

This file provides a list of nodes or surfaces where a certain force is to be applied.
 The forces can be linear or circular, and the header
 specifies the total number of constant forces to take from
 the input file, the number of linear forces onto nodes, the number of circular forces
 onto nodes, and the number of the linear forces onto faces. If any of these is zero,
 it can be ommitted.
 Then, the file has a section on linear forces, and a file on circular forces,
 needed in this order but appearing only if the correspondent number of forces
 is higher than zero.

In the linear forces on nodes section,
 each line writes Force module (in ` N `), a unit vector
 with the direction of the force, blob index, conformation index, and node
 number ` or ` string all so that the same force is applied onto every on this
 blob.

In the rotational forces on nodes section,
 each line writes Force module (in ` N `), a two letter keyword
 (see later), a point and an axis, specifying the torque rotation axis,
 and three integer specifying blob, conformation and node where the force
 is being applied. The first letter of the two-letter keyword specifies
 if the point given associated to the axis is given in Cartesian coordinates
 (`p`) or as a triad ` blob `, `conf`, `node` (`n`, **NOT SUPPORTED YET**).
 The second one
 specifies if the circular force is a torque (`t`) \f$\tau = r\times F\f$ or
 a constant force (`f`) \f$\tau = \hat{r}\times F\f$.

In the linear forces on surfaces section,
 a force F is spread evenly on a surface composed of a number of faces. Thus,
 every line specifies the magnitude of the force (in ` N ` ), a unit
 vector with the direction of the force, blob index, conformation index,
 followed by a list of face indices. At every time step the total area
 of the selected faces will be calculated \f$A = \sum_i a_i\f$, and
 every face will receive a force \f$ f_i = F a_i / A\f$ at its centre.


Thus, a valid ` ctforces ` file with force 1e-9 `N` applied onto
 every node on blob 1 towards the z axis, 3 nodes on blob 2 pushed
 towards x would look like, and all the nodes of blob 0 suffering a
 torque resulting of a force with magnitude 1e-10 `N` around the
 axis going through point (-1.9e-9, -1.4313e-9, 8.911e-10) `m`
 and direction (0,1,0) would look like:


    ffea ctforces file
    num_ctforces 5
    num_linear_forces 4
    num_rot_forces 1
    num_linear_surface_forces 1
    linear forces:
    1e-9 0 0 1 1 0 all
    1e-9 1 0 0 2 0 104
    1e-9 1 0 0 2 0 134
    1e-9 1 0 0 2 0 152
    rotational forces:
    1e-10 pt -1.9e-9 -1.4313e-9 8.911e-10 0 1 0 0 0 all  <!-- node 1402 -->
    linear surface forces:
    1e-13 1.000 0.000 0.000 0 0 280 281 282 283 284 285

As in the main input FFEA file, XML comments are supported.




Checkpoint file .fcp {#ffeaCheckpointFileIn}
----------------------------------------------
This is the same type of file than the output [checkpoint file](\ref ffeaCheckpointFileOut).


States file .states {#ffeaStatesFileIn}
-----------

This file contains data regarding the kinetic states that can be activated within this blob.
 Each state is defined through two consecutive lines. Line 1 defines the conformation index
 active in this state. Line 2 defines any binding site couplings:


    ffea kinetic states file
    num_states 2
    states:
    conformation_index 0
    binding
    conformation_index 1
    binding


Rates file .rates {#ffeaRatesFileIn}
----------

This file defines the rates (units - Hz) of kinetic switching between states. It is formatted as a matrix,
 num_states X num_states where rows are the base state, and columns are the target state:


    ffea kinetic rates file
    num_states 2
    rates:
    0 1e10
    1e9 0

These rates are converted into switching probability / timestep within the code,
 so if your rate asks for a higher than 100% probability of
 switching at each timestep, FFEArunner will return an error. The identity rates (leading diagonal) are actually not used
 by the FFEArunner; because the sum of probabilities in each row must equal 1, the code calculates this value using all of the others.

Map file: .map {#ffeaMapFileIn}
--------------

This file defines the linear map (matrix) between nodes of different conformations. Due to the size and sparsity, this matrix
 is written in the Yale Sparse matrix format. Below is the matrix defining the transformation of a cube onto itself (identity matrix):


    FFEA Kinetic Conformation Mapping File (Sparse)
    num_nodes_from 8
    num_nodes_to 8
    num_entries 8
    map:
    entries - 1 1 1 1 1 1 1 1
    key - 0 1 2 3 4 5 6 7 8
    columns - 0 1 2 3 4 5 6 7

Note: This map should not be calculated manually! FFEAtools provides a [tool](\ref kffea_implementation) to help with the
 automation of this process.

Output files {#oFiles}
======================

Trajectory file: .ftj    {#oftrajectory}
---------------------
Contains a list of `*` separated frames specifying the structure of the blob at
each outputted time-step. The frames contain a list of blobs which themselves
contain a list of node positions, velocities and forces. If the motion state is
STATIC, blob does not move, so frame is not outputted to save space, and improve
the performance. Only the positions are used in the viewer.
 In between each frame is a small segment showing whether any conformational
 changes have occurred. As such,
each frame can be a different length, which is why the blob and conformation
sizes are specified at the beginning of the file.


    FFEA trajectory file

    Initialisation:
    Number of Blobs 2
    Number of Conformations 2 1
    Blob 0: Conformation 0 Nodes 305	Conformation 1 409
    Blob 1: Conformation 0 Nodes 305

    *
    Blob 0, Conformation 0, step 0
    DYNAMIC
    x0 y0 z0 velx0 vely0 velz0 phi0 forcex0 forcey0 forcez0
    x1 y1 z1 velx1 vely1 velz1 phi1 forcex1 forcey1 forcez1
    ...
    Blob 1, Conformation 0, step 0 DYNAMIC x0 y0 z0 velx0 vely0 velz0 phi0
    forcex0 forcey0 forcez0
    x1 y1 z1 velx1 vely1 velz1 phi1 forcex1 forcey1 forcez1
    ...
    Conformation Changes:
    Blob 0: Conformation 0 -> Conformation 1
    Blob 1: Conformation 0 -> Conformation 0
    *
    Blob 0, Conformation 1, step 1001
    ...
    *



Measurement files: .fm / .fdm   {#ofmeasurement}
------------------------

Contains a list of the relevant system properties of each outputted time-step.
Two separate types of file for internal and external measurements.


### Global measurements: .fm ###

Every simulation will print this file, which contains measurements of the entire
system. Internal properties such as strain and kinetic energies, Centroids, Angular momentum, RMSD etc
are included, as well as the total interaction energies within the system from VdW, springs and Precomp potentials.
Additionally, system details and a copy of the total input parameter set are written to this file.
 The format is:


	FFEA Global Measurement File

	Simulation Details:
		Simulation Began on DD/MM/YYYY at hh:mm:ss
		Script Filename = /path/to/script/script.ffea
		Simulation Type = Full

	Parameters:
		restart = 0
		dt = 1.000000e-14
		...

	Measurements:
	Time          StrainEnergy  Centroid.x    Centroid.y    Centroid.z    RMSD
	0.000000e+00  0.000000e+00  8.095944e-09  1.491714e-09  5.000000e-10  0.000000e+00
	5.000000e-11  2.478548e-19  8.017644e-09  1.455046e-09  5.094042e-10  1.446259e-10
	1.000000e-10  2.868027e-19  7.978989e-09  1.376037e-09  4.846020e-10  2.252258e-10
	...


A small thing to note; if mass is not included, momenta and kinetic energies are not written. If VdW is not active, it is not written. Anything optional may not be written,
as seen here with KineticEnergy, SpringEnergy, SurfSurfEnergy and PreCompEnergy.
This is to avoid writing out many zeroes, gain some performance and save disk space.
The [FFEA toolkit](\ref FFEAanalysistut) is equipped to read in these files in general.

### Detailed measurements: .fdm ###

This file records the measurements on individual blobs, and between each specific pair of blobs. It is created by default but can be suppressed with '-d' and the command line.
and allows one to know specific details about individual molecules / pairs of molecules in a multi-blob simulation. Measurements specific to each blob are printed first, followed by each interacting pair.
 The format is:


	FFEA Detailed Measurement File

	Measurements:
	Time          | B0 StrainEnergy  Centroid[0]    Centroid[1]    Centroid[2]    RMSD          | B1 StrainEnergy  Centroid[1]    Centroid[1]    Centroid[2]    RMSD          | B0B1 SpringEnergy  SurfSurfEnergy
	0.000000e+00       0.000000e+00  8.095944e-09  1.491714e-09  5.000000e-10  0.000000e+00       0.000000e+00  8.095944e-09  1.491714e-09  5.000000e-10  0.000000e+00	   0.000000e+00  0.000000e+00
	5.000000e-11       2.478548e-19  8.017644e-09  1.455046e-09  5.094042e-10  1.446259e-10       2.478548e-19  8.017644e-09  1.455046e-09  5.094042e-10  1.446259e-10 	   5.123413e-21  0.000000e+00
	1.000000e-10       2.868027e-19  7.978989e-09  1.376037e-09  4.846020e-10  2.252258e-10       2.868027e-19  7.978989e-09  1.376037e-09  4.846020e-10  2.252258e-10	   3.567162e-21  0.000000e+00
	...


Again, if something is not include, it's associated column is not written. So here, Blob 0 had no mass but Blob 1 did. Additionally, 2 interaction types were included between blobs 0 and 1 but VdW is currently out of range,
hence the zeroes. No blobs were interacting with themselves.

Checkpoint file .fcp {#ffeaCheckpointFileOut}
------------------------------------------------
The checkpoint file stores the state of the Random Number Generator(s) RNG(s) at the
 last saved step. The format of this file, that should be automatically written,
 starts with a single header line specifying the number of RNGStreams dedicated to
 the thermal stress:

     RNGStreams dedicated to the thermal stress: <N>

then it follows ` <N> ` lines with 6 integers per line specifying the state of the
  stream, e. g.:

     1828815773 2063911640 1007531883 4278127031 1691073935 2862951708


Afterwards, in the case of using kinetics, the following line will be found:

    RNGStream dedicated to the thermal stress:\n

followed by a single line with 6 integers describing the state of
   the RNGStream dedicated to the kinetic stress.

