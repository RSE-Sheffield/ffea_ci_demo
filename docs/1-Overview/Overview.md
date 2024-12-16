# Overview {#overview}

Fluctuating Finite Element Analysis is a new molecular modelling algorithm, designed to support systems that are larger and more complex than those modelled by atomistic molecular dynamics. Instead of modelling biological systems as a collection of connected atoms, it models them as 3D volumes comprised of tetrahedrons. Unlike previous coarse-grained models, the models FFEA generates are visco-elastic continuum solids. Unlike other applications of Finite Element Analysis, these systems are subject to thermal fluctuations.

This technique has the potential to model large, complex systems, made of many molecules, and complex processes at the frontiers of molecular biology. As it does not not require an atomistic level of detail, it can also be used to simulate biological molecules that cannot be imaged using X-ray crystallography.

# Features  {#features}

 * Protein-protein interactions:
     * A 6-12 Lennard-Jones potential.
     * A repulsive potential that is proportional to the overlapping volume.
     * Specific interactions defined using [precomputed potentials](\ref fmApproach).
 * Kinetic state changes can be simulated together with the continuum model to account for conformational changes and binding events.
 * Conversion tools for EM density data and atomistic structures into FFEA simulations.
 * A plugin for PyMOL, allowing the visualisation of FFEA systems and trajectories.
 * Initialisation and analysis tools available on the command line and under a Python API.
 * [KOBRA model](\ref rods) for slender biological objects such as coiled-coils.
 * Lees-Edwards boundary conditions.
 * An extremely high degree of reproducity - check out our integration tests!
 * Distance restraints with harmonic potentials, such as springs.

# Getting Started  {#gettingstarted}

FFEA is free to download and use under the GPLv3 software license. We provide binary releases, and building from source is relatively painless.

* We recommend you **download the latest x86_64 [binary release](https://bitbucket.org/FFEA/ffea/downloads)** from BitBucket, otherwise, [compile from source](\ref configure) by cloning the repository. For the latest bleeding-edge features, switch to the [development branch](https://bitbucket.org/FFEA/ffea/src/superdev/) instead.
* If you're compiling from source, install FFEA according to the instructions found in the [installation guide](\ref install), which includes the complete list of [software requirements](\ref prerequisites).
* Once FFEA is installed, consult the [first-time user tutorial](\ref Tutorial). For KOBRA rods, try the [rods tutorial](\ref Tutorial) instead.

# Publications  {#publications}

## Methodology:
* Welch R. J., Harris S. A., Harlen O. G. & Read D. J. ["KOBRA: A Fluctuating Elastic Rod Model for Slender Biological Macromolecules"](https://doi.org/10.1039/D0SM00491J) (2020), Soft Matter 16, 7544-7555.
* Solernou A., Hanson B. S., Richardson R. A., Welch R., Harris S. A., Read D. J., Harlen O. G. ["Fluctuating Finite Element Analysis (FFEA): A continuum mechanics software tool for mesoscale simulation of biomolecules"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005897) (2018), PLoS Comput. Biol. 14(3): e1005897.
* Oliver R., Read D. J., Harlen O. G. & Harris S. A. ["A Stochastic finite element model for the dynamics of globular macromolecules"](http://www.sciencedirect.com/science/article/pii/S0021999112007589) (2013), J. Comp. Phys. 239, 147-165.
* Patargias G. N., Harris S. A. & Harding J. ["A demonstration of the inhomogeneity of the local dielectric response of proteins by molecular dynamics simulations"](https://www.ncbi.nlm.nih.gov/pubmed/20572740) (2010), J. Chem. Phys. 132, 235103.

## Applications:
* Richardson R. A., Hanson B. S., Read D. J., Harlen O. G. & Harris S. A. ["Exploring the dynamics of flagellar dynein within the axoneme with Fluctuating Finite Element Analysis"](https://doi.org/10.1017/S0033583520000062) (2020), Q. Rev. Biophys. 53, E9.
* Lee S. C., Collins R., Lin Y., Jamshad M., Broughton C., Harris S. A., Hanson B. S., Tognoloni C., Parslow R. A., Terry A. E., Rodger A., Smith C. J., Edler K. J., Ford R., Roper D. I. & Dafforn T. R. ["Nano-encapsulated Escherichia coli Divisome Anchor ZipA, and in Complex with FtsZ"](https://doi.org/10.1038/s41598-019-54999-x) (2019), Sci. Rep. 9, 18712.
* Richardson R., Papachristos K., Read D. J., Harlen O. G., Harrison M. A., Paci E., Muench S. P. & Harris S. A ["Understanding the apparent stator-rotor connections in the rotary ATPase family using coarse-grained computer modelling"](https://www.ncbi.nlm.nih.gov/pubmed/25174610) (2014), Proteins: Struct., Funct., Bioinf. 82, 3298-3311.

## Reviews:
* Gray A., Harlen O. G., Harris S. A., Khalid S., Leung Y. M., Lonsdale R., Mulholland A. J., Pearson A. R., Read D. J. & Richardson R. A. ["In pursuit of an accurate spatial and temporal model of biomolecules at the atomistic level: a perspective on computer simulation"](https://www.ncbi.nlm.nih.gov/pubmed/25615870) (2015), Acta Cryst. D71, 162-172.
* Hanson B., Richardson R., Oliver R., Read D. J., Harlen O. & Harris S. ["Modelling biomacromolecular assemblies with continuum mechanics"](https://www.ncbi.nlm.nih.gov/pubmed/25849915) (2015), Biochem. Soc. Trans. 43, 186-192.
* Oliver R., Richardson R. A., Hanson B., Kendrick K., Read D. J., Harlen O. G. & Harris S. A. ["Modelling the Dynamic Architecture of Biomaterials Using Continuum Mechanics"](http://link.springer.com/chapter/10.1007%2F978-3-319-09976-7_8) (2014), Protein Modelling, G. Náray-Szabó, Editor. Springer International Publishing. p. 175-197.
       

# Contribute

Do you have a research question that FFEA could help to answer?

   * Try FFEA and let us know how you're using the software.
   * Send bug reports, questions and feature requests to our [issue tracker](https://bitbucket.org/FFEA/ffea/issues?status=new&status=open)
   * [Fork us](https://bitbucket.org/FFEA/ffea/fork) and help develop FFEA!

If you have questions and comments, please contact us! For biophysics, contact Sarah Harris ([S.A.Harris@leeds.ac.uk](mailto:S.A.Harris@leeds.ac.uk)). For finite elements, contact Oliver Harlen ([O.G.Harlen@leeds.ac.uk](mailto:O.G.Harlen@leeds.ac.uk)). For KOBRA and stochastic enquiries, contact Daniel Read ([D.J.Read@leeds.ac.uk](mailto:D.J.Read@leeds.ac.uk)). For software engineering, contact Joanna Leng ([J.Leng@leeds.ac.uk](mailto:J.Leng@leeds.ac.uk)).


# FFEA Team  {#FFEAteam}

   * Ryan Cocking
   * Molly Gravett
   * Ben Hanson
   * [Oliver Harlen](https://eps.leeds.ac.uk/maths/staff/4035/dr-oliver-harlen)
   * [Sarah Harris](https://eps.leeds.ac.uk/physics/staff/4109/dr-sarah-harris)
   * [Joanna Leng](https://eps.leeds.ac.uk/computing/staff/1509/joanna-leng)
   * Robin Oliver
   * [Daniel Read](https://eps.leeds.ac.uk/maths/staff/4074/professor-daniel-read)
   * Robin Richardson
   * Tom Ridley
   * Jarvellis Rogers
   * Albert Solernou
   * [Rob Welch](http://robwel.ch/)
   * Changqing Zheng
   

# Thanks {#thanks}

We want to thank everybody who has helped in making FFEA possible, from 
 summer students, to experimental professors, we would not be here without you:
 
   * Neelofer Banglawala
   * Jana Boltersdorf
   * Jonathan Boyle
   * Stan Burgess
   * Glenn Carrington
   * Samantha Coffey
   * Toni Collis
   * Mike Croucher
   * Ash Dwarka
   * Matthew Faulkner 
   * Ashley Fenton
   * Katrina Goldman
   * Thijs van der Heijden
   * Guanhao Lu
   * Stephen Muench
   * Michelle Peckham
   * Paul van der Schoot
   * Kerrie Smith 
   * Kees Storm
   * Ondřej Vysocký
   * Christopher Woods
   
   