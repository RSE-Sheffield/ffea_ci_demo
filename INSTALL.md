# Installation {#install}

This document gives instructions on how to build and install the FFEA package,
 consisting of the FFEA runner and FFEA tools. Compiling FFEA from the source code is very 
 easy, and essentially consists of three commands: "cmake", "make", 
 and "make install". If you don't want to compile FFEA, then
 you can **download the latest x86_64 [binary release](https://bitbucket.org/FFEA/ffea/downloads/)**.
 Once FFEA has been installed, you can learn how to use it in the [tutorial](\ref Tutorial).


# Prerequisites {#prerequisites}

## Essential

* [C and C++ compilers](https://gcc.gnu.org/) (GCC >=9.4.0 or Intel >=19.1.0). There is some C++ code written using 
  the C++17 standard, and so CMake will ensure that you have a 
  recent enough compiler. FFEA will compile with GCC 9.4.0 or Intel 19.1.0,
  but we suggest using an up to date compiler to get the best performance.

* [CMake](https://cmake.org) (>=3.18). Required for building FFEA.

## Recommended

* [Python](https://www.python.org/) (3.10). Used to run the FFEA tools modules, in addition to unit and integration tests to verify that FFEA was correctly built. The [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) and [Matplotlib](https://matplotlib.org/) libraries are required.

* [PyMOL](https://pymol.org/) (>=1.8, though [1.8 is recommended](https://anaconda.org/mw/pymol)).
  Used to visualise FFEA systems and trajectories, with the plugin we provide.
    
* [NETGEN](https://sourceforge.net/projects/netgen-mesher/) 
   or [TETGEN](http://wias-berlin.de/software/tetgen/). 
     Programs which convert surface profile into volumetric meshes 
        to be used by FFEA. Essential if you want to generate meshes from
        experimental imaging data. **We suggest TETGEN!**

* [Doxygen](http://www.doxygen.org) (>=1.8). Builds the FFEA documentation. Some mathematical formulae 
     will not render correctly if [LaTeX](https://www.tug.org/texlive/) is not found.

## Optional

* [MDanalysis](https://www.mdanalysis.org/) (>=0.18.0).
    Used during rod parameterisation. If you don't plan to use KOBRA rods, you can ignore this.

* [pyPcazip](https://pypi.python.org/pypi/pyPcazip)<sup>[1](#pyPCApaper)</sup>. Some of the Python FFEA analysis tools interact with the pyPcazip 
     Principal Component Analysis libraries in order to generate standard
     PCA output(eigensystems, projections, animations etc)
     equivalent to those obtained from equivalent MD simulations.
     See [here](https://bitbucket.org/ramonbsc/pypcazip/wiki/Home) for their wiki.

* [GTS](http://gts.sourceforge.net) (>=0.7.6). The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.

* [Meshlab](http://www.meshlab.net). An open soure 
     system for processing and editing 3D triangular meshes.

## Included

* [FFEA tools](https://bitbucket.org/FFEA/ffeatools).
   A python package of command-line utilities, useful for working with FFEA, will be downloaded by CMake, but this can be adjusted with a [flag](\ref cmakeflags).

* [Boost](http://www.boost.org) (>=1.54.0). Used for ease of programming 
     at the initialisation phase. Modules "math", "algorithm" and 
     "program-options" are required. Boost 1.85.0 will be downloaded by CMake, but this can be adjusted with a [flag](\ref cmakeflags).

* [Eigen](http://eigen.tuxfamily.org) (>=3.2.10).
   FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models. Eigen 3.4.0 will be downloaded by CMake, but this can be adjusted with a [flag](\ref cmakeflags).

* [RngStreams](http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/)<sup>[2](#RngStreams1)</sup><sup>,[3](#RngStreams2)</sup>.
        Allows FFEA to safely generate random numbers when running 
        on a number of threads, as well as safe restarts, recovering the state 
        of the random number generators in the last saved time step. 

* [Tet_a_tet](https://github.com/erich666/jgt-code/blob/master/Volume_07/Number_2/Ganovelli2002/tet_a_tet.h)<sup>[4](#tetatetpaper)</sup>. 
    Used to detect element overlapping in the steric repulsion module. 

* [mtTkinter](https://github.com/RedFantom/mtTkinter/wiki) (0.4). Used in the PyMOL plugin, allowing safe threading. 

<a name="pyPCApaper">1</a>:  A Shkurti, et al., "pyPcazip: A PCA-based toolkit for compression and analysis of molecular simulation data" (2016), SoftwareX, 7:44-50. <br> 
<a name="RngStreams1">2</a>: P L'Ecuyer, "Good Parameter Sets for Combined Multiple Recursive Random Number Generators" (1999), Oper. Res., 47(1):159-164. <br> 
<a name="RngStreams2">3</a>: P L'Ecuyer et al., "An Objected-Oriented Random-Number Package with Many Long Streams and Substreams" (2002), Oper. Res., 50(6):1073-1075. <br>
<a name="tetatetpaper">4</a>:  F Ganovelli, et al., "Fast tetrahedron-tetrahedron overlap algorithm" (2002), J. Graph. Tools, 7(2):17-25.


# Configure {#configure}

FFEA uses CMake to find the compiler, dependencies and to configure files and Makefiles 
 automatically (a short introduction to CMake can be read [here](https://cmake.org/runningcmake)). 
 It is generally advisable to configure and compile 
 FFEA outside of the source tree. E.g.
 
     git clone https://bitbucket.org/FFEA/ffea.git
     cd ffea
     mkdir build
     cd build
     cmake .. [OPTIONS]

There is a list of `cmake` `[OPTIONS]` later in this section. Our favourite option
is to specify the installation directory, since the default (`/usr/local/`) 
may not be available if you do not have administrator privileges. You should also tell cmake
to use a Python 3.X interpreter if it isn't already on your `$PATH`, so that ffeatools can be correctly imported into Python:

     cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/softw/ffea -DPYTHON3_ROOT_DIR=/usr/miniconda3/envs/ffea

where `$HOME/softw/ffea` can be replaced with an installation directory of your choice.
 CMake default options seeks to build FFEA for production, and will suit most of the users.
 The following subsection gives greater detail, but if you are happy with defaults,
 you can jump to [build](\ref build).
 

# Build and Install {#build} 

After configuring you will be able to build FFEA typing:

     cmake --build . --parallel 8

Then, the following command will install FFEA:

     cmake --install .

either to the folder specified through ` -DCMAKE_INSTALL_PREFIX `
  or into a default folder (where you may need administrative privileges).


Two parallel versions of the FFEA_runner, `ffea` and `ffea_mb`, 
 as well as the ffeatools, `ffeatools`  will be found 
 in `$FFEA_HOME/bin `. Instructions on how to use them can be read 
 [here](\ref userManual) and [here](ffeamodules/html/index.html) respectively. 
 <!--In addition, the ` ffeatools ` Python package will be found in 
 `$FFEA_HOME/lib/python<version>/site-packages`.-->

If you built the documentation you will be able to read it with a web browser, 
  and so if firefox was the browser of your choice, and you installed 
  FFEA in ` $FFEA_HOME `, the command would be:

      firefox $FFEA_HOME/share/ffea/doc/html/index.html &

Install a plugin to visualise systems and trajectories in 
 [PyMOL](https://www.pymol.org). The plugin should be found in:

     $FFEA_HOME/share/ffea/plugins/pymol/FFEAplugin.tar.gz

and in order to install it, one would need to run PyMOL, and then click on
  ` Plugin ` -> ` Plugin Manager `, and on the new window, go to tab 
  ` Install New Plugin `, click ` Choose file... ` and finally find and 
  select ` FFEAplugin.tar.gz ` from your disk. You will be asked to install the 
  plugin either into a local folder or a global folder. If the 
  folder does not exist, PyMOL will ask your permission on creating the folder, 
  and will install the plugin properly. However, and just in this case,
  it will bump ` Plugin FFEAplugin has
  been installed but initialization failed `. All you need to do is 
  to restart PyMOL to use the plugin.


# Working environment {#workingEnvironment}

Executing `ffea`, `ffea_mb` and `ffeatools` is probably what most users will wish, so 
 UNIX users may find convenient to add the install folder in the ` PATH `:

      export PATH=$FFEA_HOME/bin:$PATH

In addition, in order to have direct access to the Python modules that integrate 
 the FFEA tools, and especially to those users willing to write new measure tools,
 a Python environment with the ffeatools package installed should be used.
     
ffeatools provides a suite of command-line tools to initialise FFEA systems,
via a Python API, therefore it is highly recommended for both FFEA users and developers.

When downloaded by CMake, ffeatools is automatically installed into a new Python venv with all the dependencies required by the test suite. This can be activated by navigating to `build/venv` and calling

      source bin/activate

Alternatively, you can install ffeatools into a different Python environment by entering ffeatools' source folder (e.g. `build/_deps/ffeatools-src`) and running 

      pip3 install .
      
The FFEA API may then be used within Python by running `import ffeatools`. Note: KOBRA/FFEA_rod requires this module. Information on how to use `ffeatools` can be found in the [FFEA analysis tutorial](\ref FFEAanalysistut) and the [KOBRA/rods tutorial](\ref rods).

 
# Tests {#ctest}

You may now want to check that the code was correctly compiled. 
 Do so running the provided suite of tests, either sequentially (using a single processor, 
 one tests after the other):
  
     cmake --build . --target test

or concurrently (multiple tests running independently on different processors):

     ctest -j <number-of-processes> 
 
  
# CMake options {#cmakeflags}

The following configuration flags are either fundamental to CMake or specific to FFEA:

  * `-DCMAKE_INSTALL_PREFIX=<install_dir>`       -  (default /usr/local) installation directory
  * `-DPYTHON3_ROOT_DIR=<python_env_path>`     - (e.g. /usr/miniconda3/envs/ffea) Python interpreter
  * `-DCMAKE_BUILD_TYPE=<Debug|Release|RelWithDebInfo>` -  (default Release) build type
  <!--* `-DCMAKE_CXX_COMPILER=<program>`     -  (default g++)  C++ compiler.-->

By default, CMake download versions of Boost and Eigen that it has ben tested with. Additionally, it will download the latest FFEATools (e.g. the `master` branch). If you want to use your own installations you can still do so
 through:
  
  * `-DUSE_FFEATOOLS_INTERNAL=<ON|OFF>` - default ` ON `
  * `-DUSE_BOOST_INTERNAL=<ON|OFF>` - default ` ON `
  * `-DUSE_EIGEN3_INTERNAL=<ON|OFF>` - default ` ON `
 
If you decide to do so, CMake will look for the required Boost and/or Eigen libraries.
 In the case they are not 
 installed in a standard place, you can help CMake either through: 

  * ` -DCMAKE_PREFIX_PATH="Path-to-Eigen;Path-to-Boost" `,
  * ` -DEIGEN_HOME="Path-to-Eigen" `, ` -DBOOST_ROOT="Path-to-Boost" `
  * or exporting environment variables ` EIGEN_HOME `  and ` BOOST_ROOT ` 
      to the corresponding software folders

Additional specific FFEA flags include:

  * `USE_FAST`    (default ON) will try to find the best compiler flags in terms of performance.
                               Disable ` USE_FAST ` if you intend to debug the code.
  * `USE_OPENMP`  (default ON) will enable OpenMP parallel calculations.

  * `BUILD_DOC`    (default TRY) where:

    - `NO` will not build the documentation.
    - `TRY` will try to find Doxygen and build the documentation. 
    - `YES` will try to find Doxygen (raising an error if not found) and build the documentation. 
    - `ONLY` will try to find Doxygen (raising an error if not found) and only build the documentation.




**Now that you've installed FFEA, proceed to the [tutorial](\ref Tutorial)!**