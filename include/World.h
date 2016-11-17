#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <typeinfo>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

// #include "MersenneTwister.h"
#include "RngStream.h"
#include "NearestNeighbourLinkedListCube.h"
#include "BEM_Poisson_Boltzmann.h"
#include "BiCGSTAB_solver.h"
#include "FFEA_user_info.h"
#include "FFEA_return_codes.h"
#include "FFEA_input_reader.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"
#include "VdW_solver.h"
#include "Steric_solver.h"
#include "Steric_solverX.h"
#include "LJSteric_solver.h"
#include "PreComp_solver.h"
#include "LJ_matrix.h"
#include "BindingSite.h"
#include "Spring.h"
#include "SparseMatrixFixedPattern.h"
#include "KineticState.h"

#include "dimensions.h"
using namespace std;

class World {
public:
    World();

    ~World();

    /* */
    int init(string FFEA_script_filename, int frames_to_delete, int mode, bool writeEnergy);

    /* */
    int get_smallest_time_constants();

    /* */
    int enm(set<int> blob_indices, int num_modes);

    /* */
    int dmm(set<int> blob_indices, int num_modes);

    /* */
    int dmm_rp(set<int> blob_indices, int num_modes);

    /* */
    int run();

    /* */
    int read_and_build_system(vector<string> script_vector);

    /* */
    int load_kinetic_maps(vector<string> map_fnames, vector<int> map_from, vector<int> map_to, int blob_index);

    /* */
    int build_kinetic_identity_maps();

    /* */
    int load_kinetic_states(string states_fname, int blob_index);

    /* */
    int load_kinetic_rates(string rates_fname, int blob_index);

    /* */
    void print_kinetic_rates_to_screen(int type);

    /* */
    void get_system_CoM(vector3 *system_CoM);

    /* */
    void get_system_centroid(vector3 *centroid);

    /* */
    void get_system_dimensions(vector3 *dimenstion_vector);

    /* */
    int enm(int *blob_index, int num_modes);

    /* */
    int get_num_blobs();


private:

    /** @brief 2-D Array of Blob objects (blob i, conformation j) */
    Blob **blob_array;

    /** @brief Which conformation is active in each blob */
    Blob **active_blob_array;

    /** @brief Maps for kinetic switching of conformations */
    SparseMatrixFixedPattern ***kinetic_map;
    SparseMatrixFixedPattern ****kinetic_return_map;

    //@{
    /** @brief Kinetic State and Rate objects */
    KineticState **kinetic_state;
    scalar ***kinetic_rate;
    scalar ***kinetic_base_rate;
    //@}

    /** @brief An array of springs which connect nodes if necessary */
    Spring *spring_array;

    /** @brief And how many springs are there? */
    int num_springs;

    /** @brief How many kinetic binding sites are there? */
    int num_binding_sites;

    /** @brief Check whether mass is present anywhere, to determine whether or not to write kinetic energies to files */
    bool mass_in_system;

    /** @brief How many threads are available for parallelisation */
    int num_threads;

    /** @brief An array of pointers to random number generators (for use in parallel) */
    RngStream *rng;

    /** @brief A pointer to an array of arrays, containing the seeds of the different RNGStreams */
    unsigned long **Seeds; 

    /** @brief The number of seeds stored in Seeds. */
    int num_seeds; 

    /** @brief An array of pointers to random number generators for use in kinetics */
    RngStream *kinetic_rng;

    /** @brief Parameters being used for this simulation */
    SimulationParams params;

    /** @brief
     * Data structure keeping track of which `cell' each face lies in (where the world has been discretised into a grid of cells of dimension 1.5 kappa)
     * so that the BEM matrices may be constructed quickly and sparsely.
     */
    NearestNeighbourLinkedListCube lookup;

    /** @brief Output trajectory file */
    FILE *trajectory_out;

    /** @brief Output kinetics file */
    FILE *kinetics_out;

    /** @brief Output measurement file */
    FILE *measurement_out;

    /** @brief Output detailed measurements file. May be unneccesary */
    FILE *detailed_meas_out;

    //@{
    /** Energies */
    scalar kineticenergy, strainenergy, springenergy, **springfieldenergy, vdwenergy, preCompenergy;
    //@}

    /** Momenta */
    vector3 L;

    //@{
    /** Geometries */
    vector3 CoM, CoG;
    scalar rmsd;
    //@}

    /** @brief * Output Checkpoint file */
    FILE *checkpoint_out;

    /*
     *
     */
    //		SurfaceElementLookup surface_element_lookup;

    /** @brief BEM solver for the exterior electrostatics */
    BEM_Poisson_Boltzmann PB_solver;

    /** @brief Number of surface faces in entire system */
    int total_num_surface_faces;

    /** @brief
     * Vector of the electrostatic potential on each surface in entire system
     */
    scalar *phi_Gamma;
    scalar *J_Gamma;
    scalar *work_vec;

    /** @brief
     * Biconjugate gradient stabilised solver for nonsymmetric matrices
     */
    BiCGSTAB_solver nonsymmetric_solver;

    /** Van der Waals solver */
    VdW_solver *vdw_solver;

    /** @brief LJ parameters matrix */
    LJ_matrix lj_matrix;

    /** @brief Binding Interactions matrix */
    BindingSite_matrix binding_matrix;


    /**
      * @brief stores info within the <precomp> block at the .ffea file.
      */
    PreComp_params pc_params; 
    /**
      * @brief PreComputed potentials solver
      */
    PreComp_solver pc_solver;


    vector3 box_dim;

    long long step_initial;

    int load_springs(const char *fname);

    void activate_springs();

    int apply_springs();

    scalar get_spring_field_energy(int index0, int index1);

    /** @brief calculates the kinetic rates as a function of the energy of the system*/
    int calculate_kinetic_rates();

    /** @brief randomly chooses a new kinetic state based upon the kinetic rates / switching probabilitie */
    int choose_new_kinetic_state(int blob_index, int *target);

    /** @brief changes the kinetic state based upon the kinetic rates. Maps between conformations and adds/ removes bound sites */
    int change_kinetic_state(int blob_index, int target_state);

    int get_next_script_tag(FILE *in, char *buf);

    void apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N);

    void do_es();

    void make_trajectory_from_eigenvector(string traj_out_fname, int blob_index, int mode_index, Eigen_VectorX evec, scalar step);

    void print_evecs_to_file(string fname, Eigen_MatrixX ev, int num_rows, int num_modes);

    void print_evals_to_file(string fname, Eigen_VectorX ev, int num_modes);

    void write_eig_to_files(scalar *evals_ordered, scalar **evecs_ordered, int num_modes, int num_nodes);

    void write_output_header(FILE *fout, string fname);

    void print_trajectory_and_measurement_files(int step, scalar wtime);

    void make_measurements();

    void write_measurements_to_file(FILE *fout, int step);

    void write_detailed_measurements_to_file(FILE *fout);

    void print_trajectory_conformation_changes(FILE *fout, int step, int *from_index, int *to_index);

    void print_kinetic_files(int step);

    void print_static_trajectory(int step, scalar wtime, int blob_index);
};

#endif
