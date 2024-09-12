//
//  This file is part of the FFEA simulation package
//
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file.
//
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
//
//  To help us fund FFEA development, we humbly ask that you cite
//  the research papers on the package.
//

/*
 *      Blob.h
 */

#ifndef BLOB_H_INCLUDED
#define BLOB_H_INCLUDED

#include "SparsityPattern.h"
#include "ConnectivityTypes.h"

#include <cstdio>
#include <set>
#include <memory>
#include <omp.h>
#include <algorithm>  // std::find
#include <Eigen/Sparse>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "ConjugateGradientSolver.h"
#include "MassLumpedSolver.h"
#include "NoMassCGSolver.h"
#include "Face.h"
#include "CG_solver.h"
#include "BEM_Poisson_Boltzmann.h"
#include "LJ_matrix.h"
#include "BindingSite.h"
#include "PreComp_solver.h"
#include "dimensions.h"

#ifdef USE_DOUBLE_LESS
typedef Eigen::MatrixXf Eigen_MatrixX;
typedef Eigen::VectorXf Eigen_VectorX;
typedef Eigen::Matrix3f Eigen_Matrix3;
typedef Eigen::Vector3f Eigen_Vector3;
#else
typedef Eigen::MatrixXd Eigen_MatrixX;
typedef Eigen::VectorXd Eigen_VectorX;
typedef Eigen::Matrix3d Eigen_Matrix3;
typedef Eigen::Vector3d Eigen_Vector3;
#endif


struct Blob_conf {

   // initial placement, orientation and velocity:
   int set_centroid;
   scalar centroid[3];
   int set_velocity; 
   scalar velocity[3];
   int set_rotation; 
   int rotation_type;
   scalar rotation[9];

   // kinetics:
   string states, rates;
   vector<string> maps; 
   vector<int> maps_conf_index_to, maps_conf_index_from; 
   
};


/*
 * The "Blob" class
 */
class Blob {
public:

    /**
     * Blob constructor:
     * Initialises all variables and pointers to 0 (or nullptr). Does not perform any memory allocation. Actual Blob initialisation
     * is carried out by the init() method.
     */
    Blob();

    /**
     * Blob destructor:
     * Closes all open files, deallocates memory held by arrays, solvers etc. and sets everything to zero (nullptr)
     */
    ~Blob();

    /**
     * Allocates the node, element and embedded charge arrays, initialised with the data
     * read from the given node, element and embedded charge files.
     * 'linear_solver' sets which type of linear solver is to be used: 0 for direct (forward/backward
     * substitution) and 1 for iterative (preconditioned gonjugate gradient).
     * Also takes the simulation parameters and the array of RNGs (for multiprocessor runs).
     */
    int config(const int blob_index, const int conformation_index, const string& node_filename,
               const string& topology_filename, const string& surface_filename, const string& material_params_filename,
               const string& stokes_filename, const string& ssint_filename, const string& pin_filename,
               const string& binding_filename, const string& beads_filename, scalar scale, scalar calc_compress,
             scalar compress, int linear_solver, int blob_state, const SimulationParams &params,
             const PreComp_params &pc_params, SSINT_matrix *ssint_matrix,
             BindingSite_matrix *binding_matrix, RngStream rng[]);
    int init();

    /**
     * Calculates all internal forces on the finite element mesh, storing them on the elements
     * This requires the jacobian for each element, so if the mesh has inverted anywhere we will find out within this function
     */
    int update_internal_forces();

    /**
     * Checks whether any of the elements have inverted or not
     */
    int check_inversion();

    /**
     * Solves the EOM on the finite element mesh, updating the node positions and velocities by one time step
     */
    int update_positions();

    /**
     * If the system changes mid run (binding event, say) we may need to reinitialise the solver
     */
    int reset_solver();

    /**
      * Translates the linear nodes, then linearises the secondary nodes
      */
    void translate_linear(arr3 *vec);

    /**
     * Calculates the centroid of this Blob, then brings the Blob to the origin,
     * rotates all nodes in the Blob, and brings back the Blob to the initial position.
     * If beads = 1, then it rotates its own "bead_positions" too.
     */
    void rotate(float r11, float r12, float r13, float r21, float r22, float r23, float r31, float r32, float r33, bool beads=false);

    /**
     *   Performs rotation about x axis, then y axis, then z axis
     */
    void rotate(float xang, float yang, float zang, bool beads=false);

    /**
     * Calculates the centroid of this Blob, then translates all nodes in the Blob
     * so that the new centroid position is at the given (x,y,z) position, while
     * returning a vector with the displacement (dx, dy, dz) applied to every node.
     */
    arr3 position(scalar x, scalar y, scalar z);

    /**
     * Moves the beads according to (dx, dy, dz).
     * The name is to be related to "arr3 position(scalar x, scalar y, scalar z).
     */
    void position_beads(scalar x, scalar y, scalar z);

    /**
     * Beads are only useful before PreComp_solver.init is called.
     * They can be removed later on.
     */
    int forget_beads();

    /**
     * Add nodes to the face objects if and only if this blob is STATIC
     *
     */
     void add_steric_nodes();

    /**
     * Translate the Blob by the given vector
     */
    void move(scalar dx, scalar dy, scalar dz);

    /**
     * Calculates and returns the centre of mass of this Blob
     */
    void get_CoM(arr3 &com);

    /**
     * Calculates and returns the centroid of this Blob
     */
    void get_centroid(arr3 &com);
    void calc_and_store_centroid(arr3 &com);
    arr3 calc_centroid();

    void set_pos_0();
    void kinetically_set_faces(bool state);
    /**
     * Writes a new node file for the case of an initially translated STATIC blob, so viewer doesn't read straight from node file
     */
    int create_viewer_node_file(const char *node_filename, scalar scale);

    /**
     * Dumps all the node positions (in order) from the node array to the given file stream.
     */
    void write_nodes_to_file(FILE *trajectory_out);

    /**
     * Dumps all the node positions (in order) from the node array to the given file stream in two steps.
     */
    void pre_print(); 
    void write_pre_print_to_file(FILE *trajectory_out); 
    int toBePrinted_conf[2]{}; 
    int toBePrinted_state[2]{}; 

    /**
     * Reads the node positions from the given trajectory file stream.
     * This is useful for restarting simulations from trajectory files.
     */
    int read_nodes_from_file(FILE *trajectory_out);

    /**
     * Takes measurements of system properties: KE, PE, Centre of Mass and Angular momentum etc
     */
    void make_measurements();

    /**
     * Writes only the detailed measurements local to this blob to file!
     */
    void write_measurements_to_file(FILE *fout);

    /**
     * Calculates the current jacobian and elasticity properties of the structure
     */
    int calculate_deformation();

    scalar calc_volume();

    void make_stress_measurements(FILE *stress_out, int blob_number);

    /* DEPRECATED
     *   Will be removed.
    void calculate_vdw_bb_interaction_with_another_blob(FILE *vdw_measurement_out, int other_blob_index);
     */

    /**
     * Get the centroid of all faces on the blob surface
     */
    void calc_centroids_and_normals_of_all_faces();

    /**
     * Get the centroid of all faces and elements in the blob
     */
    void calc_all_centroids();

    /*
     *
     */
    int get_num_faces();

    /**
     * Return pointer to the ith Face of this Blob's surface
     */
    Face *get_face(int i);

    Face *absolutely_get_face(int i);

    /**
     * Return pointer to the ith Element of this Blob's surface
     */
    tetra_element_linear *get_element(int i);

    /** get_bead_position i onto arr3 v [precomp] */
    void get_bead_position(int i, arr3 &v);

    /** get the pointer to "bead_type"  [precomp] */
    std::vector<int> &get_bead_types();

    /** get_bead_type [precomp] */
    int get_bead_type(int i) const;

    /**
     * @brief returns the list of nodes where bead i should be assigned to.
     *
     * @ingroup FMM
     **/
    std::vector<int> &get_bead_assignment(int i);


    scalar get_ssint_area();

    //		/*
    //		 *
    //		 */
    //		int get_num_surface_elements();


    /**
     * Solves the poisson equation inside this Blob for a given fixed surface potential, and calculates the surface flux out of the protein
     */
    int solve_poisson(scalar *phi_gamma_IN, scalar *J_Gamma_OUT);

    /**
     * Apply the constant forces onto the corresponding nodes;
     */
    int apply_ctforces();


    /**
     * Set all forces on the blob to zero
     */
    void zero_force();

    void set_forces_to_zero();

    // arr3 get_node(int index);
    // std::array<scalar,3> get_node(int index);
    void get_node(int index, arr3 &v);
    
    void get_node_0(int index, arr3 &v);

    void copy_node_positions(std::vector<arr3*> &nodes);

    std::vector<arr3*> &get_actual_node_positions();

    void set_node_positions(const std::vector<arr3*> &node_pos);

    void add_force_to_node(const arr3& f, int index);

    // void zero_vdw_bb_measurement_data(); // DEPRECATED

    // void zero_vdw_xz_measurement_data(); // DEPRECATED

    /**
     * Set all nodes on the Blob to the given velocity vector
     */
    void velocity_all(scalar vel_x, scalar vel_y, scalar vel_z);

    /**
     * Constructs the Poisson matrix for this Blob.
     */
    void build_poisson_matrices();

    //		int elements_are_connected(int e1, int e2);

    /**
     * Builds a global viscosity matrix for this blob
     */
    int build_linear_node_viscosity_matrix(Eigen::SparseMatrix<scalar> *K);

    /**
     * Builds a global diffusion matrix for this blob based on the work of Rotne and Prager (1969)
     */
    int build_linear_node_rp_diffusion_matrix(Eigen_MatrixX *D);

    /**
     * Linearises the elasticity vector and build a global elasticity matrix for this blob
     */
    int build_linear_node_elasticity_matrix(Eigen::SparseMatrix<scalar> *A);

    /**
     * Build the mass distribution matrix for this blob
     */
    int build_linear_node_mass_matrix(Eigen::SparseMatrix<scalar> *M);

    /**
     * Returns the total mass of this Blob.
     */
    scalar get_mass();

    /**
     * Applies the WALL_TYPE_HARD boundary conditions simply by finding all nodes that have "passed through" the wall,
     * then zeroing any component of those nodes' velocities that points into the wall. This allows the Blob to "wobble"
     * its way back out of the wall (and effectively prevents any penetration larger than dt * largest velocity,
     * generally very small for sensible dt).
     */
    void enforce_box_boundaries(arr3 &box_dim);

    /* DEPRECTATED
     *   Will be removed.
     * Set the interaction flag for all faces of this Blob back to false (i.e. not interacting)
    void reset_all_faces();
     */

    void linearise_elements();

    void linearise_force();



    /**compresses blob by compression factor specified in input script*/
    void compress_blob(scalar compress);

    int get_num_nodes();

    int get_num_elements();

    int get_motion_state();

    scalar get_scale();
    
    scalar get_RandU01(); 

    int get_num_linear_nodes();

    int get_num_beads();
    bool is_using_beads();

    int getNumBindingSites();

    scalar get_rmsd();

    int get_linear_solver();

    // std::array<scalar,3> get_CoG();
    // arr3 get_CoG();
    void get_stored_centroid(arr3 &cog); 

    int get_conformation_index();
    int get_previous_conformation_index();
    void set_previous_conformation_index(int index);
    int get_state_index();
    void set_state_index(int index);
    int get_previous_state_index();
    void set_previous_state_index(int index);
    BindingSite* get_binding_site(int index);

    scalar calculate_strain_energy();

    void get_min_max(arr3 &blob_min, arr3 &blob_max);

    /* Blob, conformation and state indices */
    int blob_index;
    int conformation_index, previous_conformation_index;
    int state_index, previous_state_index;

    /*
     *
     */
    //void kinetic_bind(int site_index);
    //void kinetic_unbind(int site_index);

    /** Activates binding sites by adding nodes to the list */
    void pin_binding_site(set<int> node_indices);

    /** Deactivates binding sites by removing nodes to the list */
    void unpin_binding_site(set<int> node_indices);

    void print_node_positions();
    void print_bead_positions();
    bool there_is_mass();
    void set_springs_on_blob(bool state);
    bool there_are_springs();
    bool there_are_beads();
    bool there_is_ssint();

    scalar get_kinetic_energy();
    scalar get_strain_energy();

    int pbc_count[3]{};

private:

    /** Total number of surface elements in Blob */
    int num_surface_elements;

    /** Total number of interior elements in Blob */
    int num_interior_elements;

    /** Total number of nodes on Blob surface */
    int num_surface_nodes;

    /** Total number of nodes on Blob interior */
    int num_interior_nodes;

    /** Number of ctforces to be applied to this Blob */
    int num_l_ctf;
    int num_r_ctf;
    int num_sltotal_ctf;
    int num_slsets_ctf; ///< number of surface sets, corresponding to the length of the ctf_sl_forces, num_slsurf_ctf array

    /** Number of faces in every surface set: */
    std::vector<int> ctf_slsurf_ndx;

    /** Whether this Blob is DYNAMIC (movable; dynamics simulated) or STATIC (fixed; no simulation of dynamics; Blob is a perfectly solid object fixed in space)*/
    int blob_state;

    /** Total mass of Blob */
    scalar mass;

    /** Total ssint energy between blobs */
    //scalar ssint_bb_energy{};

    /** Array of nodes */
    std::vector<mesh_node> node;

    /** Array of node positions only */
    std::vector<arr3 *> node_position;

    /** Array of elements */
    std::vector<tetra_element_linear> elem;

    /** Array of surface faces */
    std::vector<Face> surface;

    /** List of fixed ('pinned') nodes */
    std::vector<int> pinned_nodes_list;

    /** Additional pinned node list for binding processes */
    set<int> bsite_pinned_nodes_list;

    // Interacting beads within this blob, not tracked beyond PreComp_solver
    /** Array with bead positions xyzxyzxyz.... [precomp]
      *   will be nullptr after info is loaded into PreComp_solver */
    std::vector<arr3> bead_position;

    /** 2D vector with the set of nodes where every bead should be assigned to.
      *   It will be removed after PreComp_solver is initialised [precomp] */
    std::vector<std::vector<int>> bead_assignment;

    /** Array with bead types [precomp]
      *   will be empty after info is loaded into PreComp_solver */
    std::vector<int> bead_type;

    /** Array with the nodes having linear ctforces assigned */
    std::vector<int> ctf_l_nodes;
    /** Array with the nodes having rotational ctforces assigned */
    std::vector<int> ctf_r_nodes;
    /** Array with the faces having linear ctforces assigned */
    std::vector<int> ctf_sl_faces;
    /** array with the number of faces in every surface set. */
    std::vector<int> ctf_sl_surfsize;

    /** Array with the linear ctforces: FxFyFzFxFyFz...,
      * being Fx, Fy, Fz the components of the force */
    std::vector<scalar> ctf_l_forces;
    /** Array with the magnitude of the rotational ctforces: FFF..., */
    std::vector<scalar> ctf_r_forces;
    /** Array with the rotational axis (given with point + unit vector)
      *  for ctforces: XYZxyzXYZxyzFxFyFz...,
      *  or BlobConfNodeBlobConfNode,BlobConfNodeBlobConfNode,...
      *  if using two nodes to define the axis.  */
    std::vector<scalar> ctf_r_axis;
    /** Array with the type of rotation force, 2 chars per node:
      *  where the first one can be n or p, depending of the axis defined by nodes or points
      *   and the second one can be f or t, depending on applying ctforce or cttorque.*/
    std::vector<char> ctf_r_type;
    /** Array with the linear surface ctforces: FxFyFzFxFyFz...,
      * being Fx, Fy, Fz the components of the force */
    std::vector<scalar> ctf_sl_forces;

    /** Strings of all the files that contain input data: */
    string s_node_filename, s_topology_filename, s_surface_filename, 
          s_material_params_filename, s_stokes_filename, s_ssint_filename, 
          s_pin_filename, s_binding_filename, s_beads_filename;


    /** Scale of the input coordinates to m: */
    scalar scale{};

    /** Compression stuff: */
    scalar calc_compress{}, compress{}; 

    /** Class containing simulation parameters, such as the time step, dt */
    SimulationParams params;
    PreComp_params pc_params;

    /** A pointer to the same binding matrix configured in World */ 
    BindingSite_matrix *binding_matrix{};

    /** pointer to the ssint forcefield parameters (for some energy calcs) */
    SSINT_matrix *ssint_matrix{};

    /** A pointer to whatever Solver is being used for this Blob (eg SparseSubstitutionSolver
     * or ConjugateGradientSolver). The Solver solves the equation Mx = f where M is the
     * mass matrix of this Blob and f is the force vector.
     */
    std::unique_ptr<Solver> solver;

    /** Remember what type of solver we are using */
    int linear_solver;

    /** And whether or not there is mass in this system */
    bool mass_in_blob;

    /** Are there springs on this blob? */
    bool springs_on_blob;

    /** Are there ssint on this blob? */
    bool ssint_on_blob;

    /** Are the preComp beads on this blob? */
    bool beads_on_blob;

    /** The Blob force vector (an array of the force on every node) */
    std::vector<arr3> force;

    /** The array of random number generators (needed for parallel runs) */
    RngStream *rng;

    //@{
    /** Energies */
    scalar kineticenergy{}, strainenergy{};
    //@}

    /** Momenta */
    arr3 L;

    //@{
    /** Geometries */
    arr3 CoM, CoG, CoM_0, CoG_0;
    scalar rmsd{};
    //@}

    std::unique_ptr<CG_solver> poisson_solver;
    std::shared_ptr<SparseMatrixFixedPattern> poisson_surface_matrix;
    std::shared_ptr<SparseMatrixFixedPattern> poisson_interior_matrix;
    std::vector<scalar> phi_Omega;
    std::vector<scalar> phi_Gamma;
    std::vector<scalar> q;
    //std::vector<scalar> nodal_q;
    std::vector<scalar> poisson_rhs;

    std::vector<int> num_contributing_faces;

    connectivity_entry *element_connectivity_table{};

    /*
     * Mass Matrix
     * @see build_mass_matrix()
     */
    std::shared_ptr<SparseMatrixFixedPattern> M;

    std::vector<BindingSite> binding_site;
    
    /*
     */
    std::vector<scalar> toBePrinted_nodes; 

    /**
     * Opens and reads the given 'ffea node file', extracting all the nodes for this Blob.
     * Records how many of these are surface nodes and how many are interior nodes.
     */
    int load_nodes(const char *node_filename, scalar scale);

    /**
     * Opens and reads the given 'ffea topology file', extracting all the elements for this Blob.
     * Records how many of these are surface elements and how many are interior elements.
     */
    int load_topology(const char *topology_filename);

    /**
     * Opens and reads the given 'ffea surface file', extracting all the faces for this Blob.
     */
    int load_surface(const char *surface_filename);


    /**
     * Opens and reads the given 'ffea surface file', extracting all the faces for this Blob, ignoring topology.
     */
    int load_surface_no_topology(const char *surface_filename);

    /**
     * Opens and reads the given 'ffea material params file', extracting the material parameters for each element.
     */
    int load_material_params(const char *material_params_filename);

    /**
     * Opens and reads the given 'ffea stokes params file', extracting the stokes radii for each node in the Blob.
     */
    int load_stokes_params(const char *stokes_filename, scalar scale);


    /**
     * Opens and reads the given 'ffea ssint file', extracting all the van der waals species for each face of this Blob.
     */
    int load_ssint(const char *ssint_filename, int num_ssint_face_types, string ssint_method);

    /**
     * Opens and reads the given 'ffea beads file', extracting all the beads types and positions and for this Blob.
     */
    int load_beads(const char *beads_filename, scalar scale);


    /**
     * Opens and reads the given 'ffea ctforces file', and assigns the constant forces onto nodes for this Blob.
     */
    int load_ctforces(const string& ctforces_fname);


    /**
     * Opens and reads the given 'ffea binding site file', extracting all the kinetic binding sites (types and face lists) for this Blob.
     */
    int load_binding_sites(); // const char *binding_filename, int num_binding_site_types);

    /**
     * Opens and reads the given 'ffea pinned nodes file', extracting all the faces for this Blob.
     */
    int load_pinned_nodes(const char *pin_filename);

    /**
     * Creates a new pinned nodes list from a given set
     */
    int create_pinned_nodes(set<int> list);

    /**
     * Calculate some quantities such as the rest jacobian, rest volume etc.
     */
    void calc_rest_state_info();

    /*
     *
     */
    int aggregate_forces_and_solve();

    /*
     *
     */
    void euler_integrate();

    /*
     *
     */
    int calculate_node_element_connectivity();

    int build_mass_matrix();
};

#endif
