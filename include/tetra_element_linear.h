/*
 * tetra_element_linear.hpp
 */

#ifndef TETRA_ELEMENT_LINEAR_H_INCLUDED
#define TETRA_ELEMENT_LINEAR_H_INCLUDED

#include <stddef.h>
#include <stdio.h>
#include "mat_vec_fns.h"
#include "SimulationParams.h"
#include "FFEA_return_codes.h"
#include "MersenneTwister.h"
#include "SecondOrderFunctions.h"
#include "PoissonMatrixQuadratic.h"
#include "MassMatrixQuadratic.h"

class Blob;

#define NUM_NODES_QUADRATIC_TET 10

/*
 * A preprocessor macro generating code that calculates a submatrix of the bulk viscous
 * matrix. 'I' specifies what position in the matrix 'V' this 4x4 submatrix should be
 * inserted to (position I,I in the matrix).
 *
 * N.B. This macro only deals with submatrices *on the diagonal* of the bulk viscous matrix.
 */
#define BULK_VISCOUS_SUBMATRIX_DIAG(V,I) \
	V[I + 0][I + 0] = dpsi[I + 0] * dpsi[I + 0]	* (A + B); \
	V[I + 1][I + 1] = dpsi[I + 1] * dpsi[I + 1]	* (A + B); \
	V[I + 2][I + 2] = dpsi[I + 2] * dpsi[I + 2]	* (A + B); \
	V[I + 3][I + 3] = dpsi[I + 3] * dpsi[I + 3]	* (A + B); \
	V[I + 0][I + 1] = dpsi[I + 1] * dpsi[I + 0]	* (A + B); \
	V[I + 0][I + 2] = dpsi[I + 2] * dpsi[I + 0]	* (A + B); \
	V[I + 0][I + 3] = dpsi[I + 3] * dpsi[I + 0]	* (A + B); \
	V[I + 1][I + 2] = dpsi[I + 1] * dpsi[I + 2]	* (A + B); \
	V[I + 1][I + 3] = dpsi[I + 1] * dpsi[I + 3]	* (A + B); \
	V[I + 2][I + 3] = dpsi[I + 2] * dpsi[I + 3]	* (A + B);

/*
 * A preprocessor macro generating code that calculates a submatrix of the bulk viscous
 * matrix. 'I' and 'J' specify what position in the matrix 'V' this 4x4 submatrix should
 * be inserted to (position I,J in the matrix).
 *
 * N.B. This code will calculate submatrices anywhere in the bulk viscous matrix,
 * but the BULK_VISCOUS_SUBMATRIX_DIAG macro above provides simpler more efficient
 * code for the case in which the submatrix lies on the diagonal of the bulk viscous
 * matrix.
 */
#define BULK_VISCOUS_SUBMATRIX_OFFDIAG(V, I, J) \
	K[0][0] = dpsi[I] * dpsi[J]; \
	K[0][1] = dpsi[J + 1] * dpsi[I]; \
	K[0][2] = dpsi[J + 2] * dpsi[I]; \
	K[0][3] = dpsi[J + 3] * dpsi[I]; \
	\
	K[1][0] = dpsi[I + 1] * dpsi[J]; \
	K[1][1] = dpsi[I + 1] * dpsi[J + 1]; \
	K[1][2] = dpsi[I + 1] * dpsi[J + 2]; \
	K[1][3] = dpsi[I + 1] * dpsi[J + 3]; \
	\
	K[2][0] = dpsi[I + 2] * dpsi[J]; \
	K[2][1] = dpsi[I + 2] * dpsi[J + 1]; \
	K[2][2] = dpsi[I + 2] * dpsi[J + 2]; \
	K[2][3] = dpsi[I + 2] * dpsi[J + 3]; \
	\
	K[3][0] = dpsi[I + 3] * dpsi[J]; \
	K[3][1] = dpsi[I + 3] * dpsi[J + 1]; \
	K[3][2] = dpsi[I + 3] * dpsi[J + 2]; \
	K[3][3] = dpsi[I + 3] * dpsi[J + 3]; \
	\
	V[I + 0][J + 0] = (B + A) * K[0][0]; \
	V[I + 0][J + 1] = B * K[0][1] + A * K[1][0]; \
	V[I + 0][J + 2] = B * K[0][2] + A * K[2][0]; \
	V[I + 0][J + 3] = B * K[0][3] + A * K[3][0]; \
	\
	V[I + 1][J + 0] = B * K[1][0] + A * K[0][1]; \
	V[I + 1][J + 1] = (B + A) * K[1][1]; \
	V[I + 1][J + 2] = B * K[1][2] + A * K[2][1]; \
	V[I + 1][J + 3] = B * K[1][3] + A * K[3][1]; \
	\
	V[I + 2][J + 0] = B * K[2][0] + A * K[0][2]; \
	V[I + 2][J + 1] = B * K[2][1] + A * K[1][2]; \
	V[I + 2][J + 2] = (B + A) * K[2][2]; \
	V[I + 2][J + 3] = B * K[2][3] + A * K[3][2]; \
	\
	V[I + 3][J + 0] = B * K[3][0] + A * K[0][3]; \
	V[I + 3][J + 1] = B * K[3][1] + A * K[1][3]; \
	V[I + 3][J + 2] = B * K[3][2] + A * K[2][3]; \
	V[I + 3][J + 3] = (B + A) * K[3][3]; \

#define RAND(A, B) ((A) + ((B)-(A))*(rng[thread_id].rand()))

#define DPSI1_DX 0
#define DPSI2_DX 1
#define DPSI3_DX 2
#define DPSI4_DX 3
#define DPSI1_DY 4
#define DPSI2_DY 5
#define DPSI3_DY 6
#define DPSI4_DY 7
#define DPSI1_DZ 8
#define DPSI2_DZ 9
#define DPSI3_DZ 10
#define DPSI4_DZ 11

/*
 * A 10-point "quadratic" tetrahedron
 */
class tetra_element_linear {
public:

    tetra_element_linear();

    /* Properties of this element */
    scalar rho; /// density
    scalar A; /// shear viscosity
    scalar B; /// second coefficient of viscosity
    scalar G; /// shear modulus
    scalar E; /// bulk modulus
    scalar dielectric;
    scalar mass;

    /** @brief A quadratic tetrahedron has 10 nodes. Keep pointers to the actual memory location of the nodes. */
    mesh_node *n[NUM_NODES_QUADRATIC_TET];

    /** @brief The 12-vector containing the shape function derivatives for this element */
    vector12 dpsi;

    /** @brief
     * The del2 matrix for this element (in upper triangular form, since it's symmetric).
     * Used in constructing the diffusion matrix and the poisson matrix.
     */
    upper_triangular_matrix4 del2;

    PoissonMatrixQuadratic K_alpha;

    /** @brief Store the contribution from this element to the force on each of its four nodes */
    vector3 node_force[NUM_NODES_QUADRATIC_TET];

    /** @brief The rest volume of this element */
    scalar vol_0;

    /** @brief The current volume of this element */
    scalar vol;

    /** @brief The gradient deformation tensor for this element (needed for potential energy calculation) */
    matrix3 F_ij;

    /** @brief The double contraction of the internal stress tensor, including elastic and thermal stresses */
    scalar internal_stress_mag;

    /** @brief The inverse jacobian of this element at rest */
    matrix3 J_inv_0;

    /** @brief Viscosity Matrix for the internal forces of this element */
    matrix12 viscosity_matrix;

    /* Blob this element is a part of */
    Blob *daddy_blob;

    /** @brief Index of this element in the parent Blob */
    int index;

    vector3 centroid;

    /** @brief
     * Get the memory location of the specified element of K_alpha
     */
    scalar * get_K_alpha_element_mem_loc(int ni, int nj);

    /** @brief Calc the diffusion matrix for this element */
    void calculate_K_alpha();

    void construct_element_mass_matrix(MassMatrixQuadratic *M_alpha);

    void add_K_alpha(scalar *K, int num_nodes);

    /** @brief Returns the gradient of the potential at the given (s,t,u) position in the element */
    void get_grad_phi_at_stu(vector3 *grad_phi, scalar s, scalar t, scalar u);

    /** @brief Calculates the force on each node of the element due to the electrostatic potential gradient there */
    void calculate_electrostatic_forces();

    /** @brief Calculates the Jacobian matrix for this element */
    void calculate_jacobian(matrix3 J);

    /** @brief Calculates Deformation Gradient */
    void calc_deformation(matrix3 J);

    /** @brief Calculate the elastic contribution to the force */
    void calc_elastic_force_vector(vector12 F);

    /** @brief
     * Inverts the given jacobian matrix J, using this to calculate the derivatives
     * of the shape functions which are stored in dpsi. This is an array
     * of all 12 derivatives, in the following order:
     *		dpsi    =	[ d(psi)_1/dx ]
     *				[ d(psi)_2/dx ]
     *				[ d(psi)_3/dx ]
     *				[ d(psi)_4/dx ]
     *				[ d(psi)_1/dy ]
     *				[ d(psi)_2/dy ]
     *				[     ...     ]
     *				[ d(psi)_4/dz ]
     *
     * Function also (as a by-product of the inversion) calculates the volume of the
     * element whose jacobian this is, which is stored in 'vol'.
     */
    int calc_shape_function_derivatives_and_volume(matrix3 J);

    /** @brief
     *  Uses above functions to get volume. Easy
     */
    scalar calc_volume();

    /** @brief
     * Builds the viscosity matrix from the shape function derivatives, the shear and bulk
     * viscosity constants, and the element volume
     */
    void create_viscosity_matrix();

    /*
     *
     */
    void add_shear_elastic_stress(matrix3 J, matrix3 stress);

    /*
     *
     */
    void add_bulk_elastic_stress(matrix3 stress);

    /** @brief
     * Given the shape function derivatives, the element volume and a random number generator, this
     * function calculates the fluctuating stress tensor, generating a stochastic change in the
     * nodal velocities for the element under consideration. This function will add its contribution
     * to the given 12-vector du.
     *
     */
    void add_fluctuating_stress(SimulationParams *params, MTRand rng[], matrix3 stress, int thread_id);

    /** @brief
     * Applies the given stress tensor to the shape function derivatives to get the contribution to du
     */
    void apply_stress_tensor(matrix3 stress, vector12 du);

    /** @brief
     * Sets the given 12-vector to the velocities of this element's four nodes,
     */
    void get_element_velocity_vector(vector12 v);

    /** @brief
     * Add this element's nodal forces to those given in the force 12-vector
     */
    void add_element_force_vector(vector12 force);

    /** @brief Add given force to the specified node of this element */
    void add_force_to_node(int i, vector3 *f);

    /** @brief A roundabout and inefficient way of working out what node (from 0 to 9) this index corresponds to */
    int what_node_is_this(int index);

    void print();

    /** @brief
     * Applies the mass matrix (for a linear tetrahedral element of density rho and equilibrium volume vol_0)
     * to the force vector du to get the correct distribution of force between nodes.
     */
    void apply_element_mass_matrix(vector12 du);

    void volume_coord_to_xyz(scalar eta0, scalar eta1, scalar eta2, scalar eta3, vector3 *r);

    void zero_force();

    void linearise_element();

    void calc_centroid();


private:

    /** @brief
     * The last determinant of this element's transformation (used to work out whether it has inverted itself)
     */
    scalar last_det;

    /** @brief
     * Creates del2 matrix from the shape function derivatives
     */
    void calc_del2_matrix();

    void add_diffusion_matrix(matrix12 V);

    struct tetrahedron_gauss_point {
        scalar W;
        scalar eta[4];
    };

};

#endif
