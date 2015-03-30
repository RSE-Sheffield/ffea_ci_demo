#ifndef MASSLUMPEDSOLVER_HPP_INCLUDED
#define MASSLUMPEDSOLVER_HPP_INCLUDED


#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseMatrixTypes.h"

class MassLumpedSolver : public Solver {
public:

    /* Constructor */
    MassLumpedSolver();

    /* Destructor */
    ~MassLumpedSolver();

    /* Builds the diagonal mass matrix and gets reciprocal of each value */
    int init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list);

    /* Applies inverse mass matrix (since diagonal: Mx = f => x = f_i/M_i */
    int solve(vector3 *x);

    /* Applies the mass matrix to the given vector, 'in', putting the result in 'result'*/
    void apply_matrix(scalar *in, scalar *result);

private:

    /* Number of rows in original matrix */
    int num_rows;

    /* Mass matrix diagonal inversed */
    scalar *inv_M;
};

#endif
