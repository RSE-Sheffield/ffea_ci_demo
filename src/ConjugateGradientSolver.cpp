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

#include "ConjugateGradientSolver.h"

ConjugateGradientSolver::ConjugateGradientSolver() {
    num_rows = 0;
    epsilon2 = 0;
    i_max = 0;
    key = nullptr;
    entry = nullptr;
    preconditioner = nullptr;
    d = nullptr;
    r = nullptr;
    q = nullptr;
    s = nullptr;
    f = nullptr;
}

ConjugateGradientSolver::~ConjugateGradientSolver() {
    delete[] key;
    delete[] entry;
    delete[] preconditioner;
    delete[] d;
    delete[] r;
    delete[] q;
    delete[] s;
    delete[] f;
    key = nullptr;
    entry = nullptr;
    preconditioner = nullptr;
    d = nullptr;
    r = nullptr;
    q = nullptr;
    s = nullptr;
    f = nullptr;
    num_rows = 0;
    epsilon2 = 0;
    i_max = 0;
}

int ConjugateGradientSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list) {
    int n, i, j, ni, nj;

    // Store the number of rows, error threshold (stopping criterion for solver) and max
    // number of iterations, on this Solver (these quantities will be used a lot)
    this->num_rows = num_nodes;
    this->epsilon2 = params->epsilon2;
    this->i_max = params->max_iterations_cg;

    printf("\t\tAttempting to allocate %d scalars for mass_LU...\n", num_rows * num_rows);
    scalar *mass_LU = new(std::nothrow) scalar[num_rows * num_rows];
    if (!mass_LU) {
        FFEA_ERROR_MESSG("Could not allocate mass_LU\n");
    }
    printf("\t\t...success.\n");

    printf("\t\tZeroing...\n");
    for (i = 0; i < num_rows * num_rows; i++) {
        mass_LU[i] = 0;
    }
    printf("\t\t...done\n");

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    std::vector<int> is_pinned(num_nodes);
    for (i = 0; i < num_nodes; i++) {
        is_pinned[i] = 0;
    }
    for (i = 0; i < num_pinned_nodes; i++) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }

    // build the matrix
    scalar sum1 = 0.0, sum2 = 0.0;
    printf("\t\tBuilding the mass matrix...\n");
    for (n = 0; n < num_elements; n++) {
        // add mass matrix for this element
        for (i = 0; i < 10; i++) {
            for (j = 0; j < 10; j++) {
                if (i < 4 && j < 4) {
                    ni = elem[n].n[i]->index;
                    nj = elem[n].n[j]->index;
                    if (is_pinned[ni] == 0 && is_pinned[nj] == 0) {
                        if (i == j) {
                            mass_LU[ni * num_rows + nj] += .1 * elem[n].rho * elem[n].vol_0;
			    sum1 += .1 * elem[n].rho * elem[n].vol_0;
                        } else {
                            mass_LU[ni * num_rows + nj] += .05 * elem[n].rho * elem[n].vol_0;
			    sum1 += .05 * elem[n].rho * elem[n].vol_0;
                        }
                    } else {
                        if (i == j) {
                            mass_LU[ni * num_rows + nj] = 1;
                        }
                    }
                } else {
                    if (i == j) {
                        ni = elem[n].n[i]->index;
                        nj = elem[n].n[j]->index;
                        mass_LU[ni * num_rows + nj] = 1;
                    }
                }
            }
        }
    }
    printf("\t\t...done\n");

    for (n = 0; n < num_elements; n++) {
	sum2 += elem[n].vol_0 * elem[n].rho;
    }
    // Allocate memory for and initialise 'key' array
    key = new(std::nothrow) int[num_rows + 1];
    if (!key) FFEA_ERROR_MESSG("Failed to allocate 'key' in ConjugateGradientSolver\n"); 

    for (i = 0; i < num_rows; i++) {
        key[i] = 0;
    }

    // Get the number of non-zero entries in each row
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_rows; j++) {
            if (mass_LU[i * num_rows + j] != 0) {
                key[i]++;
            }
        }
    }

    // Sum up all the non-zero totals of each row to get the total number of non-zero entries
    // in the whole matrix, constructing the 'key' in the process
    int total_non_zeros = 0, last;
    for (i = 0; i < num_rows; i++) {
        last = key[i];
        key[i] = total_non_zeros;
        total_non_zeros += last;
    }
    key[num_rows] = total_non_zeros;

    // Allocate memory for the 'entry' array
    printf("\t\tNum non zero entries = %d.\n\t\tAllocating space...\n", total_non_zeros);
    entry = new(std::nothrow) sparse_entry[total_non_zeros];
    if (!entry) FFEA_ERROR_MESSG("Failed to allocate 'entry' in ConjugateGradientSolver\n"); 
    printf("\t\t...done.\n");

    // Fill the 'entry' array
    int entry_index = 0;
    scalar val;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_rows; j++) {
            val = mass_LU[i * num_rows + j];
            if (val != 0) {
                entry[entry_index].val = val;
                entry[entry_index].column_index = j;
                entry_index++;
            }
        }
    }
    double sum3 = 0.0;
    for(i = 0; i < total_non_zeros; ++i) {
	if(entry[i].val != 1) {
		sum3 += entry[i].val;
	}
    }
    cout << "Mass Solver Sums: " << mesoDimensions::mass*sum1 << " " << mesoDimensions::mass*sum2 << " " << mesoDimensions::mass*sum3 << endl;
    // Create the jacobi preconditioner matrix (diagonal)
    preconditioner = new(std::nothrow) scalar[num_rows];
    if (!preconditioner) FFEA_ERROR_MESSG("Failed to allocate 'preconditioner' in ConjugateGradientSolver\n"); 
    for (i = 0; i < num_rows; i++)
        preconditioner[i] = 1.0 / mass_LU[i * num_rows + i];

    // create the work vectors necessary for use by the conjugate gradient solver
    d = new(std::nothrow) arr3[num_rows];
    r = new(std::nothrow) arr3[num_rows];
    q = new(std::nothrow) arr3[num_rows];
    s = new(std::nothrow) arr3[num_rows];
    f = new(std::nothrow) arr3[num_rows];
    if (!d || !r || !q || !s || !f) FFEA_ERROR_MESSG(" Failed to create the work vectors necessary for ConjugateGradientSolver\n"); 

    delete[] mass_LU;

    return FFEA_OK;
}

int ConjugateGradientSolver::solve(arr3 *x) {
    int i = 0;

    scalar delta_new, delta_old, dTq, alpha;
    delta_new = conjugate_gradient_residual_assume_x_zero(x);
    for (i = 0; i < i_max; i++) {

        // Once convergence is achieved, return
        if (residual2() < epsilon2) {
            return FFEA_OK;
        }

        dTq = parallel_sparse_matrix_apply();

        alpha = delta_new / dTq;
        parallel_vector_add_self(x, alpha, d, num_rows);
        parallel_vector_add_self(r, -alpha, q, num_rows);

        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();

        parallel_vector_add(d, (delta_new / delta_old), s, num_rows);
    }

    // If desired convergence was not reached in the set number of iterations...
    FFEA_ERROR_MESSG("Conjugate gradient solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", i_max);
}

void ConjugateGradientSolver::apply_matrix(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] = 0;
        for (int j = key[i]; j < key[i + 1]; j++) {
            result[i] += entry[j].val * in[entry[j].column_index];
        }
    }
}

/* */
scalar ConjugateGradientSolver::conjugate_gradient_residual_assume_x_zero(arr3 *b) {
    int i;
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
// //#pragma omp parallel for default(none) private(i) shared(b) reduction(+:delta_new)
#endif
    for (i = 0; i < num_rows; i++) {
        r[i][0] = b[i][0];
        r[i][1] = b[i][1];
        r[i][2] = b[i][2];
        f[i][0] = b[i][0];
        f[i][1] = b[i][1];
        f[i][2] = b[i][2];
        b[i][0] = 0;
        b[i][1] = 0;
        b[i][2] = 0;
        d[i][0] = preconditioner[i] * r[i][0];
        d[i][1] = preconditioner[i] * r[i][1];
        d[i][2] = preconditioner[i] * r[i][2];
        delta_new += r[i][0] * d[i][0] + r[i][1] * d[i][1] + r[i][2] * d[i][2];
    }

    return delta_new;
}

/* */
scalar ConjugateGradientSolver::parallel_sparse_matrix_apply() {
    int i, j;
    scalar dTq = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
// //#pragma omp parallel for default(none) private(i, j) reduction(+:dTq)
#endif
    for (i = 0; i < num_rows; i++) {
        q[i][0] = 0;
        q[i][1] = 0;
        q[i][2] = 0;
        for (j = key[i]; j < key[i + 1]; j++) {
            q[i][0] += entry[j].val * d[entry[j].column_index][0];
            q[i][1] += entry[j].val * d[entry[j].column_index][1];
            q[i][2] += entry[j].val * d[entry[j].column_index][2];
        }

        dTq += d[i][0] * q[i][0] + d[i][1] * q[i][1] + d[i][2] * q[i][2];
    }

    return dTq;
}

void ConjugateGradientSolver::parallel_vector_add_self(arr3 *v1, scalar a, arr3 *v2, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, a, v2, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i][0] += a * v2[i][0];
        v1[i][1] += a * v2[i][1];
        v1[i][2] += a * v2[i][2];
    }
}

void ConjugateGradientSolver::parallel_vector_add(arr3 *v1, scalar a, arr3 *v2, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, a, v2, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i][0] = v2[i][0] + a * v1[i][0];
        v1[i][1] = v2[i][1] + a * v1[i][1];
        v1[i][2] = v2[i][2] + a * v1[i][2];
    }
}

/* */
scalar ConjugateGradientSolver::parallel_apply_preconditioner() {
    int i;
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
// //#pragma omp parallel for default(none) private(i) reduction(+:delta_new)
#endif
    for (i = 0; i < num_rows; i++) {
        s[i][0] = preconditioner[i] * r[i][0];
        s[i][1] = preconditioner[i] * r[i][1];
        s[i][2] = preconditioner[i] * r[i][2];
        delta_new += r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2];
    }

    return delta_new;
}

/* */
scalar ConjugateGradientSolver::residual2() {
    int i;
    scalar r2 = 0, f2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
// //#pragma omp parallel for default(none) private(i) shared(stderr) reduction(+:r2, f2)
#endif
    for (i = 0; i < num_rows; i++) {
        r2 += r[i][0] * r[i][0] + r[i][1] * r[i][1] + r[i][2] * r[i][2];
        f2 += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
    }
    if (f2 == 0.0) {
        return 0.0;
    } else {
        return r2 / f2;
    }
}

