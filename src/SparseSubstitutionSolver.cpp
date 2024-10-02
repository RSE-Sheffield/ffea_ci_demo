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

#include "SparseSubstitutionSolver.h"

SparseSubstitutionSolver::SparseSubstitutionSolver() {
    // Initialise everything to zero
    num_rows = 0;
}

SparseSubstitutionSolver::~SparseSubstitutionSolver() {
    num_rows = 0;
    L_key.clear();
    U_key.clear();
    inverse_diag.clear();
    L.clear();
    U.clear();
}

int SparseSubstitutionSolver::init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) {
    // Mass matrix will have as many rows as there are nodes in the mesh
    num_rows = node.size();

    printf("Allocating memory for mass matrix...\n");
    std::vector<scalar> mass = std::vector<scalar>(num_rows * num_rows, 0);
    std::vector<scalar> mass_LU = std::vector<scalar>(num_rows * num_rows, 0);
    printf("...done.\n");

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    std::vector<int> is_pinned(node.size());
    for (int i = 0; i < node.size(); i++) {
        is_pinned[i] = 0;
    }
    for (int i = 0; i < pinned_nodes_list.size(); ++i) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }

    // build the matrix
    printf("Building the mass matrix...\n");
    int ni, nj;
    for (int n = 0; n < elem.size(); n++) {
        // add mass matrix for this element
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (i < 4 && j < 4) {
                    ni = elem[n].n[i]->index;
                    nj = elem[n].n[j]->index;
                    if (is_pinned[ni] == 0 && is_pinned[nj] == 0) {
                        if (i == j) {
                            mass[ni * num_rows + nj] += .1 * elem[n].rho * elem[n].vol_0;
                        } else {
                            mass[ni * num_rows + nj] += .05 * elem[n].rho * elem[n].vol_0;
                        }
                    } else {
                        if (i == j) {
                            mass[ni * num_rows + nj] = 1;
                        }
                    }
                } else {
                    if (i == j) {
                        ni = elem[n].n[i]->index;
                        nj = elem[n].n[j]->index;
                        mass[ni * num_rows + nj] = 1;
                    }
                }
            }
        }
    }
    printf("...done\n");

    /* Perform cholesky decomposition on the calculated mass matrix, storing the result in mass_LU */
    printf("Performing Cholesky decomposition...\n");
    for (int k = 0; k < num_rows; k++) {
        for (int i = 0; i < k + 1; i++) {
            scalar sum_ij_kj = 0;
            for (int j = 0; j < i; j++) {
                sum_ij_kj += mass_LU[INDEX(i, j)] * mass_LU[INDEX(k, j)];
            }
            mass_LU[INDEX(k, i)] = (mass[INDEX(k, i)] - sum_ij_kj) / mass_LU[INDEX(i, i)];
        }

        scalar sum_kj_2 = 0;
        for (int j = 0; j < k; j++) {
            sum_kj_2 += mass_LU[INDEX(k, j)] * mass_LU[INDEX(k, j)];
        }
        mass_LU[INDEX(k, k)] = sqrt(mass[INDEX(k, k)] - sum_kj_2);
    }
    mass.clear();
    printf("...done.\n");

    // Copy the lower matrix into upper
    for (int i = 1; i < num_rows; i++) {
        for (int j = 0; j < i; j++) {
            mass_LU[j * num_rows + i] = mass_LU[i * num_rows + j];
        }
    }

    // Allocate and fill the 'inverse_diag' array
    try {
        inverse_diag = std::vector<scalar>(num_rows);
    } catch(std::bad_alloc &) {
        FFEA_ERROR_MESSG("Failed to alloc 'inverse_diag' in SparseSubstitutionSolver::init\n");
    }

    for (int i = 0; i < num_rows; i++)
        inverse_diag[i] = 1.0 / mass_LU[i * num_rows + i];

    // Allocate the 'key' arrays
    try {
        L_key = std::vector<int>(num_rows);
        U_key = std::vector<int>(num_rows);
    } catch(std::bad_alloc &) {
        FFEA_ERROR_MESSG("Failed to alloc 'key' arrays in SparseSubstitutionSolver::init\n");
    }

    // Build the lower triangular matrix key
    int total_L = 0;
    for (int i = 0; i < num_rows; i++) {
        for (int j = num_rows - 1; j >= i; j--) {
            if (mass_LU[i * num_rows + j] != 0) {
                total_L += (j - i);
                L_key[i] = (j - i);
                break;
            }
        }
    }

    // Build the upper triangular matrix key
    total_entries_in_U = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = 0; j <= i; j++)
            if (mass_LU[i * num_rows + j] != 0) {
                total_entries_in_U += (i - j);
                U_key[i] = (i - j);
                break;
            }

    // Allocate the off-diagonal entry arrays
    try {
        L = std::vector<scalar>(total_L);
        U = std::vector<scalar>(total_entries_in_U);
    } catch(std::bad_alloc &) {
        FFEA_ERROR_MESSG("Failed to alloc off-diagonal entry arrays in SparseSubstitutionSolver::init\n");
    }

    // Fill up the L and U sparse triangular matrices
    int off_diag_data_index = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = i + 1; j <= i + L_key[i]; j++) {
            L[off_diag_data_index] = mass_LU[i * num_rows + j];
            off_diag_data_index++;
        }
    off_diag_data_index = 0;
    for (int i = 0; i < num_rows; i++)
        for (int j = i - U_key[i]; j < i; j++) {
            U[off_diag_data_index] = mass_LU[i * num_rows + j];
            off_diag_data_index++;
        }

    return FFEA_OK;
}

int SparseSubstitutionSolver::solve(std::vector<arr3> &x) {
    int i, j, index;

    // Forward substitution step Ly = b :
    index = 0;
    for (i = 0; i < num_rows; i++) {

        x[i][0] *= inverse_diag[i];
        x[i][1] *= inverse_diag[i];
        x[i][2] *= inverse_diag[i];

        for (j = 0; j < L_key[i]; j++) {
            x[i + j + 1][0] -= x[i][0] * L[j + index];
            x[i + j + 1][1] -= x[i][1] * L[j + index];
            x[i + j + 1][2] -= x[i][2] * L[j + index];
        }

        index += L_key[i];
    }

    // Backward substitution step Ux = y :
    index = total_entries_in_U - 1;
    for (i = num_rows - 1; i >= 0; i--) {

        x[i][0] *= inverse_diag[i];
        x[i][1] *= inverse_diag[i];
        x[i][2] *= inverse_diag[i];

        for (j = 0; j < U_key[i]; j++) {
            x[i - j - 1][0] -= x[i][0] * U[index - j];
            x[i - j - 1][1] -= x[i][1] * U[index - j];
            x[i - j - 1][2] -= x[i][2] * U[index - j];
        }

        index -= U_key[i];
    }

    return FFEA_OK;
}

void SparseSubstitutionSolver::apply_matrix(const std::vector<scalar> &in, std::vector<scalar> &result) {
}

