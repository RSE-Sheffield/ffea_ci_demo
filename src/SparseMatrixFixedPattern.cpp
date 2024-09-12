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

#include "SparseMatrixFixedPattern.h"

#include <vector>

SparseMatrixFixedPattern::SparseMatrixFixedPattern() {
    num_rows = 0;
    num_nonzero_elements = 0;
    entry = {};
    key = {};
    source_list = {};
    diagonal = {};
}

SparseMatrixFixedPattern::~SparseMatrixFixedPattern() {
    entry.clear();
    key.clear();
    source_list.clear();
    diagonal.clear();
    num_rows = 0;
    num_nonzero_elements = 0;
}

int SparseMatrixFixedPattern::init(int num_rows, int num_nonzero_elements, std::vector<sparse_entry> &&entry, std::vector<int> &&key, std::vector<sparse_entry_sources> source_list) {
    this->num_rows = num_rows;
    this->num_nonzero_elements = num_nonzero_elements;
    this->entry = entry;
    this->key = key;
    this->source_list = source_list;

    // Work out which elements are on the diagonal
    try {
        this->diagonal = std::vector<scalar*>(num_rows);
    } catch(std::bad_alloc &) {
        FFEA_ERROR_MESSG("Failed to allocate 'diagonal' array in SparseMatrixFixedPattern::init\n");
    }

    for (int i = 0; i < num_rows; i++) {
        for (int j = key[i]; j < key[i + 1]; j++) {
            if (i == entry[j].column_index) {
                diagonal[i] = &this->entry[j].val;
            }
        }
    }

    return FFEA_OK;
}

// Initialise matrix without a source list (doesn't need rebuilding)
int SparseMatrixFixedPattern::init(int num_rows, const std::vector<scalar> &entries, std::vector<int> &&key, const std::vector<int> &col_indices) {
    this->num_rows = num_rows;
    this->num_nonzero_elements = entries.size();
    this->key = key;
    try {
        this->entry = std::vector<sparse_entry>(entries.size());
    } catch(std::bad_alloc &) {
        FFEA_ERROR_MESSG("Failed to allocate 'entry' in SparseMAtrixFixedPattern::init\n");
    }
    for(int i = 0; i < entries.size(); ++i) {
        entry[i].column_index = col_indices[i];
        entry[i].val = entries[i];
    }
    return FFEA_OK;
}

/* Reconstruct the matrix by adding up all the contributions from the sources stored in the source list */
void SparseMatrixFixedPattern::build() {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) private(i)
#endif
    for (i = 0; i < num_nonzero_elements; i++) {
        entry[i].val = source_list[i].sum_all_sources();
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result' */
void SparseMatrixFixedPattern::apply(const std::vector<scalar> &in, std::vector<scalar> &result) const {
    for (int i = 0; i < num_rows; i++) {
        result[i] = 0;
        for (int j = key[i]; j < key[i + 1]; j++) {
            result[i] += entry[j].val * in[entry[j].column_index];
        }
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
/* Designed for use in NoMassCGSolver */
/*void SparseMatrixFixedPattern::apply(const std::vector<arr3> &in, std::vector<arr3> &result) const {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(result, in)
#endif
    for (int i = 0; i < num_rows / 3; ++i) {
        result[i][0] = 0;
        result[i][1] = 0;
        result[i][2] = 0;
        for (int j = key[3 * i]; j < key[(3 * i) + 1]; ++j) {
            if (entry[j].column_index % 3 == 0) {
                result[i][0] += entry[j].val * in[entry[j].column_index / 3][0];
            } else if (entry[j].column_index % 3 == 1) {
                result[i][0] += entry[j].val * in[entry[j].column_index / 3][1];
            } else if (entry[j].column_index % 3 == 2) {
                result[i][0] += entry[j].val * in[entry[j].column_index / 3][2];
            }
        }
        for (int j = key[(3 * i) + 1]; j < key[(3 * i) + 2]; ++j) {
            if (entry[j].column_index % 3 == 0) {
                result[i][1] += entry[j].val * in[entry[j].column_index / 3][0];
            } else if (entry[j].column_index % 3 == 1) {
                result[i][1] += entry[j].val * in[entry[j].column_index / 3][1];
            } else if (entry[j].column_index % 3 == 2) {
                result[i][1] += entry[j].val * in[entry[j].column_index / 3][2];
            }
        }
        for (int j = key[(3 * i) + 2]; j < key[(3 * i) + 3]; ++j) {
            if (entry[j].column_index % 3 == 0) {
                result[i][2] += entry[j].val * in[entry[j].column_index / 3][0];
            } else if (entry[j].column_index % 3 == 1) {
                result[i][2] += entry[j].val * in[entry[j].column_index / 3][1];
            } else if (entry[j].column_index % 3 == 2) {
                result[i][2] += entry[j].val * in[entry[j].column_index / 3][2];
            }
        }
    }
}*/

/* Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
/* Designed for use in NoMassCGSolver */
void SparseMatrixFixedPattern::apply(const std::vector<arr3> &in, std::vector<arr3> &result) const {
    // To get rid of conditionals, define an array 'num_rows' long, and copy into result at end
    vector<scalar> work_in(num_rows);
    vector<scalar> work_result(num_rows);

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel default(none) shared(result, work_result, in, work_in)
    {
    #pragma omp for 
#endif
    for(int i = 0; i < num_rows / 3; ++i) {
        work_in[3 * i] = in[i][0];
        work_in[3 * i + 1] = in[i][1];
        work_in[3 * i + 2] = in[i][2];
    }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp for 
#endif
    for (int i = 0; i < num_rows; i++) {
        // Zero array first
        work_result[i] = 0.0;

        for(int j = key[i]; j < key[i + 1]; ++j) {
            work_result[i] += entry[j].val * work_in[entry[j].column_index];
        }
    }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp for 
#endif
    for(int i = 0; i < num_rows / 3; ++i) {
        result[i][0] = work_result[3 * i];
        result[i][1] = work_result[3 * i + 1];
        result[i][2] = work_result[3 * i + 2];
    }
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    }
#endif
}

/* Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
/* Each element applies to whole vector */
/* Designed to apply sparse matrix (kinetic map) to list of node positions for conformation changes */
void SparseMatrixFixedPattern::block_apply(const std::vector<arr3> &in, std::vector<arr3> &result) const {
    for(int i = 0; i < num_rows; ++i) {
        result[i][0] = 0;
        result[i][1] = 0;
        result[i][2] = 0;
        for(int j = key[i]; j < key[i + 1]; ++j) {
            result[i][0] += entry[j].val * in[entry[j].column_index][0];
            result[i][1] += entry[j].val * in[entry[j].column_index][1];
            result[i][2] += entry[j].val * in[entry[j].column_index][2];
        }
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
/* Each element applies to whole vector */
/* Designed to apply sparse matrix (kinetic map) to list of node positions for conformation changes */
void SparseMatrixFixedPattern::block_apply(const std::vector<arr3*> &in, std::vector<arr3*> &result) const {
    int i, j;
    for(i = 0; i < num_rows; ++i) {
        (*result[i])[0] = 0;
        (*result[i])[1] = 0;
        (*result[i])[2] = 0;
        for(j = key[i]; j < key[i + 1]; ++j) {
            (*result[i])[0] += entry[j].val * (*in[entry[j].column_index])[0];
            (*result[i])[1] += entry[j].val * (*in[entry[j].column_index])[1];
            (*result[i])[2] += entry[j].val * (*in[entry[j].column_index])[2];
        }
    }
}

/* Applies this matrix to the given vector 'in', also writing the result to 'in'. 'in' is made of 'arr3's */
/* Designed to apply sparse matrix (kinetic map) to list of node positions for conformation changes */
void SparseMatrixFixedPattern::block_apply(std::vector<arr3*> &in) const {
    int i, j;
    vector<arr3> result(num_rows);

    for(i = 0; i < num_rows; ++i) {
        result[i][0] = 0;
        result[i][1] = 0;
        result[i][2] = 0;
        for(j = key[i]; j < key[i + 1]; ++j) {
            result[i][0] += entry[j].val * (*in[entry[j].column_index])[0];
            result[i][1] += entry[j].val * (*in[entry[j].column_index])[1];
            result[i][2] += entry[j].val * (*in[entry[j].column_index])[2];
        }
    }
    for(i = 0; i < num_rows; ++i) {
        (*in[i])[0] = result[i][0];
        (*in[i])[1] = result[i][1];
        (*in[i])[2] = result[i][2];
    }
}

std::shared_ptr<SparseMatrixFixedPattern> SparseMatrixFixedPattern::apply(std::shared_ptr<SparseMatrixFixedPattern> &in) const {
    // Build big matrix first, sparse it up later
    int num_rows_A = num_rows;
    int num_rows_B = in->get_num_rows();
    int num_rows_result = num_rows_A;

    // Get num_columns_result
    int num_columns_result = in->get_num_columns();
    std::vector<std::vector<scalar>> result_dense(num_rows_result, std::vector<scalar>(num_columns_result, 0.0));

    // Get in matrix pointers for quick access
    const std::vector<sparse_entry> &in_entry = in->get_entries();
    const std::vector<int> &in_key = in->get_key();

    // Make big result matrix of doom
    for(int i = 0; i < num_rows_A; ++i) {
        for(int j = key[i] ; j < key[i + 1]; ++j) {
            for(int k = 0; k < num_rows_B; ++k) {
                for(int l = in_key[k] ; l < in_key[k + 1]; ++l) {
                    if(entry[j].column_index == k) {
                        result_dense[i][in_entry[l].column_index] += entry[j].val * in_entry[l].val;
                    }
                }
            }
        }
    }

    // Build sparse matrix from big matrix
    int num_entries_result = 0;
    for(int i = 0; i < num_rows_result; ++i) {
        for(int j = 0; j < num_columns_result; ++j) {
            if(fabs(result_dense[i][j]) >= 0.001) {
                num_entries_result++;
            }
        }
    }

    std::vector<scalar>entries_result = std::vector<scalar>(num_entries_result);
    std::vector<int> key_result = std::vector<int>(num_rows_result + 1);
    std::vector<int>col_indices_result = std::vector<int>(num_entries_result);

    int l = 0;
    key_result[0] = 0;
    for(int i = 0; i < num_rows_result; ++i) {
        for(int j = 0; j < num_columns_result; ++j) {
            if(fabs(result_dense[i][j]) >= 0.001) {
                entries_result[l] = result_dense[i][j];
                col_indices_result[l] = j;
                l++;
            }
        }
        key_result[i + 1] = l;
    }

    std::shared_ptr<SparseMatrixFixedPattern> result_sparse = std::make_shared<SparseMatrixFixedPattern>();
    // std::move() is used here to convert key_result to a rval, allowing it's ownership to be transferred to init()
    result_sparse->init(num_rows_result, entries_result, std::move(key_result), col_indices_result);
    return result_sparse;
}

void SparseMatrixFixedPattern::calc_inverse_diagonal(std::vector<scalar> &inv_D) const {
//#ifdef FFEA_PARALLEL_WITHIN_BLOB
//#pragma omp parallel for default(none)  shared(inv_D)
//#endif
    for (int i = 0; i < num_rows; i++) {
        inv_D[i] = 1.0 / (*(diagonal[i]));
    }
}

void SparseMatrixFixedPattern::print() const {
    for (int i = 0; i < num_nonzero_elements; i++) {
        printf("[%d %e]\n", entry[i].column_index, entry[i].val);
    }
}

void SparseMatrixFixedPattern::print_dense() const {
    for (int i = 0; i < num_rows; ++i) {
        int l = 0;
        for (int j = 0; j < num_rows; ++j) {
            if (j == entry[key[i] + l].column_index && key[i] + l < key[i + 1]) {
                printf("%e,", entry[key[i] + l].val);
                l++;
            } else {
                printf("%e,", 0.0);
            }
        }
        printf("\n");
    }
}

/* Prints dense matrix out to file for analysis. I suggest only letting this function run once (step = 1?) */
void SparseMatrixFixedPattern::print_dense_to_file(std::vector<arr3> &a) const {
    FILE *fout, *fout2;
    fout = fopen("dense_matrix.csv", "w");
    fout2 = fopen("force.csv", "w");
    for (int i = 0; i < num_rows; ++i) {
        int l = 0;
        for (int j = 0; j < num_rows; ++j) {
            if (j == entry[key[i] + l].column_index && key[i] + l < key[i + 1]) {
                fprintf(fout, "%e,", entry[key[i] + l].val);
                l++;
            } else {
                fprintf(fout, "%e,", 0.0);
            }
        }
        fprintf(fout, "\n");
    }
    for (int i = 0; i < num_rows / 3; ++i) {
        fprintf(fout2, "%e\n%e\n%e\n", a[i][0], a[i][1], a[i][2]);
    }
    fclose(fout);
    fclose(fout2);
}

void SparseMatrixFixedPattern::print_row_column() const {
    FILE *fout;
    fout = fopen("row_column.csv", "w");
    for (int i = 0; i < num_rows; ++i) {
        for (int j = key[i]; j < key[i + 1]; ++j) {
            fprintf(fout, "%d,%d\n", entry[j].column_index, i);
        }
    }
    fclose(fout);
}

void SparseMatrixFixedPattern::check_symmetry() const {
    int row;
    for (int i = 0; i < num_rows / 3; i++) {
        for (int j = key[i]; j < key[i + 1]; j++) {
            printf("%e ", entry[j].val);
            row = entry[j].column_index;
            for (int k = key[row]; k < key[row + 1]; ++k) {
                if (entry[k].column_index == i) {
                    printf("%e %e\n", entry[k].val, entry[k].val - entry[j].val);
                }
            }
        }
    }
}

void SparseMatrixFixedPattern::am_i_diagonally_dominant() const {
    int i, j;
    scalar sum;
    FILE *fout;
    fout = fopen("diag_domin.csv", "w");
    for (i = 0; i < num_rows; ++i) {
        sum = 0.0;
        for (j = key[i]; j < key[i + 1]; ++j) {
            if (entry[j].column_index == i) {
                continue;
            }
            sum += fabs(entry[j].val);
        }
        fprintf(fout, "%d,%f\n", i, sum / fabs(*diagonal[i]));
    }
    fclose(fout);
}

std::vector<sparse_entry> &SparseMatrixFixedPattern::get_entries() {
    return entry;
}

std::vector<int> &SparseMatrixFixedPattern::get_key() {
    return key;
}

int SparseMatrixFixedPattern::get_num_nonzero_elements() const {
    return num_nonzero_elements;
}

int SparseMatrixFixedPattern::get_num_rows() const {
    return num_rows;
}

int SparseMatrixFixedPattern::get_num_columns() const {
    int i, num_columns = 0;
    for(i = 0; i < num_nonzero_elements; ++i) {
        if(entry[i].column_index > num_columns) {
            num_columns = entry[i].column_index;
        }
    }
    return num_columns + 1;
}
