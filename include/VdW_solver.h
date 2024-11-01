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

#ifndef VDW_SOLVER_H_INCLUDED
#define VDW_SOLVER_H_INCLUDED

#include <cmath>

#include "FFEA_return_codes.h"
#include "NearestNeighbourLinkedListCube.h"
#include "LJ_matrix.h"
#include "Blob.h"

class VdW_solver {
public:
    /**
     * @brief Virtual destructor required for correct inheritance behaviour
     */
    virtual ~VdW_solver() = default;

    void init(NearestNeighbourLinkedListCube *surface_face_lookup, arr3 &box_size, SSINT_matrix *ssint_matrix, scalar &steric_factor, int num_blobs, int inc_self_ssint, string ssint_type_string, scalar &steric_dr, int calc_kinetics, bool working_w_static_blobs);

    void solve(scalar *blob_corr);

    /** Allow protein VdW interactions along the top and bottom x-z planes */
    void solve_sticky_wall(scalar h);

    scalar get_field_energy(int i, int j);

    void reset_fieldenergy(); 

protected:

    int total_num_surface_faces = 0;
    NearestNeighbourLinkedListCube *surface_face_lookup = nullptr;

    arr3 box_size = {};
    SSINT_matrix *ssint_matrix;

    std::vector<std::vector<scalar>> fieldenergy;
    int num_blobs = 0;
    int inc_self_ssint = 0;  ///< whether to include interactions between faces within the same blob, or not.
    int calc_kinetics = 0; 
    bool working_w_static_blobs = false;
    struct adjacent_cell_lookup_table_entry {
        int ix, iy, iz;
    };

    scalar steric_factor = 0; ///< Proportionality factor to the Steric repulsion.
    scalar steric_dr = 0; ///< Constant to calculate the numerical derivative.
    // static const scalar phi_f[4]; ///< shape function for the center of the "element"
    static const std::array<adjacent_cell_lookup_table_entry, 27> adjacent_cell_lookup_table;

    static const int num_tri_gauss_quad_points = 3; 
    struct tri_gauss_point {
        scalar W;
        arr3 eta;
    };
    // static const struct tri_gauss_point gauss_pointx[num_tri_gauss_quad_points];
    static const std::array<tri_gauss_point, 3> gauss_points;

    bool consider_interaction(Face *f1, int l_index_i, int motion_state_i, LinkedListNode<Face> *l_j, scalar *blob_corr);

    virtual void do_interaction(Face *f1, Face *f2, scalar *blob_corr);

    bool do_steric_interaction(Face *f1, Face *f2, scalar *blob_corr);

    void do_lj_interaction(Face *f1, Face *f2, scalar *blob_corr); 

    void do_gensoft_interaction(Face *f1, Face *f2, scalar *blob_corr);

    void do_sticky_xz_interaction(Face *f, bool bottom_wall, scalar dim_y);

    void calc_lj_force_pair_matrix(
              arr3 (&force_pair_matrix)[num_tri_gauss_quad_points][num_tri_gauss_quad_points],
              arr3 (&p)[num_tri_gauss_quad_points], arr3 (&q)[num_tri_gauss_quad_points], 
              scalar &Rmin, scalar &Emin, scalar &energy);

    void calc_ljinterpolated_force_pair_matrix(
              arr3 (&force_pair_matrix)[num_tri_gauss_quad_points][num_tri_gauss_quad_points],
              arr3 (&p)[num_tri_gauss_quad_points], arr3 (&q)[num_tri_gauss_quad_points], 
              scalar &Rmin, scalar &Emin, scalar &energy);

    void calc_gensoft_force_pair_matrix(arr3 (&force_pair_matrix)[num_tri_gauss_quad_points][num_tri_gauss_quad_points],
        arr3 (&p)[num_tri_gauss_quad_points], arr3 (&q)[num_tri_gauss_quad_points],
        scalar &Rmin, scalar &Emin, scalar &k0, scalar &energy);

    void calc_lj_factors(scalar &mag_r, int index_k, int index_l, scalar &Emin, scalar &Rmin_6,
                                 scalar &force_mag, scalar &e);

    void calc_ljinterpolated_factors(scalar &mag_r, int index_k, int index_l, scalar &Emin, scalar &Rmini,
                                 scalar &force_mag, scalar &e);

    void calc_gensoft_factors(scalar &mag_r, int index_k, int index_l, scalar &Emin, scalar &Rmin_2, scalar &Rmin_3, scalar &k0, 
                                 scalar &force_mag, scalar &e);

    scalar dot_with_normal(const arr3 &p, const arr3 &q, const arr3 &n);

    scalar minimum_image(scalar delta, scalar size);

    int ssint_type = SSINT_TYPE_UNDEFINED;
};

#endif
