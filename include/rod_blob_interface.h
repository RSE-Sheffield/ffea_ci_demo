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
 *      rod_blob_interface.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk
 */

#ifndef ROD_BLOB_INTERFACE
#define ROD_BLOB_INTERFACE

#include "rod_structure.h"

// this bit
//#include "SimulationParams.h"
#include "tetra_element_linear.h"
#include "Face.h"
#include "Blob.h"

namespace rod {

// functions go here

void get_tri_norm(const float3 &node0, const float3 &node1, const float3 &node2, OUT float3 &tri_norm);
void get_jacobian(const std::array<mesh_node*, NUM_NODES_LINEAR_TET> &tet_nodes, OUT float9 &J);
void float_3x3_invert(float9 &m, OUT float9 &m_inv);
void get_gradient_deformation(const float9 &J_inv_0, std::array<mesh_node*, NUM_NODES_LINEAR_TET> &nodes_curr, OUT float9 &transposed_gradient_deformation_3x3);
void QR_decompose_gram_schmidt(const float9 &matrix_3x3, OUT float9 &Q, float9 &R);
void construct_euler_rotation_matrix(float a, float b, float g, float9 &rotmat);
void rotate_tet(const float9 &rotmat, const std::array<mesh_node*, NUM_NODES_QUADRATIC_TET> &nodes, OUT std::array<mesh_node*, NUM_NODES_LINEAR_TET> &rotated_nodes);
void get_euler_angles(const float9 &rm, OUT float3 &euler_angles);
void get_rotation_matrix_from_euler(const float3 &euler_angles, OUT float9 &rm);
bool array_equal(const float3 &arr1, const float3 &arr2);
bool array_contains(const float4x3 &large_arr, const float3x3 &small_arr);
void mesh_node_null_check(mesh_node* node, std::string location);
void rescale_attachment_node(const float3 &attachment_node, const float3 &end_node, const float3 &attachment_node_equil, const float3 &end_node_equil, OUT float3 &scaled_attachment_node, float3 &scaled_attachment_node_equil);
void matrix_invert_3x3(const float9 &m, OUT float9 &inverted_matrix);
void equil_attachment_node_from_J(const float9 &J_inv_0, const int3 &face_node_indices, bool ends_at_rod, const float3 &node_weighting, const float3 &tet_origin, const float3x3 &edge_vecs, const float3 &rotation, OUT float3 &equil_attachment_node);
bool points_out_of_tet(const float3 &node1, const float3 &node2, const float3 &node3, const float3 &node4, const float3 &attachment_element, const float3 &attachment_node);
void get_attachment_node_pos(const float3 &face_node_1, const float3 &face_node_2, const float3 &face_node_3, const float3x3 &edge_vecs, const float3 &node_weighting, const float3 &tet_origin, OUT float3 &face_node_pos);

// objects go here

struct Rod_blob_interface
{
    // member variables
    Rod* connected_rod;
    Blob* connected_blob;
    Face* connected_face;
    tetra_element_linear* connected_tet;
    bool ends_at_rod = true; // if this is false, it goes rod->blob, otherwise it goes blob->rod
    int to_index;
    int from_index;
    int face_index;
    int order;
    int3 face_nodes;
    float3 node_weighting = {0.333333333333333, 0.333333333333333, 0.333333333333333};
    float3 euler_angles = {0, 0, 0};
    float3 tet_origin;
    //float tet_origin_equil[3];
    float3x3 edge_vecs;
    //float edge_vecs_equil[3][3];
    std::array<mesh_node*, 4> deformed_tet_nodes; // used for calculating the jacobian
    int3 face_node_indices;
    float9 J_inv_0; // equilibrium jacobian of the attachment node
    
    float3 attachment_node_equil;
    //float attachment_node_pos_equil[3]; // note: remove attachment_node_pos_equil, it is useless
    float3 attachment_m_equil;
    
    float3 attachment_node;
    float3 attachment_node_pos;
    float3 attachment_m;
    
    //methods
    Rod_blob_interface (Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int3 &face_nodes, float3 &rotation, float3 &node_weighting, int order);
    void set_initial_values();
    void update_internal_state(bool update_edge_vecs, bool update_tet);
    void update_J_0();
    void set_edge_vecs();
    void get_attachment_node(OUT float3 &attachment_node, float3 &attachment_node_pos, bool equil);
    void get_attachment_node_pos(float3 &attachment_node_pos, bool equil);
    void get_initial_material_axis(); // parallel transport nearest mataxis
    void get_attachment_material_axis(float3 &attachment_node, OUT float3 &attachment_material_axis);
    void get_tet_rotation_matrix(float tet_points_before[12], float tet_points_after[12]); // calls to get_gradient_deformation and qr_decompose
    
    void make_tet(); // turn nodes into tet elements
    void set_tet(tetra_element_linear *tet);
    void reorientate_connection(float3 &attachment_element_orig, float3 &attachment_material_axis_orig, OUT float3 &new_attachment_element, float3 &new_attachment_material_axis);
    void position_rod_from_blob(bool use_equil);
    void position_blob_from_rod();
    void get_face(OUT float3 &face_node_1, float3 &face_node_2, float3 &face_node_3);
    void select_face_nodes(OUT int3 &face_node_indices);
    int get_element_id(const int3 &nodes);
    void get_node_energy(int node_index, float3 &attachment_node_equil, float3 &attachment_material_axis_equil, float3 &attachment_node, float3 &attachment_material_axis, float displacement, float6 &energy);
    //void get_rod_energy(float3 &attachment_node_equil, float3 &attachment_material_axis_equil, float displacement, float energy[2][6]);
    void position_rod_ends(float3 &attachment_node_pos);
    void do_connection_timestep();
    

    // note: maybe a wrapper function for doubles that converts them to floats? see if it works
    
    // note: use vector12 class for internal nodes or something similar
    
};

bool point_inside_tetrahedron(float3 &point, float tet[12]); // work out what direction the attachment node should face (it should point inside the tet, i think?)

} //end namespace


#endif
