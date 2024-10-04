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

/** 
 * \defgroup FMM Force Match
 * @{
 */

#include "PreComp_solver.h"

using std::cout;
using std::endl; 

const int PreComp_solver::adjacent_cells[27][3] = {
        {-1, -1, -1},
        {-1, -1, 0},
        {-1, -1, +1},
        {-1, 0, -1},
        {-1, 0, 0},
        {-1, 0, +1},
        {-1, +1, -1},
        {-1, +1, 0},
        {-1, +1, +1},
        {0, -1, -1},
        {0, -1, 0},
        {0, -1, +1},
        {0, 0, -1},
        {0, 0, 0},
        {0, 0, +1},
        {0, +1, -1},
        {0, +1, 0},
        {0, +1, +1},
        {+1, -1, -1},
        {+1, -1, 0},
        {+1, -1, +1},
        {+1, 0, -1},
        {+1, 0, 0},
        {+1, 0, +1},
        {+1, +1, -1},
        {+1, +1, 0},
        {+1, +1, +1}
};

PreComp_solver::PreComp_solver() { 
  nint = 0;
  msgc = 0;
  n_beads = 0;
  num_blobs = 0;
} 

/** @brief destructor: deallocates pointers. */
PreComp_solver::~PreComp_solver() {
    U.clear();
    F.clear();
    isPairActive.clear();
    b_types.clear();

    b_elems.clear();

    fieldenergy.clear();
    num_blobs = 0;
}


int PreComp_solver::msg(int whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
}
 

int PreComp_solver::msg(string whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
} 


/** Zero measurement stuff, AKA fieldenergy */
void PreComp_solver::reset_fieldenergy() {
#ifdef USE_OPENMP
#pragma omp parallel for default(none) 
#endif
    for (int i=0; i<num_threads; i++) {
        for(int j = 0; j < num_blobs; j++) {
           for(int k = 0; k < num_blobs; k++) {
              fieldenergy[(i*num_blobs + j)*num_blobs + k] = 0.0;
           }
        }
    }
}

void PreComp_solver::init(PreComp_params *pc_params, SimulationParams *params, Blob **blob_array) {

   /* Firstly, we get the number of lines, x_range, 
    * and Dx from the first pair type potential file. 
    * We check that the x_range and Dx 
    * are the same for all the files. Secondly, 
    * we allocate the F and U arrays, and then
    * we store the y_values, re-reading the files.
    * The functions get_U and get_F will then be ready. 
    * Thirdly and fourthly, allocate and fill the 
    * elements and relative bead position arrays.
    * And finally, delete beads stuff from the blobs,
    * and prepare a linkedlist.
    */  
  
   // do two simple checks: 
   //  1: number of types is greater than 0
   if (pc_params->types.size() == 0) {
       throw FFEAException("\n Number of bead types is Zero:\n\t correct it, or change to calc_PreComp = 0.");
   } 
   //  2: the folder exists:
   fs::path p = pc_params->folder;
   if (fs::exists(p)) {
     if (!fs::is_directory(p)) {
         throw FFEAException("\n Folder %s is not a folder.", pc_params->folder.c_str());
     }
   } else {
       throw FFEAException("\n Folder %s does not exist.", pc_params->folder.c_str());
   } 
    

#ifdef USE_OPENMP
    num_threads = omp_get_max_threads(); 
#else
    num_threads = 1;
#endif 
    num_blobs = params->num_blobs;
    try {
        fieldenergy = std::vector<scalar>(num_threads * num_blobs * num_blobs);
    } catch (std::bad_alloc &) {
       throw FFEAException("Failed to allocate memory for fieldenergy in PreCompSolver.");
    }

   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;
   scalar d2, d2_0;
   
   /*--------- FIRSTLY ----------*/ 
   /* get the number of interactions */
   for (int i=1; i<pc_params->types.size() + 1; i++){
     nint += i;
   } 
   msg(nint);
   // and the number of bead types: 
   ntypes = pc_params->types.size();

   /* find a first file that exists */
   unsigned int eti, etj;
   bool findSomeFile = false; 
   for (eti=0; eti<pc_params->types.size(); eti++) {
     for (etj=eti; etj<pc_params->types.size(); etj++) { 
       // open file i-j
       ssfile << pc_params->folder << "/" << pc_params->types[eti] << "-" << pc_params->types[etj] << ".pot";
       fin.open(ssfile.str(), std::ifstream::in);
       if (fin.is_open()) {
         fin.close();
         findSomeFile = true;
         goto end_loop;
       }
       ssfile.str("");
     }
   }
   end_loop:
   // and if no file was found, protest!
   if (findSomeFile == false) {
        throw FFEAException("Failed to open any file potential file in folder: %s.", pc_params->folder.c_str());
   } 
   /* read the first file, and get the number of lines
    * (n_values), Dx, and x_range */
   // the first file:
   ssfile.str("");
   ssfile << pc_params->folder << "/" << pc_params->types[eti] << "-" << pc_params->types[etj] << ".pot";
   fin.open(ssfile.str());
   if (fin.fail()) {
        throw FFEAFileException(ssfile.str().c_str());
   }
   // get the first line that does not start with "#"
   getline(fin, line);
   while (line.find("#", 0, 1) == 0) {
     getline(fin, line); 
   } 
   // get x_0, i. e., parse the line:
   boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
   x_0 = stod(vec_line[0]); 
   // and store it: 
   x_range[0] = x_0;
   getline(fin, line);
   boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
   x = stod(vec_line[0]); 
   Dx = x - x_0;
   n_values = 2;
   x_0 = x;
   while ( getline(fin, line) ) {
      n_values += 1;
      boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
      x = stod(vec_line[0]); 
      if (fabs(Dx - x + x_0) > 1e-6) { 
          throw FFEAException("Aborting: delta x was found to be non-uniform in: %s", ssfile.str().c_str());
      } 
      x_0 = x; 
   }  
   fin.close();

   /*--------- SECONDLY ----------*/ 
   // allocate:
   try {
       U = std::vector<scalar>(n_values * nint);
       F = std::vector<scalar>(n_values * nint);
       isPairActive = std::vector<bool>(n_values * nint);
   } catch (std::bad_alloc &) {
       throw FFEAException("Failed to allocate memory for arrays in PreComp_solver::init.");
   }
      
   // and load potentials and forces:
   read_tabulated_values(*pc_params, "pot", U, pc_params->E_to_J / mesoDimensions::Energy);
   if (pc_params->inputData == 1) {
     // scalar F_to_Jm = pc_params->E_to_J / pc_params->dist_to_m;
     scalar F_scale = ( pc_params->E_to_J / pc_params->dist_to_m ) / mesoDimensions::force;
     read_tabulated_values(*pc_params, "force", F, F_scale);
   } else if (pc_params->inputData == 2) {
        calc_force_from_pot();
   } else {
        throw FFEAException("--- PreComp: invalid value for precomp->inputData");
   }
   Dx = Dx * pc_params->dist_to_m / mesoDimensions::length ;
   x_range[0] = x_range[0] * pc_params->dist_to_m / mesoDimensions::length; 
   x_range[1] = x_range[0] + Dx * n_values;
   x_range2[0] = x_range[0]*x_range[0]; 
   x_range2[1] = x_range[1]*x_range[1]; 



   /*------------ THIRDLY --------*/
   // allocate the list of elements that have a bead will use: 
    for (int i=0; i < params->num_blobs; i ++) {
        // we only consider one conformation per blob.
        if (params->num_conformations[i] > 1) {
            throw FFEAException("--- PreComp: currently PreComp only deals with a single conformation per blob");
        }
     n_beads += blob_array[i][0].get_num_beads();
    } 
   // Check that we have some beads!:
   if (n_beads == 0) {
      throw FFEAException("--- PreComp: The total number of beads is 0, but PreComp_calc was set to 1.");
   }
   try {
       b_elems = std::vector<TELPtr>(n_beads);
       // allocate the array that store the relative positions 
       //    of the beads to the elements where they belong to. 
       b_rel_pos = std::vector<scalar>(n_beads * 3);
       // allocate the array to compute the absolute positions of the beads:
       b_pos = std::vector<scalar>(n_beads * 3);
       // and allocate the bead types: 
       b_types = std::vector<int>(n_beads);
   } catch (std::bad_alloc &) {
       throw FFEAException("Failed to allocate memory for array beads in PreComp_solver::init");
   }
   

   /*------------ FOURTHLY --------*/
   // get the elements of the list of "elements" that we will use: 
   arr3 u, v, w; 
   tetra_element_linear *e;
   matrix3 J, J_inv; // will hold the Jacobian for the current element. 
   scalar det;  // determinant for J.
   int m = 0;
   int n;
   num_diff_elems = 0; // number of different elements. 
   vector<vector<int>> b_elems_map_tmp; 
   // for each Blob: 
   for (int i=0; i < params->num_blobs; i ++) {
     int num_diff_elems_in_blob = 0; 
     vector<vector<int>> sorting_pattern; 
     // store the bead types: 
     n = blob_array[i][0].get_num_beads();
     memcpy(&b_types[m], blob_array[i][0].get_bead_types().data(), n * sizeof(int));

     // for each bead within this blob (remember that we only deal with conf 0):
     for (int j=0; j < n; j++) {
       blob_array[i][0].get_bead_position(j, v);
       vector<int> b_assignment = blob_array[i][0].get_bead_assignment(j); 
       d2_0 = 1e9;
       int mj = m+j;
       // get the closest node to this bead: 
       for (int k=0; k < blob_array[i][0].get_num_elements(); k++) { 
         e = blob_array[i][0].get_element(k);
         bool work = true; 
         // check that this element has one of the chosen nodes: 
         if (!b_assignment.empty()) {
           work = false; 
           for (int l=0; l < NUM_NODES_QUADRATIC_TET; l++) {
             if (find(b_assignment.begin(), b_assignment.end(), e->n[l]->index) != b_assignment.end()) {
               work = true;
               break;
             }
           }
         } 
         if (work == false) continue;
         e->calc_centroid();
         d2 = (e->centroid[0] - v[0])*(e->centroid[0] - v[0]) + 
              (e->centroid[1] - v[1])*(e->centroid[1] - v[1]) + 
              (e->centroid[2] - v[2])*(e->centroid[2] - v[2]);
       
         if (d2 < d2_0) {
           d2_0 = d2;
           b_elems[mj] = e; 
         }
       } 
 
       // did we get a different element? 
       bool newelem = true; 
       int mj_elem_index = b_elems[mj]->index;
       for (int k = m; k<m+j; k++){
          if (mj_elem_index == b_elems[k]->index) {
            newelem = false; 
            break; 
          }
       } 
       if (newelem == true) num_diff_elems_in_blob += 1;
       // record the elem stamp to sort arrays later. 
       vector<int> elem_stamp(2);
       elem_stamp = {mj_elem_index, mj}; 
       sorting_pattern.push_back(elem_stamp);
       
       
       // and get the relative coordinates within the element
       //   as a fraction of the basis vectors length.
       b_elems[mj]->calculate_jacobian(J); 
       mat3_invert(J, J_inv, &det);
       sub(v, b_elems[mj]->n[0]->pos, w);
       vec3_mat3_mult(w, J_inv, u); 
       // now u has the relative coordinates, not under unit vectors
       //    but under full length vectors. And we store them:
       b_rel_pos[3*mj] = u[0];
       b_rel_pos[3*mj+1] = u[1];
       b_rel_pos[3*mj+2] = u[2];
       // cout << "0ead " << mj << " in: " << v[0] << ", " << v[1] << ", " << v[2] << endl;

       // ESSENTIAL printout to relate beads to nodes!! 
       //   in case of being required.
       scalar l = mesoDimensions::length;
       stringstream beadsToNodes;
       beadsToNodes << "bead type: " << b_types[m+j] <<  ", position " <<
                        v[0]*l << ":" <<  v[1]*l << ":" << v[2]*l << " in element " <<
                        b_elems[m+j]->n[0]->index << ":" << b_elems[m+j]->n[1]->index
                    <<  b_elems[m+j]->n[2]->index << ":" << b_elems[m+j]->n[3]->index
                    <<  " position " << b_elems[m+j]->centroid[0]*l << ":" <<
                        b_elems[m+j]->centroid[1]*l << ":" << b_elems[m+j]->centroid[2]*l
                    <<  " - distance: " << sqrt(d2_0)*l << ", 2ndOE: " <<
                        b_elems[m+j]->n[4]->index << ":" << b_elems[m+j]->n[5]->index 
                    << ":" << b_elems[m+j]->n[6]->index << ":" <<
                        b_elems[m+j]->n[7]->index << ":" << b_elems[m+j]->n[8]->index
                    << ":" << b_elems[m+j]->n[9]->index;
       print_high(beadsToNodes.str());
       /*
       //  prove it: v =? s
       arr3 s, e1, e2, e3;
       vec3_vec3_subs(&b_elems[mj]->n[1]->pos, &b_elems[m+j]->n[0]->pos, &e1);
       vec3_vec3_subs(&b_elems[mj]->n[2]->pos, &b_elems[m+j]->n[0]->pos, &e2);
       vec3_vec3_subs(&b_elems[mj]->n[3]->pos, &b_elems[m+j]->n[0]->pos, &e3);
       s.x = b_elems[mj]->n[0]->pos.x + u.x*e1.x + u.y*e2.x + u.z*e3.x;
       s.y = b_elems[mj]->n[0]->pos.y + u.x*e1.y + u.y*e2.y + u.z*e3.y;
       s.z = b_elems[mj]->n[0]->pos.z + u.x*e1.z + u.y*e2.z + u.z*e3.z;
       print_arr3(v);
       print_arr3(s);
       s.x = b_elems[mj]->n[0]->pos.x + u.x*J[0][0] + u.y*J[1][0] + u.z*J[2][0];
       s.y = b_elems[mj]->n[0]->pos.y + u.x*J[0][1] + u.y*J[1][1] + u.z*J[2][1];
       s.z = b_elems[mj]->n[0]->pos.z + u.x*J[0][2] + u.y*J[1][2] + u.z*J[2][2];
       print_arr3(s);
       */

       /* the following is useless but useful while testing:
       b_pos[3*mj] = v.x;
       b_pos[3*mj+1] = v.y;
       b_pos[3*mj+2] = v.z;
       */
     } 
     num_diff_elems += num_diff_elems_in_blob; 

     //  Now sort arrays: 
     std::sort(sorting_pattern.begin(), sorting_pattern.end(), [](vector<int> i, vector<int> j){ return i[0]<j[0]; });
     for (int j=0; j<n; j++) {
       int mj = m+j;
       // tetra_element_linear *e; // already defined 
       int tmp_type; 
       arr3 tmp_pos; 

       int swap_high = sorting_pattern[j][1];
       int swap_low = mj;

       e = b_elems[swap_low];
       b_elems[swap_low] = b_elems[swap_high]; 
       b_elems[swap_high] = e; 

       tmp_type = b_types[swap_low]; 
       b_types[swap_low] = b_types[swap_high];
       b_types[swap_high] = tmp_type;

       store( arr_view<scalar,3>(b_rel_pos, 3*swap_low), tmp_pos); 
       store( arr_view<scalar,3>(b_rel_pos, 3*swap_high), arr_view<scalar,3>(b_rel_pos, 3*swap_low) ); 
       store( tmp_pos, arr_view<scalar,3>(b_rel_pos, 3*swap_high)); 

       for (int k=j; k<n; k++) {
         if (sorting_pattern[k][1] == swap_low) {
           sorting_pattern[k][1] = swap_high; 
           break;
         } 
       } 
       sorting_pattern[j][1] = swap_low;
     }
     m += n;
   }

   // now set up the map and create a list of unique elements.
   try {
       b_forces = std::vector<scalar>(3 * n_beads);
       map_e_to_b = std::vector<int>(2 * num_diff_elems);
       b_unq_elems = std::vector<TELPtr>(num_diff_elems);
   } catch (std::bad_alloc &) {
       throw FFEAException("Failed to allocate memory for supplementary array beads in PreComp_solver::init.");
   }
   m = 0; 
   int cnt = 0;
   for (int i=0; i < params->num_blobs; i ++) {
     n = blob_array[i][0].get_num_beads();
     int updating = 0; 
     int prev = b_elems[m]->index; 
     for (int  j=0; j<n; j++) {
       int mj = m+j; 
       if (updating == 1) {
         if (b_elems[mj]->index != prev) {
           updating = 0; 
           map_e_to_b[2*cnt + 1] = mj -1; 
           cnt += 1;
         }
       } 
       if (updating == 0) {
           b_unq_elems[cnt] = b_elems[mj]; 
           map_e_to_b[2*cnt] = mj; 
           updating = 1;
       } 
       prev = b_elems[mj]->index;
     } 
     m += n; 
     map_e_to_b[2*cnt +1] = m - 1;
     cnt += 1;

     // and make Blob[i] forget all about beads. 
     blob_array[i][0].forget_beads();
   } 

   // even if we have the elems, pointing to blob_index, 
   //   it may better to just store the indices, 
   //   so that b_elems can be removed some day (and used only for debugging purposes). 
   try {
       b_daddyblob = std::vector<int>(n_beads); 
   } catch (std::bad_alloc &) {
       throw FFEAException("Failed to allocate memory for daddy beads array in PreComp_solver::init.");
   }
   for (int i=0; i<n_beads; i++) {
      b_daddyblob[i] = b_elems[i]->daddy_blob->blob_index;
   } 


   num_blobs = params->num_blobs;


   /*------------ FIFTHLY ---:)---*/
   // Set up the linkedlist:
   // 5.1 - Allocate the linkedlist:
   printf("Allocating memory for the preComp linked list grid...\n");
   // the total size of the box is: 
   scalar vdwVoxelSize = params->ssint_cutoff;
   scalar dimBox[3];
   dimBox[0] = params->es_N_x * vdwVoxelSize;
   dimBox[1] = params->es_N_y * vdwVoxelSize;
   dimBox[2] = params->es_N_z * vdwVoxelSize;
   // and our cell size is x_range[1], so we'll adjust the size of the voxel
   //   to fit an integer number of cells in each direction:
   if (x_range[1] > vdwVoxelSize) {
     pcVoxelSize = ceil(x_range[1] / vdwVoxelSize) * vdwVoxelSize; 
   } else { 
     pcVoxelSize = vdwVoxelSize / floor(vdwVoxelSize/x_range[1]) ;
   } 
   for (int i=0; i<3; i++) {
     pcVoxelsInBox[i] = dimBox[i] / pcVoxelSize; 
     if (pcVoxelsInBox[i] == 0) {
         throw FFEAException("Zero voxels were allocated in the %d th direction.", i);
     }
   }
   printf("... with %d x %d x %d voxels of size %f, while vdwVoxelSize was %f and pc_range %f\n",
                   pcVoxelsInBox[0], pcVoxelsInBox[1], pcVoxelsInBox[2], pcVoxelSize, vdwVoxelSize, x_range[1]); 

#ifdef FFEA_PARALLEL_FUTURE
   pcLookUp.alloc_dual(pcVoxelsInBox[0], pcVoxelsInBox[1], pcVoxelsInBox[2], n_beads);
#else
   pcLookUp.alloc(pcVoxelsInBox[0], pcVoxelsInBox[1], pcVoxelsInBox[2], n_beads);
#endif


   // 5.2 - store indices in there:
   for (int i=0; i<n_beads; i++) {
#ifdef FFEA_PARALLEL_FUTURE
     pcLookUp.add_to_pool_dual(nullptr);
#else
     pcLookUp.add_to_pool(nullptr);
#endif
   } 
   // one last check: see if all beads were stored.
   if (pcLookUp.get_pool_size() != n_beads) {
      throw FFEAException(" The number of beads in the LinkedList %d is not the same as the total number of beads %d.", pcLookUp.get_pool_size(), n_beads); 
   }

   // 5.3 - We cannot compute_bead_positions here because the system is not into the box yet.
   //        So 5 is done. 


   /*------------ SIXTHLY ---:O---*/
   //  in case we're writing the trajectory for the beads we need to keep some data:
   if (params->trajbeads_fname_set == 1) {
      stypes = pc_params->types;
       try {
           b_blob_ndx = std::vector<int>(n_beads);
           b_elems_ndx = std::vector<int>(n_beads);
       } catch (std::bad_alloc &) {
           throw FFEAException("Failed to allocate memory for beads details to write trajectory in PreComp_solver::init.");
       }
      for (int i=0; i<n_beads; i++){ 
        b_elems_ndx[i] = b_elems[i]->index;
        b_blob_ndx[i] = b_elems[i]->daddy_blob->blob_index;
      } 
   } 


   cout << "done!" << endl;
}

void PreComp_solver::solve_using_neighbours_non_critical(scalar *blob_corr/*=nullptr*/){
    scalar d, f_ij; //, f_ijk_i, f_ijk_j; 
    arr3 dx, dtemp, dxik;
    int type_i; 
    scalar phi_i[4]; 
    tetra_element_linear *e_i, *e_j;

    // 0 - clear fieldenery:
    reset_fieldenergy(); 
    //   - and reset b_forces!
    for (int i=0; i<3*n_beads; i++) {
      b_forces[i] = 0.;
    } 


    // 1 - Compute the position of the beads:
    compute_bead_positions();


    // 2 - Compute all the i-j forces:
    LinkedListNode<int> *b_i = nullptr; 
    LinkedListNode<int> *b_j = nullptr; 
    int b_index_i, b_index_j; 
    int daddy_i, daddy_j;
#ifdef USE_OPENMP
#pragma omp parallel default(none) shared(blob_corr) private(type_i,phi_i,e_i,e_j,dx,dxik,d,dtemp,f_ij,b_i,b_j,b_index_i,b_index_j,daddy_i, daddy_j)
    {
    int thread_id = omp_get_thread_num(); 
    #pragma omp for
#else
    int thread_id = 0;
#endif
    for (int i=0; i<n_beads; i++){
      b_i = pcLookUp.get_from_pool(i); 
      b_index_i = b_i->index; 

      type_i = b_types[b_index_i]; 
      e_i = b_elems[b_index_i];  
      daddy_i = b_daddyblob[b_index_i];

      for (int c=0; c<27; c++) {
        b_j = pcLookUp.get_top_of_stack(b_i->x + adjacent_cells[c][0], 
                                        b_i->y + adjacent_cells[c][1],
                                        b_i->z + adjacent_cells[c][2]);
  
        while (b_j != nullptr) {
           b_index_j = b_j->index;
           if (b_index_j == b_index_i) {
             b_j = b_j->next; 
             continue; 
           } 

           if (!isPairActive[type_i*ntypes+b_types[b_index_j]]) {
             b_j = b_j->next; 
             continue; 
           }

           daddy_j = b_daddyblob[b_index_j];
           if (!blob_corr) {
             dx[0] = (b_pos[3*b_index_j  ] - b_pos[3*b_index_i  ]);
             dx[1] = (b_pos[3*b_index_j+1] - b_pos[3*b_index_i+1]);
             dx[2] = (b_pos[3*b_index_j+2] - b_pos[3*b_index_i+2]);
           } else {
             dx[0] = (b_pos[3*b_index_j  ] - blob_corr[3*(daddy_i*num_blobs + daddy_j)   ] - b_pos[3*b_index_i  ]);
             dx[1] = (b_pos[3*b_index_j+1] - blob_corr[3*(daddy_i*num_blobs + daddy_j) +1] - b_pos[3*b_index_i+1]);
             dx[2] = (b_pos[3*b_index_j+2] - blob_corr[3*(daddy_i*num_blobs + daddy_j) +2] - b_pos[3*b_index_i+2]);
           }
           d = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; 
           if (d > x_range2[1]) {
             b_j = b_j->next; 
             continue; 
           }
           else if (d < x_range2[0]) {
             b_j = b_j->next; 
             continue; 
           }
           d = sqrt(d);
           dx[0] = dx[0] / d;
           dx[1] = dx[1] / d;
           dx[2] = dx[2] / d;
 
           f_ij = get_F(d, type_i, b_types[b_index_j]); 
           // += the force: 
           resize3(f_ij, dx, arr_view<scalar,3>(b_forces, 3*b_index_i) );

           // e_j = b_elems[b_index_j];
           // fieldenergy[(thread_id*num_blobs + e_i->daddy_blob->blob_index) * num_blobs + e_j->daddy_blob->blob_index] += 0.5*get_U(d, type_i, b_types[b_index_j]);
           fieldenergy[(thread_id*num_blobs + daddy_i) * num_blobs + daddy_j] += 0.5*get_U(d, type_i, b_types[b_index_j]);

           b_j = b_j->next; 
        } // close b_j, beads in neighbour voxel loop 
      } // close c, 27 voxels loop 
    }  // close i, n_beads loop


    #pragma omp for
    for (int i=0; i<num_diff_elems; i++) {
      e_i = b_unq_elems[i]; 
      for (int j=map_e_to_b[2*i]; j<=map_e_to_b[2*i+1]; j++) {
        b_index_i = j; 

        phi_i[1] = b_rel_pos[3*b_index_i  ];
        phi_i[2] = b_rel_pos[3*b_index_i+1];
        phi_i[3] = b_rel_pos[3*b_index_i+2];
        phi_i[0] = 1 - phi_i[1] - phi_i[2] - phi_i[3];

        // fix input force:
        for (int k=0; k<4; k++) {
          resize2(-phi_i[k], arr_view<scalar,3>(b_forces, 3*b_index_i), dxik); 
          e_i->add_force_to_node(k, dxik);
        } // close k, nodes for the elements.
      }
    } 
#ifdef USE_OPENMP
    }
#endif
}


void PreComp_solver::solve_using_neighbours() {
    scalar d, f_ij; //, f_ijk_i, f_ijk_j; 
    arr3 dx, dtemp, dxik, dxjk;
    int type_i; 
    scalar phi_i[4], phi_j[4];
    tetra_element_linear *e_i, *e_j;

    // 0 - clear fieldenery:
    reset_fieldenergy(); 


    // 1 - Compute the position of the beads:
    compute_bead_positions();


    // 2 - Compute all the i-j forces:
    LinkedListNode<int> *b_i = nullptr; 
    LinkedListNode<int> *b_j = nullptr; 
    int b_index_i, b_index_j; 
#ifdef USE_OPENMP
#pragma omp parallel default(none) private(type_i,phi_i,phi_j,e_i,e_j,dx,d,dtemp,f_ij,b_i,b_j,b_index_i,b_index_j,dxik,dxjk)
    {
    int thread_id = omp_get_thread_num(); 
    #pragma omp for
#else
    int thread_id = 0;
#endif
    for (int i=0; i<n_beads; i++){
      b_i = pcLookUp.get_from_pool(i); 
      b_index_i = b_i->index; 

      type_i = b_types[b_index_i]; 
      phi_i[1] = b_rel_pos[3*b_index_i  ];
      phi_i[2] = b_rel_pos[3*b_index_i+1];
      phi_i[3] = b_rel_pos[3*b_index_i+2];
      phi_i[0] = 1 - phi_i[1] - phi_i[2] - phi_i[3];
      e_i = b_elems[b_index_i];  

      for (int c=0; c<27; c++) {
        b_j = pcLookUp.get_top_of_stack(b_i->x + adjacent_cells[c][0], 
                                        b_i->y + adjacent_cells[c][1],
                                        b_i->z + adjacent_cells[c][2]);
  
        while (b_j != nullptr) {
           b_index_j = b_j->index;
           if (b_index_j <= b_index_i) {
             b_j = b_j->next; 
             continue; 
           } 

           if (!isPairActive[type_i*ntypes+b_types[b_index_j]]) {
             b_j = b_j->next; 
             continue; 
           }
           dx[0] = (b_pos[3*b_index_j  ] - b_pos[3*b_index_i  ]);
           dx[1] = (b_pos[3*b_index_j+1] - b_pos[3*b_index_i+1]);
           dx[2] = (b_pos[3*b_index_j+2] - b_pos[3*b_index_i+2]);
           d = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; 
           if (d > x_range2[1]) {
             b_j = b_j->next; 
             continue; 
           }
           else if (d < x_range2[0]) {
             b_j = b_j->next; 
             continue; 
           }
           d = sqrt(d);
           dx[0] = dx[0] / d;
           dx[1] = dx[1] / d;
           dx[2] = dx[2] / d;
 
           f_ij = get_F(d, type_i, b_types[b_index_j]); 

           // Add energies to record 
           e_j = b_elems[b_index_j];
           fieldenergy[(thread_id*num_blobs + e_i->daddy_blob->blob_index) * num_blobs + e_j->daddy_blob->blob_index] += get_U(d, type_i, b_types[b_index_j]);

           resize(f_ij, dx);

           phi_j[1] = b_rel_pos[3*b_index_j];
           phi_j[2] = b_rel_pos[3*b_index_j+1];
           phi_j[3]= b_rel_pos[3*b_index_j+2];
           phi_j[0]= 1 - phi_j[1] - phi_j[2] - phi_j[3];

           // and apply the force to all the nodes in the elements i and j:
           #pragma omp critical
           {
           for (int k=0; k<4; k++) {
               resize2(-phi_i[k], dx, dxik);
               resize2(phi_j[k], dx, dxjk);
               e_i->add_force_to_node(k, dxik);
               e_j->add_force_to_node(k, dxjk); 
           } // close k, nodes for the elements.
           } // close critical


           b_j = b_j->next; 
        } // close b_j, beads in neighbour voxel loop 
      } // close c, 27 voxels loop 
    }  // close i, n_beads loop
#ifdef USE_OPENMP
    }
#endif
}

void PreComp_solver::solve(scalar *blob_corr/*=nullptr*/) {
    scalar d, f_ij; //, f_ijk_i, f_ijk_j; 
    arr3 dx, dtemp;
    int type_i; 
    scalar phi_i[4], phi_j[4];
    tetra_element_linear *e_i, *e_j;
    int daddy_i, daddy_j;

    // 0 - clear fieldenery:
    reset_fieldenergy(); 

    // 1 - Compute the position of the beads:
    compute_bead_positions();

    // 2 - Compute all the i-j forces:
    /*scalar e_tot = 0.0; 
    scalar f_tot = 0.0;*/
#ifdef USE_OPENMP
#pragma omp parallel default(none) shared(blob_corr) private(type_i,phi_i,phi_j,e_i,e_j,dx,d,dtemp,f_ij,daddy_i,daddy_j)
    {
    int thread_id = omp_get_thread_num(); 
    #pragma omp for
#else
    int thread_id = 0;
#endif
    for (int i=0; i<n_beads; i++){
      type_i = b_types[i]; 
      phi_i[1] = b_rel_pos[3*i];
      phi_i[2] = b_rel_pos[3*i+1];
      phi_i[3] = b_rel_pos[3*i+2];
      phi_i[0] = 1 - phi_i[1] - phi_i[2] - phi_i[3];
      e_i = b_elems[i];  
      daddy_i = b_daddyblob[i];
      for (int j=i+1; j<n_beads; j++) {
        if (!isPairActive[type_i*ntypes+b_types[j]]) continue; 

        daddy_j = b_daddyblob[j];
        if (!blob_corr) {
            dx[0] = (b_pos[3*j  ] - b_pos[3*i]);
            dx[1] = (b_pos[3*j+1] - b_pos[3*i+1]);
            dx[2] = (b_pos[3*j+2] - b_pos[3*i+2]);
        } else {
            dx[0] = (b_pos[3*j  ] - blob_corr[3*(daddy_i*num_blobs + daddy_j)   ] - b_pos[3*i  ]);
            dx[1] = (b_pos[3*j+1] - blob_corr[3*(daddy_i*num_blobs + daddy_j) +1] - b_pos[3*i+1]);
            dx[2] = (b_pos[3*j+2] - blob_corr[3*(daddy_i*num_blobs + daddy_j) +2] - b_pos[3*i+2]);
        }

        d = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; 
        if (d > x_range2[1]) continue;
        else if (d < x_range2[0]) continue;
        d = sqrt(d);
        dx[0] = dx[0] / d;
        dx[1] = dx[1] / d;
        dx[2] = dx[2] / d;
        
 
        f_ij = get_F(d, type_i, b_types[j]); 

        /*e_tot += get_U(d, type_i, b_types[j]);
        f_tot += f_ij;*/
        /*cout << "i: " << i << " j: " << j << " type_i: " << type_i << " type_j: " << b_types[j]
                      << " i.pos: " << mesoDimensions::length*b_pos[3*i]*1e9 << ":" << mesoDimensions::length*b_pos[3*i+1]*1e9 << ":" << mesoDimensions::length*b_pos[3*i+2]*1e9
                      << " j.pos: " << mesoDimensions::length*b_pos[3*j]*1e9 << ":" << mesoDimensions::length*b_pos[3*j+1]*1e9 << ":" << mesoDimensions::length*b_pos[3*j+2]*1e9
                      << " d: " << d*mesoDimensions::length*1e9 
                      << " U: " << mesoDimensions::Energy*get_U(d, type_i, b_types[j])/0.1660539040e-20 
                      << " F: " << mesoDimensions::force*f_ij/0.1660539040e-11 << endl;*/
        e_j = b_elems[j];

	// Add energies to record 
   fieldenergy[(thread_id*num_blobs + e_i->daddy_blob->blob_index) * num_blobs + e_j->daddy_blob->blob_index] += get_U(d, type_i, b_types[j]);

        resize(f_ij, dx);
        store(dx, dtemp); 

        phi_j[1] = b_rel_pos[3*j];
        phi_j[2] = b_rel_pos[3*j+1];
        phi_j[3]= b_rel_pos[3*j+2];
        phi_j[0]= 1 - phi_j[1] - phi_j[2] - phi_j[3];
        #pragma omp critical
        {
        // and apply the force to all the nodes in the elements i and j:
        for (int k=0; k<4; k++) {
          // forces for e_i
          resize(-phi_i[k], dx);
          e_i->add_force_to_node(k, dx);
          store(dtemp, dx); 
          // forces for e_j
          resize(phi_j[k], dx);
          e_j->add_force_to_node(k, dx);
          store(dtemp, dx); 

        } 
        }
      }
    }
#ifdef USE_OPENMP
    }
#endif
    // cout << " total energy: " << e_tot*mesoDimensions::Energy/0.1660539040e-20 << endl;
}


void PreComp_solver::compute_bead_positions() {
    matrix3 J; 
#ifdef USE_OPENMP
#pragma omp parallel for default(none) private(J)
#endif
    for (int i=0; i<n_beads; i++){
       b_elems[i]->calculate_jacobian(J); 
       #pragma omp simd
       for (int j=0; j<3; j++) {
          b_pos[3*i+j] = b_elems[i]->n[0]->pos[j] + b_rel_pos[3*i]*J[0][j] + b_rel_pos[3*i+1]*J[1][j] + b_rel_pos[3*i+2]*J[2][j];
       }
       // cout << "bead " << i << " in: " << b_pos[3*i]*mesoDimensions::length << ", " << b_pos[3*i+1]*mesoDimensions::length << ", " << b_pos[3*i+2]*mesoDimensions::length << endl;
 
    }
}


/** Calculate F as the numerical derivative of U, instead of reading .force files 
  * @param[in] int total_int: total number of different interactions.
  */
void PreComp_solver::calc_force_from_pot() {
   // scalar ym, y0, yM;
   scalar twoDx = 2*Dx;

   int index = -1;
   for (int i=0; i<nint; i++) {
     index += 1; 
     F[index] = (U[index] - U[index + 1]) / Dx;
     for (int j=1; j<n_values-1; j++) {
       index += 1;
       F[index] = (U[index -1] - U[index + 1]) / twoDx; 
     }
     index += 1;
     F[index] = 0; 
   }
}

/** Read either the .pot or .force files
  * and store its contents either in U or F */
void PreComp_solver::read_tabulated_values(PreComp_params &pc_params, string kind, std::vector<scalar> &Z, scalar scale_Z){
   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;

   int m_values;
   int index = 0;
   // read the potentials/forces in while checking:
   for (unsigned int i=0; i<pc_params.types.size(); i++) {
     for (unsigned int j=i; j<pc_params.types.size(); j++) { 
       // open file i-j
       ssfile << pc_params.folder << "/" << pc_params.types[i] << "-" << pc_params.types[j] << "." << kind;
       fin.open(ssfile.str(), std::ifstream::in);
       // if the file were not found, fill the table with zeroes, and mark the pair as inactive:
       if (!fin.is_open()) {
         FFEA_CAUTION_MESSG(" failed to open %s, so filling with zeroes interaction %s:%s\n", ssfile.str().c_str(), pc_params.types[i].c_str(), pc_params.types[j].c_str());
         for (int k=index; k<index+n_values; k++) { 
           Z[k] = 0.0;
         }
         index += n_values; 
         ssfile.str("");
         isPairActive[i*ntypes+j] = false;
         isPairActive[j*ntypes+i] = false; 
         continue; 
       }
       // get the first line that does not start with "#"
       getline(fin, line);
       while (line.find("#", 0, 1) == 0) {
         getline(fin, line); 
       }  
       // get x_0, i. e., parse the line:
       boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
       x_0 = stod(vec_line[0]); 
       // check it: 
       if (x_range[0] != x_0){
           throw FFEAException("different starting point for file: %s", ssfile.str().c_str());
        }
       // and store the potential:
       Z[index] = stod(vec_line[1]) * scale_Z;
       index += 1;
       // get the next line, and check Dx:
       getline(fin, line);
       boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
       x = stod(vec_line[0]); 
       if (fabs(Dx - x + x_0) > 1e-6){
           throw FFEAException("different step size for x in file: %s", ssfile.str().c_str());
       }
       // and store the potential:
       Z[index] = stod(vec_line[1]) * scale_Z;
       index += 1;
       // Now, go for the rest of the file:
       m_values = 2;
       x_0 = x;
       while ( getline(fin, line) ) {
          m_values += 1;
          boost::split( vec_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
          x = stod(vec_line[0]); 
          // check Dx at every line:
          if (fabs(Dx - x + x_0) > 1e-6) {
              throw FFEAException("delta x was found to be non-uniform for file: %s", ssfile.str().c_str());
          } 
          x_0 = x; 
          // and check that the file is not too long:
          if (m_values > n_values) { 
            FFEA_error_text();
            throw FFEAException("too many points for file: %s", ssfile.str().c_str());
          } 
          Z[index] = stod(vec_line[1]) * scale_Z;
          index += 1;
       }  
       // and check that it is not too short:
       if (m_values != n_values) {
           throw FFEAException("too few points for file: %s", ssfile.str().c_str());
       } 
       fin.close();
       ssfile.str("");
       isPairActive[i*ntypes+j] = true;
       isPairActive[j*ntypes+i] = true; 
     }
   }
}

/**
 * @brief Use finterpolate to return the value for U at x between types i and j.
 * 
 * @details Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_U(scalar x, int typei, int typej) { 
  return finterpolate(U, x, typei, typej);
}


/**
 * @brief Use finterpolate to return the value for F at x between types i and j.
 *
 * @details Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_F(scalar x, int typei, int typej) { 
  return finterpolate(F, x, typei, typej);
}

/** @brief Get the value Z value (either U or F) at x between types typei and typej.
  * @details The function is currently private and accessed through get_U and get_F.
  * For production purposes it may be reasonable to move it to public, and 
  * remove get_U and get_F
  */
scalar PreComp_solver::finterpolate(std::vector<scalar> &Z, scalar x, int typei, int typej){
   //scalar y0, y1;
   scalar x0, x1;
   int index = 0;
   int index_l = x/Dx;
#ifdef DEBUG
   if (index_l < 0) 
     cout << "WTF?!" << endl; 
#endif

   // check that the index is not too high (all the tables are equally long): 
   if (index_l > n_values -2) {
      // cout << "returned zero for x: " << x << endl; 
      return 0.; 
   } 

   // sort so that typei <= typej
   if (typei > typej){
     int tmp = typei;
     typei = typej;
     typej = tmp;
   } 

   // get the index to read from:
   index = typei * ntypes;
   if (typei > 1) {
     index -= (typei*typei - typei)/2;
   }
   index += (typej - typei);

   index = index*n_values;
   index += index_l;


   // interpolate:
   x0 = index_l * Dx; 
   x1 = x0 + Dx; 
   return Z[index] + (Z[index +1] - Z[index])*(x - x0)/Dx;
}  

scalar PreComp_solver::get_field_energy(int index0, int index1) {
	// Sum over all field
   scalar energy = 0.0;
	if(index0 == -1 || index1 == -1) {
      for (int i = 0; i < num_threads*num_blobs*num_blobs; i++) {
        energy += fieldenergy[i]; 
      }
	} else if (index0 == index1) {
      for (int i = 0; i < num_threads; i++) {
        energy += fieldenergy[ (i*num_blobs + index0)*num_blobs + index1]; 
      } 
		// return fieldenergy[index0][index1];
	} else {
		// Order of blob indices is unknown in the calculations, so must add
      for (int i = 0; i < num_threads; i++) {
        energy += fieldenergy[ (i*num_blobs + index0)*num_blobs + index1]; 
        energy += fieldenergy[ (i*num_blobs + index1)*num_blobs + index0];
      } 
		// return fieldenergy[index0][index1] + fieldenergy[index1][index0];
	}

   return energy;
}


void PreComp_solver::build_pc_nearest_neighbour_lookup() {
   pcLookUp.clear(); 
   int x, y, z; 
   for (int i=0; i<n_beads; i++) {
     x = (int) floor(b_pos[3*i  ] / pcVoxelSize);
     y = (int) floor(b_pos[3*i+1] / pcVoxelSize);
     z = (int) floor(b_pos[3*i+2] / pcVoxelSize);

     pcLookUp.add_node_to_stack(i, x, y, z); 
   }
}


void PreComp_solver::prebuild_pc_nearest_neighbour_lookup_and_swap() { 
   pcLookUp.clear_shadow_layer();
   int x, y, z; 
   for (int i=0; i<n_beads; i++) {
     x = (int) floor(b_pos[3*i  ] / pcVoxelSize);
     y = (int) floor(b_pos[3*i+1] / pcVoxelSize);
     z = (int) floor(b_pos[3*i+2] / pcVoxelSize);

     pcLookUp.add_node_to_stack_shadow(i, x, y, z); 
   }
   pcLookUp.safely_swap_layers();
}


void PreComp_solver::prebuild_pc_nearest_neighbour_lookup() { 
   pcLookUp.clear_shadow_layer();
   pcLookUp.forbid_swapping();
   int x, y, z; 
   for (int i=0; i<n_beads; i++) {
     x = (int) floor(b_pos[3*i  ] / pcVoxelSize);
     y = (int) floor(b_pos[3*i+1] / pcVoxelSize);
     z = (int) floor(b_pos[3*i+2] / pcVoxelSize);

     pcLookUp.add_node_to_stack_shadow(i, x, y, z); 
   }
   pcLookUp.allow_swapping();
}


void PreComp_solver::safely_swap_pc_layers() {
   pcLookUp.safely_swap_layers();
}
  
void PreComp_solver::write_beads_to_file(FILE *fout, int timestep){
   const float toA = mesoDimensions::length * 1e10;
   fprintf(fout, "MODEL %12d\n", timestep); 
   for (int i=0; i<n_beads; i++){ 
      fprintf(fout, "ATOM %6d %4s %3s %1s%4i    %8.3f%8.3f%8.3f  %8d %3d\n",
               i+1, "CA", stypes[b_types[i]].c_str(), "A", i+1, 
               b_pos[3*i]*toA, b_pos[3*i+1]*toA, b_pos[3*i+2]*toA,
               b_elems_ndx[i], b_blob_ndx[i]); 
   } 
   fprintf(fout, "ENDMDL\n");
} 


/**@}*/

