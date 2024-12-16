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

#include <iostream>
#include <cmath>
#include <iomanip>
#include "VolumeIntersection.h"

using std::cout;
using std::endl; 

template <typename T>
bool exists(const std::array<T,3> &p, int ips, const std::array<std::array<T,3>, 56>& W){

  for (int i=0; i<ips; i++){ 
    if (distance(p, W[i]) < ffea_const::threeErr) {
      // cout << "discarded" << endl;
      return true;
    }
  } 
  return false;

}

/**  max volume for a number of points computed as the size of the cube 
 *    with the same side size as the double of the maximum distance found 
 *    between the given points.
 */
template <typename T>
T maxVolume(int ips, std::array<std::array<T,3>, 56> &W){

  std::array<T,3> c = {0,0,0};
  for (int i=0; i<ips; i++){
    add(c, W[i], c); 
  } 
  T s = 1./ips; 
  resize(s, c);
  T d0 = ffea_const::zero;
  T di; 
  for (int i=0; i<ips; i++){
    di = distance(c, W[i]); 
    if (di > d0) d0 = di;
  } 
  return ffea_const::eight*di*di*di; 
  // return ffea_const::volFactor*di*di*di; 

}

/**  max volume and area for a number of points computed as the size of the cube 
 *    with the same side size as the double of the maximum distance found 
 *    between the given points.
 */
template <typename T>
void maxVolumeAndArea(int ips, std::array<std::array<T, 3>, 56> &W, T &volume, T &area){

  std::array<T,3> c = {0,0,0};
  for (int i=0; i<ips; i++){
    add(c, W[i], c); 
  } 
  T s = 1./ips; 
  resize(s, c);
  T d0 = ffea_const::zero;
  T di; 
  for (int i=0; i<ips; i++){
    di = distance(c, W[i]); 
    if (di > d0) d0 = di;
  } 
  area = ffea_const::twentyfour*di*di;
  volume = ffea_const::eight*di*di*di; 
  // return ffea_const::volFactor*di*di*di; 

}


/**  CM for the first ips points stored in std::array<T,3> W[56] */ 
template <typename T>
void findCM(int ips, const std::array<std::array<T, 3>, 56> &W, std::array<T,3> &cm){
  initialise(cm); 
  for (int i=0; i<ips; i++){
    add(cm, W[i], cm); 
  } 
  T s = 1./ips; 
  resize(s, cm);

}


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 * the code assumes that tetA has its vertices in the right order. 
 */
template <typename T>
void getBAndN_Order(const std::array<std::array<T, 3>, 4> &tetA, int n0, int n1, int n2, const std::array<T,3> &t, std::array<T,3> &b, std::array<T,3> &n){

  // std::array<T,3>  t_j; 
  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        getNormal(tetA[0], tetA[2], tetA[1], b);
      } else {
        getNormal(tetA[3], tetA[0], tetA[1], b);
      }
    } else {
        getNormal(tetA[3], tetA[2], tetA[0], b);
    }
  } else {
     getNormal(tetA[3], tetA[1], tetA[2], b);
  } 
  getUnitNormal(b, t, n); 
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  std::array<T,3> aux; 
  add(tetA[n0], n, aux); 
  if (!sameSideLine(aux, tetA[n2], tetA[n1], tetA[n0])) resize(ffea_const::mOne, n); 
  // if we write the function we save a call to substract 
}


/** Given the face formed by tetA[0]:tetA[1]:tetA[2] and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 */
template <typename T>
void getBAndN(const std::array<std::array<T,3>, 4> &tetA, int n0, int n1, int n2, const std::array<T,3> &t, std::array<T,3> &b, std::array<T,3> &n){
  
   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   //   but t == pl1, so: 
   std::array<T,3> pl2, aux; 
   sub(tetA[n2], tetA[n0], pl2);
   // aux is a vector normal to the face:
   arr3arr3VectorProduct(t, pl2, aux);
   normalize(aux, b);
   add(aux, tetA[n0], aux); 
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   int n3 = getMissingNode(n0, n1, n2); 
   T d = - dot(b, tetA[n0]);
   T t1 = dot(tetA[n3], b) + d;
   T t2 = dot(aux, b) + d; // 1 + d; // std::array<T,3>std::array<T,3>DotProduct(b, b) + d   
   if (!sameSign(t1,t2)) resize(ffea_const::mOne, b);
 
  // Now we need N: 
  getUnitNormal(b, t, n);
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  add(tetA[n0], n, aux);
  if (!sameSideLine(aux, tetA[n2], tetA[n1], tetA[n0])) resize(ffea_const::mOne, n);
}

/** Given the face formed by f0, f1, f2 (knowing p3) and the tangent unit vector t:
 * get b: the normal to a face pointing inwards.
 *     n: the normal to t, on the face 
 */
template <typename T>
void getBAndN(const std::array<T,3> &f0, const std::array<T,3> &f1, const std::array<T,3> &f2, const std::array<T,3> &p3, const std::array<T,3> &t, std::array<T,3> &b, std::array<T,3> &n){
  
   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   //   but t == pl1, so: 
   std::array<T,3> pl2, aux; 
   sub(f2, f0, pl2);
   // aux is a vector normal to the face:
   arr3arr3VectorProduct(t, pl2, aux);
   normalize(aux, b);
   add(aux, f0, aux); 
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   T d = - dot(b, f0);
   T t1 = dot(p3, b) + d;
   T t2 = dot(aux, b) + d; 
   if (!sameSign(t1,t2)) resize(ffea_const::mOne, b);
 
  // Now we need N: 
  getUnitNormal(b, t, n);
  // and we'll check whether n is:
  //    at the same side of the plane tetA[n0]:tetA[n1]:tetA[n3] than tetA[n2]:
  // or even better inline sameSide2, and check whether:
  //    tetA[n2] and n lay on the same side 
  //       of the line given by t = unit(tetA[n1] - tetA[n0])
  add(f0, n, aux);
  if (!sameSideLine(aux, f2, f1, f0))
      resize(ffea_const::mOne, n);

}

/** Get the volume contribution of a node of type 1, i. e.,
 *    a node of the intersection polyhedron that already existed in one of the 
 *    intersectant thetrahedra. 
 */
template <typename T>
T volumeForNode(const std::array<std::array<T,3>, 4> &tetA, int node) {
  T volume = 0.0;

  // do the loop: 
  // for each edge:
  std::array<T,3> t_i, n_j, b_j; // t_j, t_ij
  for (int i=0; i<4; i++) { 
    if (i == node) continue; 
    tangent(tetA[i], tetA[node], t_i); 
    T pt = dot(tetA[node], t_i);
    // for each of the faces that this edge is adjacent to:
    for (int j=0; j<4; ++j) {
      if (j == node) continue;
      if (j == i) continue; 
      // getBAndN_Order(tetA, node, i, j, t_i, b_j, n_j); // both work well. 
      getBAndN(tetA, node, i, j, t_i, b_j, n_j);  // it works, too. 
      /////////////////
      /* // CHECK ! // 
      std::array<T,3> aux; 
      cout << "node: " << node << " i: " << i << " j: " << j << endl; 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      resize2(ffea_const::ten,t_i,aux); 
      add(tetA[node], aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      resize2(ffea_const::ten,b_j,aux); 
      add(tetA[node], aux, aux); 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      resize2(ffea_const::ten,n_j,aux); 
      add(tetA[node], aux, aux); 
      cout << "ip: " << tetA[node][0] << ", " << tetA[node][1] << ", " << tetA[node][2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! // 
      /////////////////
      volume += pt * dot(tetA[node], n_j) * dot(tetA[node], b_j);
    } 
  } 

  // cout << " volume: " << volume << endl; 
  return volume; 

}

template <typename T>
T volumeForNode(const std::array<T,3> &tet0, const std::array<T,3> &tet1, const std::array<T,3> &tet2, const std::array<T,3> &tetN) {
  T volume = 0.0;
  T pt; 
  std::array<T,3> t_i, n_j, b_j;

  // unfold the loop:
  // for each edge:
  tangent(tet0, tetN, t_i);
  pt = dot(tetN, t_i);

  getBAndN(tetN, tet0, tet1, tet2, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);
  getBAndN(tetN, tet0, tet2, tet1, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);



  tangent(tet1, tetN, t_i);
  pt = dot(tetN, t_i);

  getBAndN(tetN, tet1, tet0, tet2, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);
  getBAndN(tetN, tet1, tet2, tet0, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);


  tangent(tet2, tetN, t_i);
  pt = dot(tetN, t_i);

  getBAndN(tetN, tet2, tet0, tet1, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);
  getBAndN(tetN, tet2, tet1, tet0, t_i, b_j, n_j); 
  volume += pt * dot(tetN, n_j) * dot(tetN, b_j);



  return volume; 
}

/** Get the volume and area contribution of a node of type 1, i. e.,
 *    a node of the intersection polyhedron that already existed in one of the 
 *    intersectant thetrahedra. 
 */
template <typename T>
void volumeAndAreaForNode(const std::array<std::array<T,3>, 4> &tetA, int node, T &volume, T &area) {
  volume = 0.0;
  area = 0.0; 

  // do the loop: 
  // for each edge:
  std::array<T,3> t_i, n_j, b_j; // t_j, t_ij
  for (int i=0; i<4; i++) { 
    if (i == node) continue; 
    tangent(tetA[i], tetA[node], t_i); 
    T pt = dot(tetA[node], t_i);
    // for each of the faces that this edge is adjacent to:
    for (int j=0; j<4; ++j) {
      if (j == node) continue;
      if (j == i) continue; 
      // getBAndN_Order(tetA, node, i, j, t_i, b_j, n_j); // both work well. 
      getBAndN(tetA, node, i, j, t_i, b_j, n_j);  // it works, too. 
      T pt_pn = pt * dot(tetA[node], n_j);
      area += pt_pn; 
      volume += pt_pn * dot(tetA[node], b_j);
    } 
  } 
}


/** Calculate the volume contribution of the intersection point ip,
 *    resulting of the intersection between edge tetA[e1]->tetA[e2]
 *    and the face given by the points tetB[f1]:tetB[f2]:tetB[f3]. 
 */ 
template <typename T>
T volumeForIntPoint(const std::array<T,3> &ip, const std::array<std::array<T,3>, 4> &tetA, int e1, int e2, const std::array<std::array<T,3>, 4>& tetB, int f1, int f2, int f3){

   T volume = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // std::array<T,3> B[3] will store the unit normals to these faces. 
   std::array<std::array<T,3>,3> B; 
   //
   // ip has three edges:
   //   tetA[e1]->tetA[e2],
   //   the intersection between face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   //   the intersection between face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   // std::array<T,3> tngt[3] will store the unit vectors of these edges.
   std::array<std::array<T, 3>, 3> tngt;
   //
   // Finally, we will go through the loop edges x faces. 
   //
   // Let's start:
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards(tetB, f1, f2, f3, B[0]); 
   //  and the easy tangent:  
   tangent(tetA[e2], tetA[e1], tngt[0]);
   // but check whether we are entering or escaping: 
   int f4 = getMissingNode(f1,f2,f3); 
   //   so if escaping:
   if (sameSidePlane(tetA[e1], tetB[f4], tetB[f1], tetB[f2], tetB[f3])) 
     resize(ffea_const::mOne, tngt[0]);
 
   // Secondly, normals and tangents for the tetA faces: 
   int cnt = 0; 
   std::array<int, 2> F;
   for (int i=0; i<4; i++){ 
     if ((i != e1) && (i != e2)) {
       // store "i" so we'll know that B[cnt] will correspond to the face triad e1:e2:i
       F[cnt] = i;
       cnt++; 
       getNormalInwards(tetA, e1, e2, i, B[cnt]);
       getUnitNormal(B[0],B[cnt], tngt[cnt]);
       // but obviously tngt[cnt] may not be in the right sense: 
       //   we know the intersection point, 
       //           and that it comes from the edge tetA[e1]:tetA[e2]
       //   Hence, for
       //   the intersection between face ( tetA[e1]->tetA[i] ) : (tetA[e1]->tetA[e2]) 
       //       and face tetB[f1]:tetB[f2]:tetB[f3];
       //   tngt[1] has to point towards... 
       int e4 = getMissingNode(i, e1, e2);
       std::array<T,3> vaux;
       add(ip, tngt[cnt], vaux);
       if (!sameSidePlane( vaux , tetA[i], tetA[e4], tetA[e1], tetA[e2])) 
                 resize(ffea_const::mOne, tngt[cnt]);
     }
   }

  /////////////////////
  ////// CHECK ////////
  /*
  std::array<T,3> aux; 
  std::array<T,3> C[3]; 
  faceCentroid(tetB[f1], tetB[f2], tetB[f3], C[0]);
  faceCentroid(tetA[e1], tetA[e2], tetA[F[0]], C[1]);
  faceCentroid(tetA[e1], tetA[e2], tetA[F[1]], C[2]);
  add(C[0], B[0], aux); 
  cout << "C: " << C[0][0] << ", " << C[0][1] << ", " << C[0][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 

  add(C[1], B[1], aux); 
  cout << "C: " << C[1][0] << ", " << C[1][1] << ", " << C[1][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 
  add(C[2], B[2], aux); 
  cout << "C: " << C[2][0] << ", " << C[2][1] << ", " << C[2][2] << endl; 
  cout << "B: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
  cout << endl << endl; 
  */ 
  ///////////////////////////
  /*
  for (int i=0; i<3; i++) { 
    add(ip, tngt[i], aux); 
    cout << "P: " << ip[0] << " " << ip[1] << " " << ip[2] << endl; 
    cout << "T: " << aux[0] << " " << aux[1] << " " << aux[2] << endl; 
    cout << endl << endl; 
  } 
  */ 
  ////// CHECK ////////
  /////////////////////
  
   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., tngt[0] with faces B[1] and B[2]
   std::array<T,3> n; 
   std::array<T, 3> ipT, ipB;
   for (int i=0; i<3; i++) {
      ipT[i] = dot(ip,tngt[i]); 
      ipB[i] = dot(ip,B[i]); 
   } 
   for (int i=0; i<2; i++) { 
      // for B[1] we need N
      getUnitNormal(B[i+1], tngt[0], n); 
      // and we'll check whether n lays on the same side as tetA[F[i]]
      //    of the line given by tngt[0] = unit(tetA[e2] - tetA[e1])
      std::array<T,3> vaux;
      add(ip, n, vaux);
      if (!sameSideLine(vaux, tetA[F[i]], tetA[e2], tetA[e1])) resize(ffea_const::mOne, n); 
      /////////////////
      /*// CHECK ! // 
      std::array<T,3> aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      std::array<T,3>Resize2(ffea_const::ten,tngt[0],aux); 
      add(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,B[i+1],aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,n,aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      // CHECK ! // */
      /////////////////
      volume += ipT[0] * dot(ip, n) * ipB[i+1]; //  std::array<T,3>std::array<T,3>DotProduct(ip, B[i+1]);
   } 
  
   // Now we need the contribution of the edges tngt[1] and tngt[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   std::array<std::array<T,3>, 3> V;
   add(ip, tngt[1], V[0]);
   add(ip, tngt[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First tngt[i+1] with B[0];
      // get the normal:
      getUnitNormal(B[0], tngt[i+1], n);
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + tngt[i]*t  
      std::array<T,3> v;
      add(ip, n, v);
      if (!sameSideLine(v, V[i+1], V[i], ip)) resize(ffea_const::mOne, n); 
      volume += ipT[i+1] * dot(ip, n) * ipB[0];   

      /////////////////
      /*// CHECK ! // 
      std::array<T,3> aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      std::array<T,3>Resize2(ffea_const::ten,tngt[i+1],aux); 
      add(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,B[0],aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,n,aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! //
      /////////////////

      // and now tngt[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal(B[i+1], tngt[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (dot(n,B[0]) < 0) resize(ffea_const::mOne, n); 
      volume += ipT[i+1] * dot(ip, n) * ipB[i+1];
      /////////////////
      /* // CHECK ! // 
      std::array<T,3> aux; 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      std::array<T,3>Resize2(ffea_const::ten,tngt[i+1],aux); 
      add(ip, aux, aux); 
      cout << "t: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,B[i+1],aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "b: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      std::array<T,3>Resize2(ffea_const::ten,n,aux); 
      add(ip, aux, aux); 
      cout << "ip: " << ip[0] << ", " << ip[1] << ", " << ip[2] << endl; 
      cout << "n: " << aux[0] << ", " << aux[1] << ", " << aux[2] << endl; 
      cout << endl << endl; 
      */ // CHECK ! //
      /////////////////

   } 
      
   // cout << " volume: " << volume << endl; 
   return volume;
}

/** Calculate the volume contribution of the intersection point ip,
 *    resulting of the intersection between edge tetAe1->tetAe2
 *    and the face given by the points tetBf1:tetBf2:tetBf3,
 *    while knowing the rest of the tetrahedra points.
 */
template <typename T>
T volumeForIntPointII(const std::array<T,3> &ip, const std::array<T,3> &tetAe1, const std::array<T,3> &tetAe2, const std::array<T,3> &tetAe3, const std::array<T,3> &tetAe4, const std::array<T,3> &tetBf1, const std::array<T,3> &tetBf2, const std::array<T,3> &tetBf3, const std::array<T,3> &tetBf4){

   T volume = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // std::array<T,3> B[3] will store the unit normals to these faces. 
   std::array<std::array<T,3>, 3> B; 
   //
   // ip has three edges:
   //   tetA[e1]->tetA[e2],
   //   the intersection between face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   //   the intersection between face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   // std::array<T,3> tngt[3] will store the unit vectors of these edges. 
   std::array<std::array<T,3>, 3> tngt;
   //
   // Finally, we will go through the loop edges x faces. 
   //
   // Let's start:
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards(tetBf1, tetBf2, tetBf3, tetBf4, B[0]); 
   //  and the easy tangent:  
   tangent(tetAe2, tetAe1, tngt[0]); 
   // but check whether we are entering or escaping: 
   //   so if escaping:
   if (sameSidePlane(tetAe1, tetBf4, tetBf1, tetBf2, tetBf3)) 
     resize(ffea_const::mOne, tngt[0]); 
 
   std::array<T,3> vaux;
   // Secondly, normals and tangents for the tetA faces: 
   getNormalInwards(tetAe1, tetAe2, tetAe3, tetAe4, B[1]);
   getUnitNormal(B[0],B[1], tngt[1]);
   add(ip, tngt[1], vaux);
   if (!sameSidePlane( vaux , tetAe3, tetAe4, tetAe1, tetAe2)) 
             resize(ffea_const::mOne, tngt[1]);


   getNormalInwards(tetAe1, tetAe2, tetAe4, tetAe3, B[2]);
   getUnitNormal(B[0],B[2],tngt[2]);
   add(ip, tngt[2], vaux);
   if (!sameSidePlane( vaux , tetAe4, tetAe3, tetAe1, tetAe2)) 
             resize(ffea_const::mOne, tngt[2]);


  
   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., tngt[0] with faces B[1] and B[2]
   std::array<T,3> n;
   std::array<T,3> ipT, ipB;
   for (int i=0; i<3; i++) {
      ipT[i] = dot(ip,tngt[i]); 
      ipB[i] = dot(ip,B[i]); 
   } 
   getUnitNormal(B[1], tngt[0], n); 
   // and we'll check whether n lays on the same side as tetA[F[i]]
   //    of the line given by tngt[0] = unit(tetA[e2] - tetA[e1])
   add(ip, n, vaux);
   if (!sameSideLine(vaux, tetAe3, tetAe2, tetAe1)) resize(ffea_const::mOne, n); 
   volume += ipT[0] * dot(ip, n) * ipB[1];

   getUnitNormal(B[2], tngt[0], n); 
   // and we'll check whether n lays on the same side as tetA[F[i]]
   //    of the line given by tngt[0] = unit(tetA[e2] - tetA[e1])
   add(ip, n, vaux);
   if (!sameSideLine(vaux, tetAe4, tetAe2, tetAe1)) resize(ffea_const::mOne, n); 
   volume += ipT[0] * dot(ip, n) * ipB[2];

  
   // Now we need the contribution of the edges tngt[1] and tngt[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   std::array<std::array<T,3>, 3> V;
   add(ip, tngt[1], V[0]);
   add(ip, tngt[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First tngt[i+1] with B[0];
      // get the normal:
      getUnitNormal(B[0], tngt[i+1], n); 
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + tngt[i]*t  
      std::array<T,3> v;
      add(ip, n, v);
      if (!sameSideLine(v, V[i+1], V[i], ip)) resize(ffea_const::mOne, n); 
      volume += ipT[i+1] * dot(ip, n) * ipB[0];   

      // and now tngt[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal(B[i+1],tngt[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (dot(n,B[0]) < 0) resize(ffea_const::mOne, n); 
      volume += ipT[i+1] * dot(ip, n) * ipB[i+1];

   } 
      
   // cout << " volume: " << volume << endl; 
   return volume;
}


/** Calculate the volume and area contribution of the intersection point ip,
 *    resulting of the intersection between edge tetA[e1]->tetA[e2]
 *    and the face given by the points tetB[f1]:tetB[f2]:tetB[f3]. 
 */ 
template <typename T>
void volumeAndAreaForIntPoint(const std::array<T,3> &ip, const std::array<std::array<T,3>, 4> &tetA, int e1, int e2, const std::array<std::array<T,3>, 4> &tetB, int f1, int f2, int f3, T &volume, T &area){
   volume = 0; 
   area = 0; 
   // ip has three faces:
   //   tetB[f1]:tetB[f2]:tetB[f3] 
   //   face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2])
   //   face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2])
   // std::array<T,3> B[3] will store the unit normals to these faces. 
   std::array<std::array<T,3>, 3> B;
   //
   // ip has three edges:
   //   tetA[e1]->tetA[e2],
   //   the intersection between face ( tetA[e1]->tetA[e3] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   //   the intersection between face ( tetA[e1]->tetA[e4] ) : (tetA[e1]->tetA[e2]) 
   //       and face tetB[f1]:tetB[f2]:tetB[f3];
   //      Its tangent unit vector is the cross product of the normals of the planes.
   // std::array<T,3> tngt[3] will store the unit vectors of these edges.
   std::array<std::array<T,3>, 3> tngt;
   //
   // Finally, we will go through the loop edges x faces. 
   //
   // Let's start:
   // Firstly, the normal to tetB[f1]:tetB[f2]:tetB[f3]:
   getNormalInwards(tetB, f1, f2, f3, B[0]); 
   //  and the easy tangent:  
   tangent(tetA[e2], tetA[e1], tngt[0]); 
   // but check whether we are entering or escaping: 
   int f4 = getMissingNode(f1,f2,f3); 
   //   so if escaping:
   if (sameSidePlane(tetA[e1], tetB[f4], tetB[f1], tetB[f2], tetB[f3])) 
     resize(ffea_const::mOne, tngt[0]); 
 
   // Secondly, normals and tangents for the tetA faces: 
   int cnt = 0; 
   int F[2];
   for (int i=0; i<4; i++){ 
     if ((i != e1) && (i != e2)) {
       // store "i" so we'll know that B[cnt] will correspond to the face triad e1:e2:i
       F[cnt] = i;
       cnt++; 
       getNormalInwards(tetA, e1, e2, i, B[cnt]);
       getUnitNormal(B[0],B[cnt],tngt[cnt]);
       // but obviously tngt[cnt] may not be in the right sense: 
       //   we know the intersection point, 
       //           and that it comes from the edge tetA[e1]:tetA[e2]
       //   Hence, for
       //   the intersection between face ( tetA[e1]->tetA[i] ) : (tetA[e1]->tetA[e2]) 
       //       and face tetB[f1]:tetB[f2]:tetB[f3];
       //   tngt[1] has to point towards... 
       int e4 = getMissingNode(i, e1, e2);
       std::array<T,3> vaux;
       add(ip, tngt[cnt], vaux);
       if (!sameSidePlane( vaux , tetA[i], tetA[e4], tetA[e1], tetA[e2])) 
                 resize(ffea_const::mOne, tngt[cnt]);
     }
   }

   // Now we'll get the (partial) volumes: 
   //  Firstly the contribution of the edge tetA[e1]:tetA[e2],
   //     i. e., tngt[0] with faces B[1] and B[2]
   std::array<T,3> n;
   std::array<T,3> ipT, ipB;
   for (int i=0; i<3; i++) {
      ipT[i] = dot(ip,tngt[i]); 
      ipB[i] = dot(ip,B[i]); 
   } 
   for (int i=0; i<2; i++) { 
      // for B[1] we need N
      getUnitNormal(B[i+1], tngt[0], n); 
      // and we'll check whether n lays on the same side as tetA[F[i]]
      //    of the line given by tngt[0] = unit(tetA[e2] - tetA[e1])
      std::array<T,3> vaux;
      add(ip, n, vaux);
      if (!sameSideLine(vaux, tetA[F[i]], tetA[e2], tetA[e1])) resize(ffea_const::mOne, n); 
      T ipT_ipN = ipT[0] * dot(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[i+1];
   } 
  
   // Now we need the contribution of the edges tngt[1] and tngt[2] with
   //           - face tetB[f1]:tetB[f2]:tetB[f3] (B[0]);
   //           - either face tet[e1]->tet[e3] 
   //           - or face tet[e1]->tet[e2];
   //
   std::array<std::array<T,3>, 3> V;
   add(ip, tngt[1], V[0]);
   add(ip, tngt[2], V[1]);
   for (int i=0; i<3; i++) V[2][i] = V[0][i]; 
   for (int i=0; i<2; i++) {
      // First tngt[i+1] with B[0];
      // get the normal:
      getUnitNormal(B[0], tngt[i+1], n); 
      // and we'll check whether n lays on the same side as tetA[F[i+1]] 
      //    of the line given by R(t) = ip + tngt[i]*t  
      std::array<T,3> v;
      add(ip, n, v);
      if (!sameSideLine(v, V[i+1], V[i], ip)) resize(ffea_const::mOne, n); 
      T ipT_ipN = ipT[i+1] * dot(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[0];   

      // and now tngt[i+1] with face tetA[e1]:tetA[e2]:tetA[F[i]]
      // get the normal: 
      getUnitNormal(B[i+1],tngt[i+1],n);
      // and we'll check whether n points towards B[0]: 
      if (dot(n,B[0]) < 0) resize(ffea_const::mOne, n); 
      ipT_ipN = ipT[i+1] * dot(ip, n);
      area += ipT_ipN; 
      volume += ipT_ipN * ipB[i+1];
   } 
}


/** Return the volume intersection between two tetrahedra */ 
template <typename T>
T volumeIntersection(const std::array<std::array<T,3>, 4> &tetA, const std::array<std::array<T,3>, 4> &tetB, bool calcCM, std::array<T,3> &cm){
  T vol = 0.0; 
  std::array<std::array<T,3>, 56> W; 

  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists(tetA[i], ips, W) && nodeInTet(tetA[i], tetB)) { 
      store(tetA[i], W[ips]);
      ips += 1; 
      vol += volumeForNode(tetA, i);
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists(tetB[i], ips, W) && nodeInTet(tetB[i], tetA)) { 
      store(tetB[i], W[ips]);
      ips += 1;
      vol += volumeForNode(tetB, i);
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  std::array<T,3> ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /* cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetA, i, j, tetB, 0, 1, 2); 
        } 
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetA, i, j, tetB, 1, 2, 3); 
        } 
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetA, i, j, tetB, 2, 3, 0); 
        }
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetA, i, j, tetB, 3, 0, 1); 
        }
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetB, i, j, tetA, 0, 1, 2); 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetB, i, j, tetA, 1, 2, 3); 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetB, i, j, tetA, 2, 3, 0); 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          vol += volumeForIntPoint(ip, tetB, i, j, tetA, 3, 0, 1); 
        }
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) 
    return 0.; 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  T w = maxVolume(ips, W); 
  if ((fabs(vol) > w) || (vol < 0)) { 
    vol = w;
  }

  if (calcCM)
      findCM(ips, W, cm);

  return vol; 
 

} 

template <typename T>
void contribVolForNode(T &vol, const std::array<T,3> &n0, const std::array<T,3> &n1, const std::array<T,3> &n2, const std::array<T,3> &n3, std::array<std::array<T,3>,56> &W, int &ips){
    store(n0, W[ips]);
    ips += 1; 
    vol += volumeForNode(n1, n2, n3, n0);
}


template <typename T>
void contribVolForIntPoint(T &vol, const std::array<T,3> &ip, const std::array<T,3> &tetAi, const std::array<T,3> &tetAj, const std::array<T,3> &tetAe3, const std::array<T,3> &tetAe4, const std::array<T,3> &tetB0, const std::array<T,3> &tetB1, const std::array<T,3> &tetB2, const std::array<T,3> &tetB3, std::array<std::array<T,3>,56> &W, int &ips){
    store(ip, W[ips]);
    ips +=1;
    vol += volumeForIntPointII(ip, tetAi, tetAj, tetAe3, tetAe4,
        tetB0, tetB1, tetB2, tetB3);
}

template <typename T>
void contribVolForIntersections(T &vol, const std::array<T,3> &tetAi, const std::array<T,3> &tetAj, const std::array<T,3> &tetAe3, const std::array<T,3> &tetAe4, const std::array<T,3> &tetB0, const std::array<T,3> &tetB1, const std::array<T,3> &tetB2, const std::array<T,3> &tetB3, std::array<std::array<T,3>,56> &W, int &ips) {
      std::array<T,3> ip; 
      if ( (intersectionPoint(ip, tetAi, tetAj, tetB0, tetB1, tetB2)) && (!exists(ip, ips, W)) ) 
         contribVolForIntPoint(vol, ip, tetAi, tetAj, tetAe3, tetAe4, 
                                   tetB0, tetB1, tetB2, tetB3, W, ips);

      if ( (intersectionPoint(ip, tetAi, tetAj, tetB1, tetB2, tetB3)) && (!exists(ip, ips, W)) )
         contribVolForIntPoint(vol, ip, tetAi, tetAj, tetAe3, tetAe4, 
                                   tetB1, tetB2, tetB3, tetB0, W, ips);
      
      if ( (intersectionPoint(ip, tetAi, tetAj, tetB2, tetB3, tetB0)) && (!exists(ip, ips, W)) )
         contribVolForIntPoint(vol, ip, tetAi, tetAj, tetAe3, tetAe4, 
                                   tetB2, tetB3, tetB0, tetB1, W, ips);
      
      if ( (intersectionPoint(ip, tetAi, tetAj, tetB3, tetB0, tetB1)) && (!exists(ip, ips, W)) )
         contribVolForIntPoint(vol, ip, tetAi, tetAj, tetAe3, tetAe4, 
                                   tetB3, tetB0, tetB1, tetB2, W, ips);
}

/** Return the volume intersection between two tetrahedra */ 
template <typename T>
T volumeIntersectionII(const std::array<T,3> &tetA0, const std::array<T,3> &tetA1, const std::array<T,3> &tetA2, const std::array<T,3> &tetA3,
    const std::array<T,3> &tetB0, const std::array<T,3> &tetB1, const std::array<T,3> &tetB2, const std::array<T,3> &tetB3, bool calcCM, std::array<T,3> &cm){

  T vol = 0.0; 
  std::array<std::array<T,3>,56> W; 

  int ips = 0;
  int an0, an1, an2; 
  // Check for interior points. 
  //   Unfold the loop: 
  // if point tetA[i] is inside tetB -> account for its contribution. 
  if (!exists(tetA0, ips, W) && nodeInTet(tetA0, tetB0, tetB1, tetB2, tetB3)) contribVolForNode(vol, tetA0, tetA1, tetA2, tetA3, W, ips); 

  if (!exists(tetA1, ips, W) && nodeInTet(tetA1, tetB0, tetB1, tetB2, tetB3)) contribVolForNode(vol, tetA1, tetA0, tetA2, tetA3, W, ips); 
  
  if (!exists(tetA2, ips, W) && nodeInTet(tetA2, tetB0, tetB1, tetB2, tetB3)) contribVolForNode(vol, tetA2, tetA0, tetA1, tetA3, W, ips); 
  
  if (!exists(tetA3, ips, W) && nodeInTet(tetA3, tetB0, tetB1, tetB2, tetB3)) contribVolForNode(vol, tetA3, tetA0, tetA1, tetA2, W, ips);


  // if point tetB[i] is inside tetA -> account for its contribution. 
  // PENDING: what happens if a node belongs to both tetA and tetB? 
  if (!exists(tetB0, ips, W) && nodeInTet(tetB0, tetA0, tetA1, tetA2, tetA3)) contribVolForNode(vol, tetB0, tetB1, tetB2, tetB3, W, ips); 
  
  if (!exists(tetB1, ips, W) && nodeInTet(tetB1, tetA0, tetA1, tetA2, tetA3)) contribVolForNode(vol, tetB1, tetB0, tetB2, tetB3, W, ips); 
  
  if (!exists(tetB2, ips, W) && nodeInTet(tetB2, tetA0, tetA1, tetA2, tetA3)) contribVolForNode(vol, tetB2, tetB0, tetB1, tetB3, W, ips); 
  
  if (!exists(tetB3, ips, W) && nodeInTet(tetB3, tetA0, tetA1, tetA2, tetA3)) contribVolForNode(vol, tetB3, tetB0, tetB1, tetB2, W, ips); 


  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge tetX[i]-tetX[j] i<j:
  contribVolForIntersections(vol, tetA0, tetA1, tetA2, tetA3, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB0, tetB1, tetB2, tetB3, tetA0, tetA1, tetA2, tetA3, W, ips); 
  
  contribVolForIntersections(vol, tetA0, tetA2, tetA1, tetA3, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB0, tetB2, tetB1, tetB3, tetA0, tetA1, tetA2, tetA3, W, ips); 
  
  contribVolForIntersections(vol, tetA0, tetA3, tetA1, tetA2, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB0, tetB3, tetB1, tetB2, tetA0, tetA1, tetA2, tetA3, W, ips); 
  
  contribVolForIntersections(vol, tetA1, tetA2, tetA0, tetA3, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB1, tetB2, tetB0, tetB3, tetA0, tetA1, tetA2, tetA3, W, ips); 
  
  contribVolForIntersections(vol, tetA1, tetA3, tetA0, tetA2, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB1, tetB3, tetB0, tetB2, tetA0, tetA1, tetA2, tetA3, W, ips); 
  
  contribVolForIntersections(vol, tetA2, tetA3, tetA0, tetA1, tetB0, tetB1, tetB2, tetB3, W, ips); 
  contribVolForIntersections(vol, tetB2, tetB3, tetB0, tetB1, tetA0, tetA1, tetA2, tetA3, W, ips); 
  


  // Done! Now multiply the current value by the magic factor:
  vol *= - ffea_const::oneOverSix;
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) 
    return 0.; 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  T w = maxVolume(ips, W); 
  if ((fabs(vol) > w) || (vol < 0)) { 
    vol = w;
  }

  if (calcCM) findCM(ips, W, cm);

  return vol; 
 

} 


/** Return the volume and area of intersection between two tetrahedra */ 
//  input: std::array<T,3> (&tetA)[4], std::array<T,3> (&tetB)[4])
//  output: T vol and T area. 
template <typename T>
void volumeAndAreaIntersection(const std::array<std::array<T,3>,4> &tetA, const std::array<std::array<T,3>,4> &tetB, T &vol, T &area){


  vol = 0.0; 
  T v_aux, a_aux;
  area = 0.0;
  std::array<std::array<T,3>, 56> W; 

  int ips = 0;
  // Check for interior points. 
  for (int i=0; i<4; i++) {
    // if point tetA[i] is inside tetB -> account for its contribution. 
    if (!exists(tetA[i], ips, W) && 
            nodeInTet(tetA[i], tetB[0], tetB[1], tetB[2], tetB[3])) { 
      // cout << "tetA-node[" << i << "] found int tetB " << endl; 
      store(tetA[i], W[ips]);
      ips += 1; 
      volumeAndAreaForNode(tetA, i, v_aux, a_aux); 
      vol += v_aux;
      area += a_aux; 
    } 
    // if point tetB[i] is inside tetA -> account for its contribution. 
    // PENDING: what happens if a node belongs to both tetA and tetB? 
    if (!exists(tetB[i], ips, W) && nodeInTet(tetB[i], tetA)) { 
      // cout << "tetB-node[" << i << "] found int tetA " << endl; 
      store(tetB[i], W[ips]);
      ips += 1;
      volumeAndAreaForNode(tetB, i, v_aux, a_aux); 
      vol += v_aux;
      area += a_aux; 
    } 
  } 

  
  // Check for new points coming from intersections: 
  //    We need to check every edge-face pair belonging to tetA and tetB. 
  // Thus, for every edge:
  
  std::array<T,3> ip; 
  for (int i=0; i<4; i++) { 
    for (int j=i+1; j<4; j++) {
      // check intersection for edge tetA:ij and every tetB face: 
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 0, 1, 2)) {
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[0]:tetB[1]:tetB[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetA, i, j, tetB, 0, 1, 2, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 1, 2, 3)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[1]:tetB[2]:tetB[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetA, i, j, tetB, 1, 2, 3, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 2, 3, 0)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[2]:tetB[3]:tetB[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetA, i, j, tetB, 2, 3, 0, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetA[i], tetA[j], tetB, 3, 0, 1)){
        /*cout << "intersection between edge tetA[" << i << "]:tetA[" << j << "] and " 
             << "face tetB[3]:tetB[0]:tetB[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetA, i, j, tetB, 3, 0, 1, v_aux, a_aux); 
          vol += v_aux; 
          area += a_aux; 
        }
      }

      // check intersection for edge tetB:ij and every tetA face: 
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 0, 1, 2)) {
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[0]:tetA[1]:tetA[2] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetB, i, j, tetA, 0, 1, 2, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 1, 2, 3)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[1]:tetA[2]:tetA[3] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetB, i, j, tetA, 1, 2, 3, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 2, 3, 0)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[2]:tetA[3]:tetA[0] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetB, i, j, tetA, 2, 3, 0, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
      if (intersectionPoint(ip, tetB[i], tetB[j], tetA, 3, 0, 1)){
        /*cout << "intersection between edge tetB[" << i << "]:tetB[" << j << "] and " 
             << "face tetA[3]:tetA[0]:tetA[1] " << endl; 
        cout << "   happening at: " << ip[0] << ":" << ip[1] << ":" << ip[2] << endl;*/
        if (!exists(ip, ips, W)) {
          store(ip, W[ips]);
          ips +=1;
          volumeAndAreaForIntPoint(ip, tetB, i, j, tetA, 3, 0, 1, v_aux, a_aux); 
          vol += v_aux;
          area += a_aux; 
        }
      }
    } 
  } 

  vol *= - ffea_const::oneOverSix;
  area *= ffea_const::half; 
  // CHECKing /// 
  // three points cannot have a volume.
  if (ips <= 3) {
    vol = 0.;
    area = 0.;
    return; 
  } 
 
  // the volume cannot be larger than the volume of the sphere with maximum radius. 
  maxVolumeAndArea(ips, W, v_aux, a_aux);
  if ((fabs(vol) > v_aux) || (vol < 0)) { 
    vol = v_aux;
    area = a_aux;
  }

} 


// INSTANTIATE everything to arr3 and eventually to grr3: 

template bool exists<scalar>(const std::array<scalar, 3>& p, int ips, const std::array<std::array<scalar, 3>, 56>& W);
template scalar maxVolume<scalar>(int ips, std::array<std::array<scalar, 3>, 56>& W);
template void maxVolumeAndArea<scalar>(int ips, std::array<std::array<scalar, 3>, 56>& W, scalar& volume, scalar& area);
template void findCM<scalar>(int ips, const std::array<std::array<scalar, 3>, 56>& W, std::array<scalar, 3>& cm);
template void getBAndN_Order<scalar>(const  std::array<std::array<scalar, 3>, 4>& tetA, int n0, int n1, int n2, const std::array<scalar, 3>& t, std::array<scalar, 3>& b, std::array<scalar, 3>& n);
template void getBAndN<scalar>(const std::array<std::array<scalar, 3>, 4>& tetA, int n0, int n1, int n2, const std::array<scalar, 3>& t, std::array<scalar, 3>& b, std::array<scalar, 3>& n);
template void getBAndN<scalar>(const std::array<scalar, 3>& f0, const std::array<scalar, 3>& f1, const std::array<scalar, 3>& f2, const std::array<scalar, 3>& p3, const std::array<scalar, 3>& t, std::array<scalar, 3>& b, std::array<scalar, 3>& n);
template scalar volumeForNode<scalar>(const std::array<std::array<scalar, 3>, 4>& tetA, int node);
template scalar volumeForNode<scalar>(const std::array<scalar, 3>& tet0, const std::array<scalar, 3>& tet1, const std::array<scalar, 3>& tet2, const std::array<scalar, 3>& tetN);
template void volumeAndAreaForNode<scalar>(const std::array<std::array<scalar, 3>, 4>& tetA, int node, scalar& volume, scalar& area);
template scalar volumeForIntPoint<scalar>(const std::array<scalar, 3>& ip, const std::array<std::array<scalar, 3>, 4>& tetA, int e1, int e2, const std::array<std::array<scalar, 3>, 4>& tetB, int f1, int f2, int f3);
template scalar volumeForIntPointII<scalar>(const std::array<scalar, 3>& ip, const std::array<scalar, 3>& tetAe1, const std::array<scalar, 3>& tetAe2, const std::array<scalar, 3>& tetAe3, const std::array<scalar, 3>& tetAe4, const std::array<scalar, 3>& tetBf1, const std::array<scalar, 3>& tetBf2, const std::array<scalar, 3>& tetBf3, const std::array<scalar, 3>& tetBf4);
template void volumeAndAreaForIntPoint<scalar>(const std::array<scalar, 3>& ip, const std::array<std::array<scalar, 3>, 4>& tetA, int e1, int e2, const std::array<std::array<scalar, 3>, 4>& tetB, int f1, int f2, int f3, scalar& volume, scalar& area);
template scalar volumeIntersection<scalar>(const std::array<std::array<scalar, 3>, 4>& tetA, const std::array<std::array<scalar, 3>, 4>& tetB, bool calcCM, std::array<scalar, 3>& cm);
template scalar volumeIntersectionII<scalar>(const std::array<scalar, 3>& tetA0, const std::array<scalar, 3>& tetA1, const std::array<scalar, 3>& tetA2, const std::array<scalar, 3>& tetA3, const std::array<scalar, 3>& tetB0, const std::array<scalar, 3>& tetB1, const std::array<scalar, 3>& tetB2, const std::array<scalar, 3>& tetB3, bool calcCM, std::array<scalar, 3>& cm);

template void volumeAndAreaIntersection<scalar>(const std::array<std::array<scalar, 3>, 4>& tetA, const std::array<std::array<scalar, 3>, 4>& tetB, scalar& vol, scalar& area);
template void contribVolForNode<scalar>(scalar& vol, const std::array<scalar, 3>& n0, const std::array<scalar, 3>& n1, const std::array<scalar, 3>& n2, const std::array<scalar, 3>& n3, std::array<std::array<scalar, 3>, 56>& W, int& ips);
template void contribVolForIntPoint<scalar>(scalar& vol, const std::array<scalar, 3>& ip, const std::array<scalar, 3>& tetAi, const std::array<scalar, 3>& tetAj, const std::array<scalar, 3>& tetAe3, const std::array<scalar, 3>& tetAe4, const std::array<scalar, 3>& tetB0, const std::array<scalar, 3>& tetB1, const std::array<scalar, 3>& tetB2, const std::array<scalar, 3>& tetB3, std::array<std::array<scalar, 3>, 56>& W, int& ips);
template void contribVolForIntersections<scalar>(scalar& vol, const std::array<scalar, 3>& tetAi, const std::array<scalar, 3>& tetAj, const std::array<scalar, 3>& tetAe3, const std::array<scalar, 3>& tetAe4, const std::array<scalar, 3>& tetB0, const std::array<scalar, 3>& tetB1, const std::array<scalar, 3>& tetB2, const std::array<scalar, 3>& tetB3, std::array<std::array<scalar, 3>, 56>& W, int& ips);

#ifndef USE_DOUBLE
template bool exists<geoscalar>(const std::array<geoscalar,3> &p, int ips, const std::array<std::array<geoscalar,3>,56> &W);
template geoscalar maxVolume<geoscalar>(int ips, std::array<std::array<geoscalar,3>,56> &W);
template void maxVolumeAndArea<geoscalar>(int ips, std::array<std::array<geoscalar,3>,56> &W, geoscalar &volume, geoscalar &area);
template void findCM<geoscalar>(int ips, const std::array<std::array<geoscalar,3>,56> &W, std::array<geoscalar,3> &cm);
template void getBAndN_Order<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, int n0, int n1, int n2, const std::array<geoscalar,3> &t, std::array<geoscalar,3> &b, std::array<geoscalar,3> &n);
template void getBAndN<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, int n0, int n1, int n2, const std::array<geoscalar,3> &t, std::array<geoscalar,3> &b, std::array<geoscalar,3> &n);
template void getBAndN<geoscalar, std::array<geoscalar,3>>(const std::array<geoscalar,3> &f0, const std::array<geoscalar,3> &f1, const std::array<geoscalar,3> &f2, const std::array<geoscalar,3> &p3, const std::array<geoscalar,3> &t, std::array<geoscalar,3> &b, std::array<geoscalar,3> &n);
template geoscalar volumeForNode<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, int node);
template geoscalar volumeForNode<geoscalar>(const std::array<geoscalar,3> &tet0, const std::array<geoscalar,3> &tet1, const std::array<geoscalar,3> &tet2, const std::array<geoscalar,3> &tetN);
template void volumeAndAreaForNode<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, int node, geoscalar &volume, geoscalar &area);
template geoscalar volumeForIntPoint<geoscalar>(const std::array<geoscalar,3> &ip, const std::array<std::array<geoscalar,3>,4> &tetA, int e1, int e2, const std::array<std::array<geoscalar,3>,4> &tetB, int f1, int f2, int f3);
template geoscalar volumeForIntPointII<geoscalar>(const std::array<geoscalar,3> &ip, const std::array<geoscalar,3> &tetAe1, const std::array<geoscalar,3> &tetAe2, const std::array<geoscalar,3> &tetAe3, const std::array<geoscalar,3> &tetAe4, const std::array<geoscalar,3> &tetBf1, const std::array<geoscalar,3> &tetBf2, const std::array<geoscalar,3> &tetBf3, const std::array<geoscalar,3> &tetBf4);
template void volumeAndAreaForIntPoint<geoscalar>(const std::array<geoscalar,3> &ip, const std::array<std::array<geoscalar,3>,4> &tetA, int e1, int e2, const std::array<std::array<geoscalar,3>,4> &tetB, int f1, int f2, int f3, geoscalar &volume, geoscalar &area);
template geoscalar volumeIntersection<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, const std::array<std::array<geoscalar,3>,4> &tetB, bool calcCM, std::array<geoscalar,3> &cm);
template geoscalar volumeIntersectionII<geoscalar>(const std::array<geoscalar,3> &tetA0, const std::array<geoscalar,3> &tetA1, const std::array<geoscalar,3> &tetA2, const std::array<geoscalar,3> &tetA3, const std::array<geoscalar,3> &tetB0, const std::array<geoscalar,3> &tetB1, const std::array<geoscalar,3> &tetB2, const std::array<geoscalar,3> &tetB3, bool calcCM, std::array<geoscalar,3> &cm);

template void volumeAndAreaIntersection<geoscalar>(const std::array<std::array<geoscalar,3>,4> &tetA, const std::array<std::array<geoscalar,3>,4> &tetB, geoscalar &vol, geoscalar &area);
template void contribVolForNode<geoscalar>(geoscalar &vol, const std::array<geoscalar,3> &n0, const std::array<geoscalar,3> &n1, const std::array<geoscalar,3> &n2, const std::array<geoscalar,3> &n3, std::array<std::array<geoscalar,3>,56> &W, int &ips);
template void contribVolForIntPoint<geoscalar>(geoscalar &vol, const std::array<geoscalar,3> &ip, const std::array<geoscalar,3> &tetAi, const std::array<geoscalar,3> &tetAj, const std::array<geoscalar,3> &tetAe3, const std::array<geoscalar,3> &tetAe4, const std::array<geoscalar,3> &tetB0, const std::array<geoscalar,3> &tetB1, const std::array<geoscalar,3> &tetB2, const std::array<geoscalar,3> &tetB3, std::array<std::array<geoscalar,3>,56> &W, int &ips);
template void contribVolForIntersections<geoscalar>(geoscalar &vol, const std::array<geoscalar,3> &tetAi, const std::array<geoscalar,3> &tetAj, const std::array<geoscalar,3> &tetAe3, const std::array<geoscalar,3> &tetAe4, const std::array<geoscalar,3> &tetB0, const std::array<geoscalar,3> &tetB1, const std::array<geoscalar,3> &tetB2, const std::array<geoscalar,3> &tetB3, const std::array<std::array<geoscalar,3>,56> &W, int &ips);
#endif 
