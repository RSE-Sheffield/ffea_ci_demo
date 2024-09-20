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

#include "mat_vec_fns_II.h"
#include <iostream>
using std::cout;
using std::endl;


////////////////// SECTION 1 ///////////////////////
///  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////
/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
//template <class scalar, class arr3>
//void arr3arr3Add(arr_view<scalar,3> vecA, arr_view<scalar,3> vecB, arr_view<scalar,3> res){
//
//   for (int i=0; i<3; i++) {
//     res[i] = vecA[i] + vecB[i];
//    }
//
//}

///** res = vecA - vecB */
///** res can be either vecA or vecB */
//// void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res){
//// template <class arr3> void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res){
//template <class scalar, class arr3>
//void arr3arr3Substract(arr_view<scalar,3> vecA, arr_view<scalar,3> vecB, arr_view<scalar,3> res){
//
//   for (int i=0; i<3; i++) {
//     res[i] = vecA[i] - vecB[i];
//    }
//
//}
///** w = u x v */
///**  (w != u) && (w != v) */
//// void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w)){
//template <class scalar, class arr3>
//void arr3arr3VectorProduct(arr_view<scalar,3> u, arr_view<scalar,3> v, arr_view<scalar,3> w){
//
//    w[0] = u[1]*v[2] - v[1]*u[2];
//    w[1] = -u[0]*v[2] + v[0]*u[2];
//    w[2] = u[0]*v[1] - v[0]*u[1];
//
//
//}

///** return the dot product for arrays vecA and vecB */
//// scalar arr3arr3DotProduct(arr3 &vecA, arr3 &vecB) {
//template <class scalar, class arr3>
//scalar arr3arr3DotProduct(arr_view<scalar,3> vecA, arr_view<scalar,3> vecB) {
//
//  scalar result = 0.0;
//  for (int i=0; i<3; i++) {
//     // cout << "vecA[" << i << "]: " << vecA[i] << " vecB[" << i << "]: " << vecB[i] << endl; 
//     result += vecA[i] * vecB[i];
//  }
//  return result;
//}
//
///** Normalise vector arr3 e */
//template <class scalar, class arr3>
//void arr3Normalise(arr_view<scalar,3> e){
//
//   scalar norm = 0.0;
//   for (int i=0; i<3; i++) {
//     norm += e[i]*e[i];
//   }
//   norm = sqrt(norm);
//   for (int i=0; i<3; i++) {
//     e[i] /= norm;
//   }
//
//}
//
///** get the normalised vector of arr3 e into arr3 n */
//// template <class scalar, class arr3> void arr3Normalise2(arr3 &e, arr3 &n){
//template <class scalar, class arr3>
//void arr3Normalise2(arr_view<scalar,3> e, arr_view<scalar,3> n){
//
//   scalar norm = 0.0;
//   for (int i=0; i<3; i++) {
//     norm += e[i]*e[i];
//   }
//   if (norm == 0.0) throw -1; 
//   norm = sqrt(norm);
//   for (int i=0; i<3; i++) {
//     n[i] = e[i]/norm;
//   }
//
//}
//
///** resize vector u, given scalar f */
//template <class scalar, class arr3>
//void arr3Resize(scalar f, arr_view<scalar,3> u){
//
//   for (int i=0; i<3; i++) {
//      u[i] = f*u[i];
//    }
//
//}
//
///** resize vector u into vector v, given scalar f */
//template <class scalar, class arr3>
//void arr3Resize2(scalar f, arr_view<scalar,3> u, arr_view<scalar,3> v){
//
//   #pragma omp simd
//   for (int i=0; i<3; i++) {
//      v[i] = f*u[i];
//    }
//
//}
//
///** Given a scalar f, change v so that v += f*u */
//template <class scalar, class arr3>
//void arr3Resize3(scalar f, arr_view<scalar,3> u, arr_view<scalar,3> v){
//
//   #pragma omp simd
//   for (int i=0; i<3; i++) {
//      v[i] += f*u[i]; 
//   }
//}
//
//
///** cp arr3 u into arr3 v */ 
//template <class scalar, class arr3>
//void arr3Store(arr_view<scalar,3> u, arr_view<scalar,3> v){
//
//   #pragma omp simd
//   for (int i=0; i<3; i++) {
//     v[i] = u[i];
//    }
//
//}
//
///** return the distance from vecA to vecB */
//template <class scalar, class arr3>
//scalar arr3arr3Distance(arr_view<scalar,3> vecA, arr_view<scalar,3> vecB){
// 
//    scalar d=0.0;
//    for (int i=0; i<3; i++){ 
//      d  += (vecA[i] - vecB[i])*(vecA[i] - vecB[i]); 
//    } 
//    return sqrt(d); 
//}
//
//
///** Return the length of a vector v */
//template <class scalar, class arr3>
//scalar mag(arr_view<scalar,3> v) {
//
//   scalar s=0.0;
//   #pragma omp simd reduction(+:s)
//   for (int i=0; i<3; i++) {
//      s += v[i] * v[i];
//   }
//   return sqrt(s);
//
//}
//
//
///** Return the squared length of a vector v */
//template <class scalar, class arr3>
//scalar mag2(arr_view<scalar,3> v) {
//
//   scalar s=0.0;
//   for (int i=0; i<3; i++) {
//      s += v[i] * v[i];
//   }
//   return s;
//
//}
//
//
//template <class arr3>
//void arr3Initialise(arr3 &v){
//
//    for (int i=0; i<3; i++) {
//      v[i] = ffea_const::zero;
//    }
//
//}





////////////////////////////////////////////////
////////////// END OF SECTION 1 ////////////////
////////////////////////////////////////////////





////////////////////////////////////////////////
////////////// END OF SECTION 2 ////////////////
////////////////////////////////////////////////
