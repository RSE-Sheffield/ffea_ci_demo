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

#ifndef MAT_VEC_FNS_II_H_INCLUDED
#define MAT_VEC_FNS_II_H_INCLUDED

#include <cmath>
#include <cstring>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////
//
///  check whether two scalars have the same sign 
template <class scalar>
bool sameSign(scalar a, scalar b) {
    // return a <= ffea_const::zero == b <= ffea_const::zero;
    if (a*b >= 0) return true;
    return false;
}

/** Given 3 integers n0, n1, n2
 *    return the index missing in the list [0,1,2,3] 
 */
inline int getMissingNode(int n0, int n1, int n2) {
  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        return 3;
      } else {
        return 2;
      }
    } else {
      return 1;
    }
  } else {
    return 0;
  }
}

/** Given 1 integers iN,
 *    return the index missing indices of the list [0,1,2,3] 
 */
inline void getRestOfNodes(int iN, int &iO0, int &iO1, int &iO2) {
  if (iN == 0) {
    iO0 = 1; 
    iO1 = 2; 
    iO2 = 3; 
  } else if (iN == 1) {
    iO0 = 0; 
    iO1 = 2; 
    iO2 = 3; 
  } else if (iN == 2) {
    iO0 = 0; 
    iO1 = 1; 
    iO2 = 3; 
  } else {
    iO0 = 0; 
    iO1 = 1; 
    iO2 = 2; 
  }
}

/** Given 1 integers iN,
 *    return the index missing indices of the list [0,1,2,3] 
 */
inline void getMissingPair(int in0, int in1, int &on0, int &on1) {
  for (int i=0; i<4; i++) {
    if ((i != in0) && (i != in1)) {
      on0 = i;
      on1 = getMissingNode(in0, in1, on0);
      break; 
    }
  }
}


///////////////// SECTION 1 ////////////////////
////  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////

/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
void add(const Array1<T, N>&vecA, const Array2<T, N> &vecB, Array3<T, N>&res) {
    for (int i = 0; i < N; ++i) {
        res[i] = vecA[i] + vecB[i];
    }
}
template <typename T>
void add(const std::vector<T> &vecA, const std::vector<T> &vecB, std::vector<T> &res) {
    for (int i = 0; i < res.size(); ++i) {
        res[i] = vecA[i] + vecB[i];
    }
}


/** res = vecA - vecB */
/** res can be either vecA or vecB */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
void sub(const Array1<T, N> &vecA, const Array2<T, N>&vecB, Array3<T, N>&res) {
    for (int i = 0; i < N; ++i) {
        res[i] = vecA[i] - vecB[i];
    }
}
template <typename T>
void sub(const std::vector<T> &vecA, const std::vector<T> &vecB, std::vector<T> &res) {
    for (int i = 0; i < res.size(); ++i) {
        res[i] = vecA[i] - vecB[i];
    }
}

/** res = vecA - vecB */
/** res can be either vecA or vecB */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
T dot(const Array1<T,N> &vecA, const Array2<T, N> &vecB) {
    T result = 0;
    for (int i = 0; i < N; ++i) {
        result += vecA[i] * vecB[i];
    }
    return result;
}
template <typename T>
T dot(const std::vector<T> &vecA, const std::vector<T> &vecB) {
    T result = 0;
    for (int i = 0; i < vecA.size(); ++i) {
        result += vecA[i] * vecB[i];
    }
    return result;
}
/** Normalise vector e into result*/
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
void normalize(const Array1<T, N> &e, Array2<T, N> &result) {
    T sum = 0;
    for (int i = 0; i < N; ++i) {
        sum += e[i] * e[i];
    }
    if (sum == 0) throw FFEAException("Unable to normalize a vector with length zero");
    sum = std::sqrt(sum);
    for (int i = 0; i < N; ++i) {
        result[i] = e[i] / sum;
    }
}
template <typename T>
void normalize(const std::vector<T> &e, std::vector<T> &res) {
    T sum = 0;
    for (int i = 0; i < e.size(); ++i) {
        sum += e[i] * e[i];
    }
    if (sum == 0) throw FFEAException("Unable to normalize a vector with length zero");
    sum = std::sqrt(sum);
    for (int i = 0; i < e.size(); ++i) {
        res[i] = e[i] / sum;
    }
}
template <typename T, size_t N,
    template<typename, size_t> typename Array1>
void normalize(Array1<T, N>& e) {
    normalize(e, e);
}
/** resize vector u into vector v, given scalar f v*/
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
void resize2(T f, const Array1<T, N> &u, Array2<T, N> &v) {
    for (int i = 0; i < N; ++i) {
        v[i] = f * u[i];
    }
}
template <typename T, size_t N,
    template<typename, size_t> typename Array1>
void resize(T f, Array1<T, N>&u) {
    resize2(f, u, u);
}
template <typename T>
void resize2(T f, const std::vector<T> &u, std::vector<T> &v) {
    for (int i = 0; i < u.size(); ++i) {
        v[i] = f * u[i];
    }
}
/** Given a scalar f, change v so that v += f*u */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
void resize3(T f, const Array1<T, N> &u, Array2<T, N> &v) {
    for (int i = 0; i < N; ++i) {
        v[i] += f * u[i];
    }
}
template <typename T>
void resize3(T f, const std::vector<T>& u, std::vector<T> &v) {
    for (int i = 0; i < u.size(); ++i) {
        v[i] += f * u[i];
    }
}
// Extra version for arr_view created in-place
template <typename T, size_t N,
    template<typename, size_t> typename Array1>
void resize3(T f, const Array1<T, N>& u, arr_view<T, N> v) {
    for (int i = 0; i < N; ++i) {
        v[i] += f * u[i];
    }
}
/** copy arr3 src into arr3 dest
 *
 *  @warning For backwards-compatibility purposes, the src/dest arg order is reverse to standard C++ memcpy etc
 */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
void store(const Array2<T, N> &src, Array1<T, N> &dest) {
    memcpy(dest.data(), src.data(), sizeof(T) * N);
}
template <typename T>
void store(const std::vector<T> &src, std::vector<T> &dest) {
    memcpy(dest.data(), src.data(), sizeof(T) * dest.size());
}
// Extra version for arr_view created in-place
template <typename T, size_t N,
    template<typename, size_t> typename Array2>
void store(const Array2<T, N>& src, arr_view<T, N> dest) {
    memcpy(dest.data(), src.data(), sizeof(T) * N);
}
/** Given a scalar f, change v so that v += f*u */
template <typename T, size_t N,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
T distance(const Array1<T, N> &a, const Array2<T, N> &b) {
    T result = 0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < N; ++i) {
        const T t = a[i] - b[i];
        result += t * t;
    }
    return std::sqrt(result);
}
template <typename T>
T distance(const std::vector<T> &a, const std::vector<T> &b) {
    T result = 0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < a.size(); ++i) {
        const T t = a[i] - b[i];
        result += t * t;
    }
    return std::sqrt(result);
}
/** Return the length of a vector */
template <typename T, size_t N, template<typename, size_t> typename Array>
T magnitude(const Array<T, N>& a) {
    T result = 0.0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < N; i++) {
        result += a[i] * a[i];
    }
    return std::sqrt(result);
}
template <typename T>
T magnitude(const std::vector<T>& a) {
    T result = 0.0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < a.size(); i++) {
        result += a[i] * a[i];
    }
    return  std::sqrt(result);
}
/** Given a scalar f, change v so that v += f*u */
template <typename T, size_t N, template<typename, size_t> typename Array>
T magnitude2(const Array<T, N> &a) {
    T result = 0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < N; ++i) {
        const T t = a[i] - a[i];
        result += t * t;
    }
    return result;
}
template <typename T>
T magnitude2(const std::vector<T> &a) {
    T result = 0;
#pragma omp simd reduction(+:result)
    for (int i = 0; i < a.size(); ++i) {
        const T t = a[i] - a[i];
        result += t * t;
    }
    return result;
}
/** Fill the iterable with zeros */
template <typename T, size_t N, template<typename, size_t> typename Array>
void initialise(Array<T, N> &a) {
    memset(a.data(), 0, sizeof(T) * N);
}
template <typename T>
void initialise(std::vector<T> &a) {
    memset(a.data(), 0, sizeof(T) * a.size());
}
/** w = u x v */
/** w must be different from u and v */
// void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w));
// template <class brr3> void arr3arr3VectorProduct(brr3 (&u), brr3 (&v), brr3 (&w));
template<
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
void arr3arr3VectorProduct(const Array1<scalar, 3> &u, const Array2<scalar, 3> &v, Array3<scalar, 3> &w) {
    w[0] = u[1] * v[2] - v[1] * u[2];
    w[1] = -u[0] * v[2] + v[0] * u[2];
    w[2] = u[0] * v[1] - v[0] * u[1];
}

template<
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
scalar detByRows(const Array1<scalar, 3> &a, const Array2<scalar, 3> &b, const Array3<scalar, 3> &c){
  scalar det;
  det  = a[0] * (b[1] * c[2] - b[2] * c[1]);
  det += a[1] * (b[2] * c[0] - b[0] * c[2]);
  det += a[2] * (b[0] * c[1] - b[1] * c[0]);
  return det;
}
template<
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
scalar detByCols(const Array1<scalar, 3> &a, const Array2<scalar, 3> &b, const Array3<scalar, 3> &c){
  scalar det;
  det  = a[0] * (b[1] * c[2] - b[2] * c[1]);
  det += b[0] * (c[1] * a[2] - c[2] * a[1]);
  det += c[0] * (a[1] * b[2] - a[2] * b[1]);
  return det; 
}


///////////////// SECTION 2 ////////////////////
///// Geometric  functions for arr3 types ////// 
////////////////////////////////////////////////
/** t = unit(vecA - vecB) */
/** t can be either vecA or vecB */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
void tangent(const Array1<T, 3> &vecA, const Array2<T, 3> &vecB, Array3<T, 3> &t) {
    T w = 0;
    for (int i = 0; i < 3; i++) {
        t[i] = vecA[i] - vecB[i];
        w += t[i] * t[i];
    }
    w = sqrt(w);
    for (int i = 0; i < 3; i++) {
        t[i] /= w;
    }
}

/** w = unit(u x v) */
/**  (w != u) && (w != v) */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
void getUnitNormal(const Array1<T, 3> &u, const Array2<T, 3> &v, Array3<T, 3> &w) {
    w[0] = u[1] * v[2] - v[1] * u[2];
    w[1] = -u[0] * v[2] + v[0] * u[2];
    w[2] = u[0] * v[1] - v[0] * u[1];
    T l;
    l = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    for (int i = 0; i < 3; i++) {
        w[i] /= l;
    }
}

/** calculate the unit normal vector n to the plane defined by the three points */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
void getNormal(const Array1<T, 3> &v1, const Array2<T, 3> &v2, const Array3<T, 3> &v3, Array4<T, 3> &n) {
    arr3 pl1, pl2;
    sub(v2, v1, pl1);
    sub(v3, v1, pl2);
    getUnitNormal(pl1, pl2, n);
}

/** Given the face formed by tetA[0]:tetA[1]:tetA[2] 
 * get n, the normal to a face pointing inwards.
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
void getNormalInwards(const std::array<Array1<T, 3>, 4> &tetA, int n0, int n1, int n2, Array2<T, 3> &n) {
    // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
    arr3 pl1, pl2, aux;
    sub(tetA[n1], tetA[n0], pl1);
    sub(tetA[n2], tetA[n0], pl2);
    // n is a unit vector normal to the face:
    arr3arr3VectorProduct(pl1, pl2, aux);
    normalize(aux, n);
    // but it must be inwards, i. e., on the same side of the plane than n3. 
    int n3 = getMissingNode(n0, n1, n2);
    add(aux, tetA[n0], aux);
    scalar d = - dot(n, tetA[n0]);
    scalar t1 = dot(tetA[n3], n) + d;
    scalar t2 = dot(aux, n) + d;
    if (!sameSign(t1,t2)) resize2(ffea_const::mOne, n, n);
}
/** Given the face formed by f0, f1, and f2,
 *     and knowing the remaining p3 for a tetrahedron,
 * get n, the normal to a face pointing inwards.
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5>
void getNormalInwards(const Array1<T, 3> &f0, const Array2<T, 3> &f1, const Array3<T, 3> &f2,
    const Array4<T, 3> &p3, Array5<T, 3> &n) {
    // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
    arr3 pl1, pl2, aux;
    sub(f1, f0, pl1);
    sub(f2, f0, pl2);
    // n is a unit vector normal to the face:
    arr3arr3VectorProduct(pl1, pl2, aux);
    normalize(aux, n);
    // but it must be inwards, i. e., on the same side of the plane than n3.
    add(aux, f0, aux);
    scalar d = -dot(n, f0);
    scalar t1 = dot(p3, n) + d;
    scalar t2 = dot(aux, n) + d;
    if (!sameSign(t1, t2)) resize2(ffea_const::mOne, n, n);
}

/** check if points vec and test are at the same side
 *  of the plane formed by p1, p2 and p3 
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5>
bool sameSidePlane(const Array1<T, 3> &vec, const Array2<T, 3> &test,
    const Array3<T, 3> &p1, const Array4<T, 3> &p2, const Array5<T, 3> &p3) {
    arr3 pl1, pl2, n;
    sub(p2, p1, pl1);
    sub(p3, p1, pl2);
    arr3arr3VectorProduct(pl1, pl2, n);
    normalize(n);
    scalar d = -dot(n, p1);
    scalar t1 = dot(vec, n) + d;
    scalar t2 = dot(test, n) + d;
    // if ((t1 * t2) >= 0 ) return true;
    // cout << "t1: " << t1 << " t2: " << t2 << endl; 
    if (sameSign(t1, t2)) return true;
    return false;    
}

/**  Given 4 co-planar points, check if ip and p1 lay on the same side 
 *     of the of the line formed by p2 and p3.
 *
 *   More specifically we check whether pl23 x pl21 and pl23 x pl2e 
 *     are parallel or antiparallel. 
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
bool sameSideLine(const Array1<T, 3> &e, const Array2<T, 3> &p1, const Array3<T, 3> &p2, const Array4<T, 3> &p3) {
    // if (!samePlane(e, p1, p2, p3)) cout << "alarm" << endl;
    arr3 pl23, pl21, pl2e, v1, ve;
    sub(p3, p2, pl23);
    sub(p1, p2, pl21);
    sub(e, p2, pl2e);
    arr3arr3VectorProduct(pl21, pl23, v1);
    arr3arr3VectorProduct(pl2e, pl23, ve);
    if (dot(v1,ve) >= 0) return true;
    return false;
}

/**  check whether vector vec is in tetrahedron B.
 *
 *  more specifically, it will be there if 
 *     for each plane of the tetrahedron, 
 *     the point is on the same side as the remaining vertex */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5>
bool nodeInTet(const Array1<T, 3> &vec,
    const Array2<T, 3> &tet0, const Array3<T, 3> &tet1, const Array4<T, 3> &tet2, const Array5<T, 3> &tet3) {
    if (!sameSidePlane(vec, tet0, tet1, tet2, tet3)) return false;
    if (!sameSidePlane(vec, tet1, tet2, tet3, tet0)) return false;
    if (!sameSidePlane(vec, tet2, tet3, tet0, tet1)) return false;
    if (!sameSidePlane(vec, tet3, tet0, tet1, tet2)) return false;

    return true;
}
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2>
bool nodeInTet(const Array1<T, 3>& vec, const std::array<Array2<T, 3>, 4>& tet) {
    return nodeInTet(vec, tet[0], tet[1], tet[2], tet[3]);
}

/** Find the intersection point of the line that passes through the points e1 and e2, 
 *   and the plane defined by points p1, p2 and p3.
 *  \warning This function should be called ONLY in the case 
 *            that intersection is known to occur. 
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5,
    template<typename, size_t> typename Array6>
void linePlaneIntersectionPoint(Array1<T,3> &ip,
    const Array2<T, 3> &e1, const Array3<T, 3> &e2,
    const Array4<T, 3> &p1, const Array5<T, 3> &p2, const Array6<T, 3> &p3) {
   // v is the vector of the line L(t) = e1 + v*t
   arr3 v;
   sub(e2, e1, v);

   // now we need the unit vector that defines the plane, pn:
   arr3 pl1, pl2, pn;
   sub(p2, p1, pl1);
   sub(p3, p2, pl2);
   getUnitNormal(pl1, pl2, pn);

   // the plane is defined through: ax + by + cz + d = 0; 
   //   (a,b,c) = pn
   // so to find d we simply:
   scalar d = - dot(pn, p1);

   // now find t and the point:
   scalar t = - (dot(pn, e1) + d) / dot(pn, v);
   resize(t, v);
   add(e1, v, ip);
}

/** Return true and 
 *         the intersection point of the line that passes through the points e1 and e2
 *           and the plane defined by points p1, p2, p3 if this intersection actually occurs,
 *         and false otherwise. 
 *         \warning p1, p2, and p3 are assumed to be non-colinear.
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5,
    template<typename, size_t> typename Array6>
bool safeLinePlaneIntersectionPoint(Array1<T,3> &ip,
    const Array2<T, 3> &e1, const Array3<T, 3> &e2,
    const Array4<T, 3> &p1, const Array5<T, 3> &p2, const Array6<T, 3> &p3) {
    // v is the vector of the line L(t) = e1 + v*t
    arr3 v;
    sub(e2, e1, v);

    // now we need the unit vector that defines the plane, pn:
    arr3 pl1, pl2, pn;
    sub(p2, p1, pl1);
    sub(p3, p2, pl2);
    getUnitNormal(pl1, pl2, pn);

    // CHECK! 
    scalar pnv = dot(pn, v);
    if (abs(pnv) < ffea_const::threeErr) return false; // the line won't intersect the plane!

    // the plane is defined through: ax + by + cz + d = 0; 
    //   (a,b,c) = pn
    // so to find d we simply:
    scalar d = -dot(pn, p1);

    // now find t and the point:
    scalar t = -(dot(pn, e1) + d) / dot(pn, v);
    resize(t, v);
    add(e1, v, ip);

    return true;
}

/** Return true and 
 *         the intersection point ip of the line that passes through the points e1 and e2
 *           and face defined by points p1, p2, p3 if this intersection actually occurs,
 *         and false otherwise. 
 *         \warning p1, p2, and p3 are assumed to be non-colinear.
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5,
    template<typename, size_t> typename Array6>
bool lineFaceIntersectionPoint(Array1<T,3> &ip,
    const Array2<T, 3> &e1, const Array3<T, 3> &e2,
    const Array4<T, 3> &p1, const Array5<T, 3> &p2, const Array6<T, 3> &p3) {
    // look for the intersection point... if it exists 
    if (!safeLinePlaneIntersectionPoint(ip, e1, e2, p1, p2, p3)) return false;

    // and finally check whether this point ip belongs to the triangular face:
    if ((isPointInFace(ip, p1, p2, p3))) return true;

    return false;
}


/** Check whether point ip is  
  *    inside of the three half-planes formed by the triangle's edges p1, p2, p3.
  */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
bool isPointInFace(Array1<T,3> &ip, const Array2<T,3> &p1, const Array3<T,3> &p2, const Array4<T,3> &p3) {
    if (!sameSideLine(ip, p1, p2, p3)) return false;
    if (!sameSideLine(ip, p3, p1, p2)) return false;
    if (!sameSideLine(ip, p2, p3, p1)) return false;
    return true;
}

/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * More specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (tet[f2] - tet[f1]) and (tet[f3] - tet[f1]).
 *    - the intersection of a line is a point in the plane 
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
bool intersectionPoint(Array1<T,3> &ip, const Array2<T,3> &e1, const Array3<T,3> &e2,
    const std::array<Array4<T, 3>, 4> &tet, int f1, int f2, int f3) {
    // check it e1 and e2 are on the same side of the plane:
    if (sameSidePlane(e1, e2, tet[f1], tet[f2], tet[f3])) return false;

    // given that they are on different sides of the plane look for the intersection point.
    linePlaneIntersectionPoint(ip, e1, e2, tet[f1], tet[f2], tet[f3]);

    // and finally check whether this point ip belongs to the triangular face:
    if (isPointInFace(ip, tet[f1], tet[f2], tet[f3])) return true;
    return false;
}
/** Check whether an edge and a plane intersect,
 *    and return the intersection point ip and true if found, false otherwise.
 * more specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (f2 - f1) and (f3 - f1).
 *    - the intersection of a line is a point in the plane
 */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5,
    template<typename, size_t> typename Array6>
bool intersectionPoint(Array1<T,3> &ip,
    const Array2<T,3> &e1, const Array3<T,3> &e2,
    const Array4<T,3> &f1, const Array5<T,3> &f2, const Array6<T,3> &f3) {
    // check it e1 and e2 are on the same side of the plane:
    if (sameSidePlane(e1, e2, f1, f2, f3)) return false;

    // given that they are on different sides of the plane look for the intersection point.
    linePlaneIntersectionPoint(ip, e1, e2, f1, f2, f3);

    // and finally check whether this point ip belongs to the triangular face:
    if ((isPointInFace(ip, f1, f2, f3))) return true;    
    return false;
}

/** Given a line defined by point p1 and vector p2p1, get the intersecting point p3,
 *   in that line from a third point p0.
 * Essentially implementing "Intersection of two lines in three-space",
 *  by Ronald Goldman, in Graphics Gems I. */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
void intersectingPointToLine(const Array1<T,3> &p0, const Array2<T,3> &p1, const Array3<T,3> &p2p1, Array4<T,3> &p3) {
    arr3 p0p1, cp3p0, v2, v1v2;
    scalar c;

    sub(p0, p1, p0p1);
    arr3arr3VectorProduct(p0p1, p2p1, v2);
    // CHECK! if the cross product is null, then distance is zero, and p3 is p0 itself.
    if (magnitude(v2) < ffea_const::threeErr) {
        store(p0, p3);
        return;
    }
    // get cp3p0, the vector to (or from) to the line p2p1:
    arr3arr3VectorProduct(p2p1, v2, cp3p0);

    // now calculate c, the amount of cp3p0 from p0 to the intersection point p3.
    arr3arr3VectorProduct(p2p1, cp3p0, v1v2);
    c = detByRows(p0p1, p2p1, v1v2);
    c /= (v1v2[0] * v1v2[0] + v1v2[1] * v1v2[1] + v1v2[2] * v1v2[2]);

    // and finally calculate the intersection point: 
    resize(c, cp3p0);
    add(p0, cp3p0, p3);

    /*
    // CHECKS START HERE! //
    arr3 p3p0;
    sub(p3, p0, p3p0);
    c = mag<scalar,arr3>(p3p0);

    // check the distance:
    scalar d = distanceFromPointToLine(p0, p1, p2);
    if (abs (d - abs(c)) > ffea_const::threeErr) {
      cout << "something is wrong: d = " << d << " differs from c = " << c << endl;
    }
    */
}

/** Given a line defined by points p1 and p2, 
 *    return the distance from p0, to this line. */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3>
T distanceFromPointToLine(const Array1<T,3> &p0, const Array2<T,3> &p1, const Array3<T,3> &p2) {
   arr3 p2p1, p0p1, p0p2, tmp1;
   sub(p2, p1, p2p1);
   sub(p0, p2, p0p2);
   sub(p0, p1, p0p1);
   arr3arr3VectorProduct(p0p1, p0p2, tmp1);
   scalar d = magnitude(tmp1) / magnitude(p2p1);
   return d; 
}

template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4>
T getTetrahedraVolume(const Array1<T,3> &p0, const Array2<T,3> &p1, const Array3<T,3> &p2, const Array4<T,3> &p3){
  return  (p1[0] - p0[0])*( (p1[1] - p2[1])*(p2[2] - p3[2]) - (p2[1] - p3[1])*(p1[2] - p2[2]) ) +
          (p2[0] - p1[0])*( (p2[1] - p3[1])*(p0[2] - p1[2]) - (p0[1] - p1[1])*(p2[2] - p3[2]) ) +
          (p3[0] - p2[0])*( (p0[1] - p1[1])*(p1[2] - p2[2]) - (p1[1] - p2[1])*(p0[2] - p1[2]) );
  
}

/** Return the center of coordinates for four points p1, p2, p3, p4 in c */
template <typename T,
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5>
void getTetrahedraCM(
    const Array1<T, 3>& p1,
    const Array2<T, 3>& p2,
    const Array3<T, 3>& p3,
    const Array4<T, 3>& p4,
    Array5<T, 3>& c) {
    for (int i = 0; i < 3; ++i) {
        c[i] = 0.25 * (p1[i] + p2[i] + p3[i] + p4[i]);
    }
}


template <
    template<typename, size_t> typename Array1,
    template<typename, size_t> typename Array2,
    template<typename, size_t> typename Array3,
    template<typename, size_t> typename Array4,
    template<typename, size_t> typename Array5,
    template<typename, size_t> typename Array6>
void getLocalCoordinatesForLinTet(
    const Array1<scalar, 3> &t0,
    const Array2<scalar, 3> &t1,
    const Array3<scalar, 3> &t2,
    const Array4<scalar, 3> &t3,
    const Array5<scalar, 3> &p,
    Array6<scalar, 4> &phi) {
    phi[0] = getTetrahedraVolume(p, t1, t2, t3);
    phi[1] = -getTetrahedraVolume(p, t0, t2, t3);
    phi[2] = getTetrahedraVolume(p, t0, t1, t3);
    phi[3] = -getTetrahedraVolume(p, t0, t1, t2);
    scalar v = getTetrahedraVolume(t0, t1, t2, t3);
    for (int i = 0; i < 4; i++) {
        phi[i] /= v;
    }

    /*  // CHECK!! for the time:
    arr3 tmp;
    arr3 pc = {0,0,0}; // WT*?
    resize2(phi[0], t0, tmp);
    add(pc, tmp, pc);
    resize2(phi[1], t1, tmp);
    add(pc, tmp, pc);
    resize2(phi[2], t2, tmp);
    add(pc, tmp, pc);
    resize2(phi[3], t3, tmp);
    add(pc, tmp, pc);
    sub(p, pc, tmp);
    if (magnitude(tmp) > 1e-6) {
        cout << "local coordinates were not correctly calculated!" << endl;
        cout << "p: " << p[0] << ", " << p[1] << ", " << p[2] << endl;
        cout << "pc: " << pc[0] << ", " << pc[1] << ", " << pc[2] << endl;
        cout << "diff: " << magnitude(tmp) << endl;
    } */
}


#endif 
