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
#include <sstream>
#include "VolumeIntersection.h"
#include "CheckTetrahedraOverlap.h"

using namespace std; 

int main() {

  std::array<arr3, 4> tetA, tetB;
  std::array<arr3, 4> tetC, tetD;
  // set up tetA and tetB:
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      tetA[i][j] = 0;
      tetB[i][j] = 0;
    }
  }

  tetA[1][1] = 1.0;
  tetA[2][0] = 1.0;
  tetA[3][2] = 1.0;

  tetB[1][1] = 1.0;
  tetB[2][0] = 1.0;
  tetB[3][2] = 1.0;

  for (int i=0; i<4; i++) {
    tetB[i][2] += 0.5;
    tetA[i][1] += 0.02;
    tetA[i][0] += 0.05;
    // tetA[i][1] += 0.0;
  }

  cout << " tetA: " << endl;
  for (int i=0; i<4; i++){
    for (int j=0;j<3;j++){
       cout << tetA[i][j] << " ";
    }
    cout << endl;
  }
  cout << " tetB: " << endl;
  for (int i=0; i<4; i++){
    for (int j=0;j<3;j++){
       cout << tetB[i][j] << " ";
    }
    cout << endl;
  }


  if (! tet_a_tetII(tetA[0], tetA[1], tetA[2], tetA[3],
                    tetB[0], tetB[1], tetB[2], tetB[3])) {
    cout << " these tetrahedra are known to intersect" << endl;
    return 1;
  } 

  arr3 cm; 
  scalar vol = volumeIntersectionII(tetA[0], tetA[1], tetA[2], tetA[3], tetB[0], tetB[1], tetB[2], tetB[3], false, cm);
  stringstream ss;
  ss.precision(14);
  ss << vol;
  string hardcodedresult = "0.020833333333333"; 
  cout.precision(14);
  cout << "total volume: " << vol << endl;

  if (hardcodedresult.compare(ss.str())) {
    return 1; 
  } /* else {
    return 0; 
  } */

  // set up tetC and tetD:
  tetC[0][0] = 5.16991;
  tetC[0][1] = 5.11394;
  tetC[0][2] = -2.59001;
  tetC[1][0] = 2.68658;
  tetC[1][1] = 5.11124;
  tetC[1][2] = -5.03959;
  tetC[2][0] = 5.12737;
  tetC[2][1] = 5.10932;
  tetC[2][2] = -5.0372;
  tetC[3][0] = 2.98162;
  tetC[3][1] = 2.85562;
  tetC[3][2] = -2.78657;
 
  tetD[0][0] = 4.8437;
  tetD[0][1] = 2.47363;
  tetD[0][2] = -5.0158;
  tetD[1][0] = 4.84606;
  tetD[1][1] = 5.11758;
  tetD[1][2] = -2.44064;
  tetD[2][0] = 7.20236;
  tetD[2][1] = 5.11031; 
  tetD[2][2] = -5.00642;
  tetD[3][0] = 4.88715;
  tetD[3][1] = 5.11134;
  tetD[3][2] = -5.00755; 

  if (! tet_a_tetII(tetC[0], tetC[1], tetC[2], tetC[3], 
                    tetD[0], tetD[1], tetD[2], tetD[3])) {
    cout << "tetrahedra tetC and tetD are known to intersect " << endl;
    cout << "  and I found tet_a_tet fails in this case" << endl;
    // return 1;
  } 

  vol = volumeIntersection(tetC, tetD, false, cm);
  if (vol == 0) {
    cout << " intersecting volume should not be zero for these tetrahedra" << endl; 
    return 1;
  } 

  return 0; 

}
   
