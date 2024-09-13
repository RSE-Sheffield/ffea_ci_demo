/***********************************************************************\
 *
 * This file is part of the RngStreams package
 * 
 * File:           testRngStream.cpp for multiple streams of Random Numbers
 * Language:       C++ (ISO 1998)
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Date:           14 August 2001
 *
 * RngStreams is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * RngStreams is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RngStreams.  If not, see <http://www.gnu.org/licenses/>.
 *
 * We would appreciate that
 * if you use this software for work leading to publications, 
 * you cite the following relevant articles in which MRG32k3a 
 * and the package with multiple streams were proposed:

 * P. L'Ecuyer, ``Good Parameter Sets for Combined Multiple Recursive Random Number Generators'', 
 * Operations Research, 47, 1 (1999), 159--164.

 * P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton, 
 * ``An Objected-Oriented Random-Number Package with Many Long Streams and Substreams'', 
 * Operations Research, 50, 6 (2002), 1073--1075
 *
\***********************************************************************/

//  Program to test the random number streams file:    RngStream.cpp

#include <iostream>
#include <sstream>
#include <array>
#include "RngStream.h"
using namespace std;

int main ()
{
   double sum;
   int  i;
   RngStream g1 ("g1");
   RngStream g2 ("g2");
   RngStream g3 ("g3");

   sum = g2.RandU01 () + g3.RandU01 ();

   g1.AdvanceState (5, 3);   
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   for (i = 0;  i < 35; i++)
      g1.AdvanceState (0,1);
   sum += g1.RandU01 ();

   g1.ResetStartStream ();
   long sumi = 0;
   for (i = 0;  i < 35; i++)
      sumi += g1.RandInt (1, 10);
   sum += sumi / 100.0;

   double sum3 = 0.0;
   for (i = 0;  i < 100;  i++) {
      sum3 += g3.RandU01 ();
   }
   sum += sum3 / 10.0;

   g3.ResetStartStream ();
   for (i=1; i<=5; i++)
      sum += g3.RandU01 ();

   for (i=0; i<4; i++)
      g3.ResetNextSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g3.ResetStartSubstream ();
   for (i=0; i<5; i++)
      sum += g3.RandU01 ();

   g2.ResetNextSubstream ();
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g2.RandU01 ();
   sum += sum3 / 10000.0;

   g3.SetAntithetic (true);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g3.RandU01 ();
   sum += sum3 / 10000.0;

   std::array<uint32_t, 6> germe = { 1, 1, 1, 1, 1, 1 };
   RngStream::SetPackageSeed (germe.data());

   RngStream gar[4] = { "Poisson", "Laplace", "Galois", "Cantor" };
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   gar[2].AdvanceState (-127, 0);
   sum += gar[2].RandU01 ();

   gar[2].ResetNextSubstream ();
   gar[2].IncreasedPrecis (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;

   gar[2].SetAntithetic (true);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += gar[2].RandU01 ();
   sum += sum3 / 10000.0;
   gar[2].SetAntithetic (false);

   gar[2].IncreasedPrecis (false);
   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   cout.precision(14);
   string result = "39.697547445251";
   stringstream ss; 
   ss.precision(14);
   ss << sum;
   if ( result.compare(ss.str()) ) {
     cout << "wrong result: " << sum << " " << ss.str() << endl;
     return 1;
   } else {
     // cout << "right result: " << sum << " " << ss.str() << endl;
     return 0; 
   }

   //   cin >> i; // for Visual C++
   return 0;
}
