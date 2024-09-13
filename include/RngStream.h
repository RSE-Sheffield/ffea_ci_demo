/***********************************************************************\
 *
 * This file is part of the RngStreams package
 * 
 * File:           RngStream.h for multiple streams of Random Numbers
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
 * f you use this software for work leading to publications, 
 * you cite the following relevant articles in which MRG32k3a 
 * and the package with multiple streams were proposed:

 * P. L'Ecuyer, ``Good Parameter Sets for Combined Multiple Recursive Random Number Generators'', 
 * Operations Research, 47, 1 (1999), 159--164.

 * P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton, 
 * ``An Objected-Oriented Random-Number Package with Many Long Streams and Substreams'', 
 * Operations Research, 50, 6 (2002), 1073--1075
 *
\***********************************************************************/

 
#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string>
#include <cstdint>

class RngStream
{
  public:
    RngStream (const char *name = "");

    static bool SetPackageSeed (const uint32_t seed[6]);

    void ResetStartStream ();

    void ResetStartSubstream ();

    void ResetNextSubstream ();

    void SetAntithetic (bool a);

    void IncreasedPrecis (bool incp);

    bool SetSeed (const uint32_t seed[6]);

    void AdvanceState (int32_t e, int32_t c);

    void GetState (uint32_t seed[6]) const;

    void WriteState () const;

    void WriteStateFull () const;

    double RandU01 ();

    int RandInt (int i, int j);

  private:
    double Cg[6], Bg[6], Ig[6];

    bool anti, incPrec;

    std::string name;

    static double nextSeed[6];

    double U01 ();

    double U01d ();
};
 
#endif
 

