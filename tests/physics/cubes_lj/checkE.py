#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

import sys


def readAndStore(iFile):
  H = []
  dH = {}
  with open(iFile, 'r') as sta:
    while (sta.readline() != "Measurements:\n"):
      continue

    blob = ""
    readblob = False
    for i in sta.readline().split():
      if i == "|":
        readblob = True
        continue
      if readblob:
        blob = i
        readblob = False
        continue
      H.append([blob+i])
    for line in sta:
      cnt = 0
      if line.count("RESTART"): continue
      for i in line.split():
        H[cnt].append(float(i))
        cnt += 1

  for i in range(len(H)):
    dH[H[i][0]] = i
    H[i].pop(0)

  return H, dH


# Doing the double surface integral
#      \int_S1 \int_S2 E = 1/r**12 - 2/r**6 dS1 dS2
# numerically for 2 cubes of side 1 separated 1 gives -0.704384545,
# the result obtained using CUBA (http://www.feynarts.de/cuba/)
# Thus, being the well depth and equilibrium distance values
#    in the .lj file: 2e537e17 J/m^4 and 1e-9m, the energy of the
#    first conformation is:
knownValue = -1.7870e-19
errTol = 1e-2
fail = 0

## ## ## STRAIGHT ## ## ##
iFile = ["veryFinefine1face1face.fm"]
H, dH = readAndStore(iFile[0])

errRel = abs((H[dH["SurfSurfEnergy"]][-1] - knownValue)/knownValue)
if ( errRel > errTol ):
  print( "SurfSurfEnergy energy should be %e, but was found to be %e" % (knownValue, (H[dH["SurfSurfEnergy"]][-1])))
  print( "The relative error: %e was larger than the tolerance: %e" % (errRel, errTol))
  fail = 1

sys.exit(fail)
