
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

import numpy as np
import sys

def readAndAverage(iFile, ini, end, th):
  H = []
  A = {}
  Fields = th.keys()
  with open(iFile, 'r') as sta:
    while (sta.readline() != "Measurements:\n"):
      continue

    for i in sta.readline().split():
      H.append([i.strip()])
    for line in sta:
      cnt = 0
      if line.count("RESTART"): continue
      for i in line.split():
        H[cnt].append(float(i))
        cnt += 1

  if (ini == 0): ini = 1
  for i in H:
    if i[0] not in Fields:
      continue
    print(i[0], abs(np.mean(i[ini:end])/th[i[0]] - 1), np.std(i[ini:end]), len(i[ini:end]))
    A[i[0]] = [abs(np.mean(i[ini:end])/th[i[0]] - 1), np.std(i[ini:end])]
  return A

nodes = 63
KbT = 4.11e-21
E = {"KineticEnergy":KbT*(3*nodes)/2, "StrainEnergy":KbT*(3*nodes - 6)/2}
print("K_th: ", E["KineticEnergy"])
print("E_th: ", E["StrainEnergy"])
Tol = {"KineticEnergy": 0.04, "StrainEnergy":0.03}
ini = 40
end = -1

iFile = ["sphere_63_120_mass_measurement.out"]

err = 0
for f in iFile:
  print(f)
  A = readAndAverage(f, ini, end, E)
  for s in Tol.keys():
    print(A[s])
    if A[s][0] < Tol[s]: 
      print(s, ": correct ", A[s][0], " < ", Tol[s])
    else:
      print(s, ": failed ", A[s][0], " > ", Tol[s])
      err = 1

  print("\n\n")


sys.exit(err)
