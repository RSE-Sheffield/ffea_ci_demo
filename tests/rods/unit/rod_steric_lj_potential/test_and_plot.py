#!/usr/bin/env python2
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

# This script runs the FFEA test, plots some stuff, and passes judgement (0|1)
# based     on analysis of the data

# Stop matplotlib complaining about $DISPLAY not being set
#from site import enablerlcompleter
import os
from pickletools import uint1
import sys
import matplotlib.pyplot as plt
import subprocess
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

# Allows saving figure without X-server
plt.switch_backend('agg')

def plot(surf_surf_distance, energy_ffea, r_min, kT, energy_anal, radius, eps):

    plt.plot(surf_surf_distance, energy_ffea, 'k-', lw=2.5, label="FFEA potential", zorder=1)
    plt.plot(surf_surf_distance, energy_anal, '-', color="#648FFF", lw=1.3, label="Analytical potential", zorder=2)

    plt.plot(surf_surf_distance, np.zeros(surf_surf_distance.size), 'k:', lw=0.8, zorder=0)
    plt.plot(np.zeros(energy_ffea.size), energy_ffea, 'k:', lw=0.8, zorder=0)

    plt.plot(r_min*np.ones(energy_ffea.size), energy_ffea, '--', color="#DC267F",  lw=1.5, label="$r = r_{min}$", zorder=0)

    plt.title("$r_{min} = $" + f"{r_min / radius:.2f}*R")
    plt.xlabel("Surface-surface distance ($r_{min}$)")
    plt.ylabel("Energy ($kT$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("Steric_VDW_Modified_Potential.png", dpi=300)

def lj_6_12(r, sigma, eps):
	"""
	Typical Lennard-Jones potential

	r:     surface-surface distance
	sigma: distance where U = 0
	r_min:  distance where U = -eps
	eps:   potential well depth
	"""
	r_6 = np.power(sigma/r, 6)
	r_12 = np.power(sigma/r, 12)
	return 4 * eps * (r_12 - r_6)

def lj_interp(r, r_min, eps):
	"""Interpolative region, ben Hanson thesis p74-77"""
	r_2 = 3 * np.power(r/r_min, 2)
	r_3 = 2 * np.power(r/r_min, 3)
	return eps * (r_3 - r_2)

def steric(r, k):
	"""Quadratic repulsive potential"""
	return k * r * r

def auto():
    test_name = os.getcwd().split("/")[-1]
    # The FFEA test is called from within python, to allow plotting of results
    subprocess.call([sys.argv[1], test_name + ".ffeatest"])

    const = pd.read_csv("const.csv").to_dict("list")
    radius = const["R"][0]
    r_min = const["r_min"][0]
    eps = const["eps"][0]
    k = const["k"][0]
    kT = const["kT"][0]
    sigma = const["sigma"][0]

    data = pd.read_csv("log.csv")
    r = data["r"].values
    u_ffea = data["u"].values
    u_anal = np.zeros(r.size)

    for i, ri in enumerate(r):
        if ri < 0:
            u_anal[i] = steric(ri, k)
        elif ri < r_min and ri >= 0:
            u_anal[i] = lj_interp(ri, r_min, eps)
        elif ri > r_min:
            u_anal[i] = lj_6_12(ri, sigma, eps)

    err = np.max(np.abs(u_anal - u_ffea)) / kT
    print(f"\nLargest error in FFEA vs. analytical rod-rod energies = {err:e}")

    plot(r / r_min, u_ffea / kT, r_min / r_min, kT, u_anal / kT, radius / r_min, eps / kT)

    if err < 1e-3:
        return 0
    else:
        return 1


if __name__ == "__main__":
    error_status = auto()
    sys.exit(error_status)
