"""
geotechTools.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
Visit https://cavroc.com/ for more information on IUCM and StopeX

This file is part of geotechTools (https://github.com/SMArndt/geotechTools).

geotechTools is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation.

geotechTools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with geotechTools.
If not, see <https://www.gnu.org/licenses/>.
"""

# ---------------------------------------------------------------------------
# notes
# ---------------------------------------------------------------------------
# - this file is a Python script, not a module, showcasing example use of the libraries.
# - to install Python and packages like numpy, https://www.anaconda.com/ is recommended.
# - download all files in one folder, run geotechTools.py (Python prompt or iPython)

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------

from xyzData import *
from gridData import *
from plot3D import *
from stressUtils import *

# ---------------------------------------------------------------------------
# example: convert principal stress to cartesian and back, plot as cube in 3D
# ---------------------------------------------------------------------------

# StopeX interface: input complex stress as gradient in MPa/m, positive values
# Using SI units in Abaqus [Pa, kg, m] compressive stress is negative
# Principals is a list of 3 principal stresses [[S_n, plunge_n, trend_n]]

Principals = [[-0.041, 5., 351.], [-0.032, 50., 255.], [-0.022, 40., 85.]]

S=getCartesianStress(Principals)

E=getPrincipalStress(S)

for orientation in range(3):
    plot3DCube(Principals[orientation][1],Principals[orientation][2])

P=getStressOrientation(E[1]) # E[1] is the eigenvector matrix

for dip,trend in P:
    plot3DCube(dip,trend)

# ---------------------------------------------------------------------------
# Update on the extractStress() example from 'Don't Stress' using new methods
# ---------------------------------------------------------------------------

# read x,y,z data with stress columns
x = xyzData(r'regular_stress.csv')

# extract stress components using the new method extractStress()
stress=x.extractStress(indices='xyz')

# perform eigenvalue analysis
eigens=[]
for s_vec in stress:
    T=unpackStress(s_vec)                # new method unpackStress()
    e_val, e_vec = getPrincipalStress(T) # new method getPrincipalStress()
    eigens.append(e_val[0])              # S1

orientations=getStressOrientation(e_vec) # new method getStressOrientation()

# quick way to plot S1 (in last column, t.shape[1]-1)
t=np.hstack((x.current,np.array(eigens).reshape(-1,1)))
plot3D(t,t.shape[1]-1)