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
# example: read stress from file, perform eigenvalue analysis, plot S1
# ---------------------------------------------------------------------------

# read x,y,z data with stress columns
x = xyzData(r'regular_stress.csv')
print(x.index)

# extract stress components using the new method extractStress()
stress=x.extractStress(indices='xyz')

# perform eigenvalue analysis
eigens=[]
for s in stress:
    T=unpackStress(s)                    # new method unpackStress()
    e_val, e_vec = getPrincipalStress(T) # new method getPrincipalStress()
    eigens.append(e_val[0]) # S1

orientations=getStressOrientation(e_vec) # new method getStressOrientation()

# quick way to plot S1 (in last column, t.shape[1]-1)
t=np.hstack((x.current,np.array(eigens).reshape(-1,1)))
plot3D(t,t.shape[1]-1)