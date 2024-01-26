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

# ---------------------------------------------------------------------------
# example: read CSV file
# ---------------------------------------------------------------------------

f = xyzData(r'SampleDataset.csv')
f.filterIPR(1.0)
f.filterNaN('local magnitude')
plot3D(f.current,f.index['local magnitude'])

# ---------------------------------------------------------------------------
# create grid with cell size 50m
# ---------------------------------------------------------------------------

g = gridData(f.current, cellSize=50.)
plot3DVoxel(g)