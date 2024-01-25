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
from xyzData import xyzData
from plot3D import plot3D

# ---------------------------------------------------------------------------
# example: read CSV file
# ---------------------------------------------------------------------------

f = xyzData(r'SampleDataset.csv')

# ---------------------------------------------------------------------------
# mapData() example: mapping of 'colour-xz' from random points to exampleData
# ---------------------------------------------------------------------------

g = xyzData(r'random_colour.csv')
plot3D(g.current,g.index['colour-yz'])

f.mapData(g)
plot3D(f.current,f.index['colour-yz'])

# ---------------------------------------------------------------------------
# optional: STL geometry, requires manual install of numpy-stl package
# ---------------------------------------------------------------------------

# exampleMine = r'SampleMineGeo.stl'

# from stlGeom import stlGeom
# from plot3Dgeo import plot3Dgeo

# s = stlGeom(exampleMine)

# f.filterIPR(1.0)
# f.filterNaN('local magnitude')

# plot3Dgeo(f.current, s.stlMesh, f.index['local magnitude'])