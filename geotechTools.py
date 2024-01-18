"""
geotechTools.py

- to install Python, https://www.anaconda.com/ is recommended
- download all files to same folder and run geotechTools.py
"""

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