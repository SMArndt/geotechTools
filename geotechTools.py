"""
geotechTools.py

Example usage of xyzData and stlGeom classes, 
download all files to same folder and run geotechTools.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------
from xyzData import xyzData
from plot3D import plot3D

# ---------------------------------------------------------------------------
# example: read CSV file
# ---------------------------------------------------------------------------

from xyzData import xyzData

exampleData = r'SampleDataset.csv'

f = xyzData(exampleData)

f.filterIPR(1.0)
f.filterNaN('local magnitude')

plot3D(f.current, 10)

# ---------------------------------------------------------------------------
# optional: STL geometry, requires manual install of numpy-stl package
# ---------------------------------------------------------------------------

# exampleMine = r'SampleMineGeo.stl'

# from stlGeom import stlGeom
# from plot3Dgeo import plot3Dgeo

# s = stlGeom(exampleMine)
# plot3Dgeo(f.current, s.stlMesh, 10)