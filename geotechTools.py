"""
geotechTools.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------
from stlGeom import stlGeom
from plot3D import plot3D
from xyzData import xyzData

# ---------------------------------------------------------------------------
# example
# ---------------------------------------------------------------------------

exampleData = r'SampleDataset.csv'
exampleMine = r'SampleMineGeo.stl'

# ---------------------------------------------------------------------------
# read CSV file
# ---------------------------------------------------------------------------

f = xyzData(exampleData)

f.filterIPR(1.0)
f.filterNaN('local magnitude')

s = stlGeom(exampleMine)

plot3D(f.current, s.stlMesh, 10)