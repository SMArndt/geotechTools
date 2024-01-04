"""
geotechTools.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------
from stlGeom import stlGeom

# ---------------------------------------------------------------------------
# example
# ---------------------------------------------------------------------------

s1 = stlGeom()
s1.readPath(r'C:\Projects\AMIRA\Objective 1\Mine Geometry STL\Test_Flat', combine=True)
s1.write(r'C:\Projects\AMIRA\Objective 1\Mine Geometry STL\New_Flat.stl')

s2 = stlGeom()
s2.readPath(r'C:\Projects\AMIRA\Objective 1\Mine Geometry STL\Test_Walk', combine=True)
s2.write(r'C:\Projects\AMIRA\Objective 1\Mine Geometry STL\New_Walk.stl')