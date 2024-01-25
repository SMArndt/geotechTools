"""
plot3Dgeo.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
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
# imports
# ---------------------------------------------------------------------------

import numpy as np
from stl import mesh
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# ---------------------------------------------------------------------------
# class plot3Dgeo()
# ---------------------------------------------------------------------------

class plot3Dgeo:

    def __init__(self, points: np.ndarray=None, stlMesh=None, var=3):
        """
        constructor for plot3D()
        
        points: np.ndarray of shape (N, 3 or more), with color in points[:,var]
        """

        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        #usescientificcolourmaps suggested by https://www.linkedin.com/in/lindsey-smith-17665622a/
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis'][2]

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        plt_colours = points[:,var]
    
        fax = ax.scatter(points[:,0],points[:,1],points[:,2], c=plt_colours , cmap=cmaps)

        if stlMesh:
            your_mesh = stlMesh
            ax.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.view_init(elev=12.5, azim=-22.5) # updates ax.azim, ax.elev on close

        plt.colorbar(fax)
        plt.show()

    # ~def __init__(self, points: np.ndarray, stlMesh, var)
