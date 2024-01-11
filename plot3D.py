"""
plot3D.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------

import numpy as np
from stl import mesh
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# ---------------------------------------------------------------------------
# class plot3D()
# ---------------------------------------------------------------------------

class plot3D:

    def __init__(self, points: np.ndarray, stlMesh, var):
        """
        constructor for plot3D()
        """

        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        cmap_rainbow = 'rainbow'

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        plt_colours = points[:,var]
    
        fax = ax.scatter(points[:,0],points[:,1],points[:,2], c=plt_colours , vmax=0, cmap=cmap_rainbow)

        your_mesh = stlMesh
        ax.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.view_init(elev=12.5, azim=-22.5) # updates ax.azim, ax.elev on close

        plt.colorbar(fax)
        plt.show()

    # ~def __init__(self, points: np.ndarray, stlMesh, var)