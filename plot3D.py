"""
plot3D.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# ---------------------------------------------------------------------------
# class plot3D()
# ---------------------------------------------------------------------------

class plot3D:

    def __init__(self, points: np.ndarray, var, elev=12.5, azim=-22.5):
        """
        constructor for plot3D()
        
        points: np.ndarray of shape (N, 3 or more), with color in points[:,var]
        """

        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        #usescientificcolourmaps suggested by https://www.linkedin.com/in/lindsey-smith-17665622a/
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo'][-1]

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        plt_colours = points[:,var]
    
        fax = ax.scatter(points[:,0],points[:,1],points[:,2], c=plt_colours , cmap=cmaps)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.view_init(elev=elev, azim=azim, roll=0) # updates ax.azim, ax.elev on close

        plt.colorbar(fax)
        plt.show()

    # # ~def __init__()