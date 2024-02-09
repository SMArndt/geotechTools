"""
plot3D.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
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

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from stressUtils import rot_x, rot_y, rot_z
from gridData import *
import config

# ---------------------------------------------------------------------------
# class plot3DVoxel()
# ---------------------------------------------------------------------------

class plot3DVoxel:

    def __init__(self, g: gridData, elev=12.5, azim=-22.5):
        """
        constructor for plot3DVoxel()
        """

        if not isinstance(g,gridData):
            print ('Error: plot3DVoxel needs type gridData')
            raise TypeError

        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo'][-1]
        cmap = matplotlib.colormaps[cmaps]

        i,j,k = np.indices(g.shape) # for all indices in grid, not g.cells.keys()
        n_vox = np.zeros(g.shape, dtype=bool)
        for ijk in g.cells.keys(): n_vox[ijk]=True

        cmax = 0 # max number of data points in a cell for colormap
        for ijk in g.cells.keys():
            cmax = max(cmax,len(g.cells[ijk]))

        # initialise arrays for face, edge colours 
        n_map = np.zeros(g.shape+(1,), dtype=object)
        n_edge= np.zeros(g.shape+(1,), dtype=object)

        # face, edge colours and alpha (transparencey) for each voxel
        for ijk in g.cells.keys():
            cval = min(0.0+len(g.cells[ijk])/cmax*2,1.0)
            calp = min(0.0+len(g.cells[ijk])/cmax*2,0.75)
            chex = hex(256+int(calp*255))[-2:]
            n_map[ijk] = matplotlib.colors.rgb2hex(cmap(cval))+chex
            n_edge[ijk]= matplotlib.colors.rgb2hex(cmap(cval))+chex

        # shrinking the voxels 
        # https://matplotlib.org/stable/gallery/mplot3d/voxels_numpy_logo.html

        def explode(data):
            size = np.array(data.shape)*2
            data_e = np.zeros(size - 1, dtype=data.dtype)
            data_e[::2, ::2, ::2] = data
            return data_e

        # upscale the above voxel image, leaving gaps
        n_vox2 = explode(n_vox)
        n_map2 = explode(n_map)
        n_edge2 = explode(n_edge)

        # shrink the gaps
        x, y, z = np.indices(np.array(n_vox2.shape) + 1).astype(float) // 2
        x[0::2, :, :] = (x[0::2, :, :]+0.1)*g.cellSize + g.bBox[0][0]
        y[:, 0::2, :] = (y[:, 0::2, :]+0.1)*g.cellSize + g.bBox[0][1]
        z[:, :, 0::2] = (z[:, :, 0::2]+0.1)*g.cellSize + g.bBox[0][2]
        x[1::2, :, :] = (x[1::2, :, :]+0.9)*g.cellSize + g.bBox[0][0]
        y[:, 1::2, :] = (y[:, 1::2, :]+0.9)*g.cellSize + g.bBox[0][1]
        z[:, :, 1::2] = (z[:, :, 1::2]+0.9)*g.cellSize + g.bBox[0][2]

        # create and show plot
        ax = plt.figure().add_subplot(projection='3d')
        ax.voxels(x, y, z, n_vox2, facecolors=n_map2, edgecolors=n_edge2) # ax.voxels(n_vox, facecolors=n_map)

        ax.view_init(elev=12.5, azim=-22.5)
        ax.set_aspect('equal')

        plt.show()

    # ~def __init__()

    def __str__(self):
        return f"plot3DVoxel"

# ---------------------------------------------------------------------------
# class plot3D()
# ---------------------------------------------------------------------------

class plot3D:

    def __init__(self, points: np.ndarray, var, vmin=None, vmax=None, elev=12.5, azim=-22.5):
        """
        constructor for plot3D()
        """

        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo'][-1]

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        plt_colours = points[:,var]
    
        fax = ax.scatter(points[:,0],points[:,1],points[:,2], vmin=vmin, vmax=vmax, c=plt_colours , cmap=cmaps)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        limits = array3D_BBox(points)
        
        ax.axes.set_xlim3d(left=limits[0][0], right=limits[1][0])
        ax.axes.set_ylim3d(bottom=limits[0][1], top=limits[1][1])
        ax.axes.set_zlim3d(bottom=limits[0][2], top=limits[1][2])

        ax.view_init(elev=elev, azim=azim, roll=0) # updates ax.azim, ax.elev on close
        ax.set_aspect('equal')
        plt.colorbar(fax)
        plt.show()

        if config.verbose:
            print(f"elev {ax.elev}, azim {ax.azim}")
            
    # ~def __init__()

    def __str__(self):
        return f"plot3D"

# ---------------------------------------------------------------------------
# class plot3DCube()
# ---------------------------------------------------------------------------

class plot3DCube:
    
    def __init__(self, dip, trend):
        """
        constructor for plot3DCube() 
        
        plots a rotated cube (dip and trend) and a reference cube
        """
        vertices=True

        # define a cube (nodes and faces)
        cNodes=[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],[0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]]
        cFaces=[(0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)]

        cNodes=np.array(cNodes)
        pNodes=np.array(cNodes)

        # rotation of the cube around center (0.5,0.5,0.5)
        R1 = rot_z(np.radians(-trend))
        R2 = rot_x(np.radians(-dip))
        R = np.dot(R1,R2)
        for i in range(len(cNodes)):
            pNodes[i]=np.dot(R,cNodes[i]-[0.5,0.5,0.5])+[0.5,0.5,0.5]

        # create plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # create list of edges and faces
        fList=[]
        for f in cFaces:
            fList.append([cNodes[f[0]],cNodes[f[1]],cNodes[f[2]],cNodes[f[3]]])
        
        eList=[]
        for f in cFaces:
            for i in range(4):
                if (f[i]<f[(i+1)%4]):
                    eList.append([cNodes[f[i]],cNodes[f[(i+1)%4]]])

        # plot edges
        for e in eList:
            ax.plot([e[0][0],e[1][0]],[e[0][1],e[1][1]],[e[0][2],e[1][2]],'k-')
 
        # plot faces
        faces = Poly3DCollection(fList, alpha=0.05)
        faces.set_facecolor('k')
        ax.add_collection3d(faces)

        # plot rotated cube
        eList=[]
        for f in cFaces:
            for i in range(4):
                if (f[i]<f[(i+1)%4]):
                    eList.append([pNodes[f[i]],pNodes[f[(i+1)%4]]])
        # plot edges
        for e in eList:
            ax.plot([e[0][0],e[1][0]],[e[0][1],e[1][1]],[e[0][2],e[1][2]],'k-')

        # color faces for each orientation (x,y,z)
        orientations = [['r',5],['g',2],['b',0]]
        for o in orientations:
            f=cFaces[o[1]]
            fList=[[pNodes[f[0]],pNodes[f[1]],pNodes[f[2]],pNodes[f[3]]]]
            faces = Poly3DCollection(fList, alpha=0.75)
            faces.set_facecolor(o[0])
            ax.add_collection3d(faces)

        ax.set_aspect('equal')
        plt.show()

    # ~def __init__()

# ---------------------------------------------------------------------------
# class plot3Dgeo() - requires stl, a non-standard package: pip install numpy-stl
# ---------------------------------------------------------------------------

try:
    from stl import mesh
    from mpl_toolkits import mplot3d

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

except ImportError:
    print('stl import failed. plot3Dgeo() requires numpy-stl')