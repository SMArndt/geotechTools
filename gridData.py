"""
gridData.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
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

# ---------------------------------------------------------------------------
# global variables
# ---------------------------------------------------------------------------

verbose = True

# ---------------------------------------------------------------------------
# functions
# ---------------------------------------------------------------------------

def array3D_BBox(a):
    """
    bounding box of np.array() of shape (N, 3 or more) with a[,0],[,1],[,2] = x,y,z
    
    returns tuple ((x0,y0,z0),(x1,y1,z1))
    """
    
    return ((min(a[:,0]), min(a[:,1]), min(a[:,2])), (max(a[:,0]), max(a[:,1]), max(a[:,2])))

def array3D_IPR(a, p_IPR=25.0):
    """
    filter in interpercentile range of np.array() of shape (N, 3 or more)
    """
    
    xyz_pmin, xyz_pmax = p_IPR, 100 - p_IPR # symmetric interpercentile range

    p_limits = ( \
        (np.percentile(a[:,0],xyz_pmin), np.percentile(a[:,0],xyz_pmax)), \
        (np.percentile(a[:,1],xyz_pmin), np.percentile(a[:,1],xyz_pmax)), \
        (np.percentile(a[:,2],xyz_pmin), np.percentile(a[:,2],xyz_pmax)))
        
    return a[( \
        (a[:,0] > p_limits[0][0]) & (a[:,0] < p_limits[0][1]) & \
        (a[:,1] > p_limits[1][0]) & (a[:,1] < p_limits[1][1]) & \
        (a[:,2] > p_limits[2][0]) & (a[:,2] < p_limits[2][1]) )]

# ---------------------------------------------------------------------------
# class gridData()
# ---------------------------------------------------------------------------

class gridData:

    def __init__(self, data: np.ndarray, cellSize=None, sparse=True):
        """
        constructor for gridData()
        """

        self.cells = {}             # dict {(i,j,k): data}
        self.shape = (0,0,0)        # (nx,ny,nz)
        self.sparse = sparse        # True if not every i,j,k defined
        self.cellSize = cellSize    # edge length of grid cell - (None):auto
        self.order = 2.0            # target points per cell for auto == order^3
                                    # examples 2:8, 2.75:20.8, 3.7:50
        self.data = data            # np.ndarray, columns 0,1,2 must be x,y,z

        if not(isinstance(data,np.ndarray) and (data.ndim==2) and (data.shape[1]>2)):
            print ('Error: data needs to be of type np.ndarray')
            raise TypeError

        # determine bounding box
        # ----------------------
        self.bBox = array3D_BBox(self.data)

        # auto cell size algorithm (IQR)
        # ------------------------------
        if self.cellSize==None:

            bBox = array3D_BBox(data)
            
            rLen = pow(pow((bBox[1][0]-bBox[0][0]),2) + \
                       pow((bBox[1][1]-bBox[0][1]),2) + \
                       pow((bBox[1][2]-bBox[0][2]),2), 0.5) / pow(3,.5)

            nXYZ = pow(len(data),1./3.) # avg data points per axis in rLen, total is nXYZ^3
            
            self.cellSize = round(self.order/(nXYZ/rLen),1)

        # grid dimensions
        # ---------------
        (nx,ny,nz) = (int( (self.bBox[1][0]-self.bBox[0][0]) / self.cellSize ) + 1, \
                      int( (self.bBox[1][1]-self.bBox[0][1]) / self.cellSize ) + 1, \
                      int( (self.bBox[1][2]-self.bBox[0][2]) / self.cellSize ) + 1)
        self.shape = (nx,ny,nz)

        # assign data
        # -----------
        for r in range(len(self.data)):
            x,y,z = self.data[r][0], self.data[r][1], self.data[r][2]
            i = int((x-self.bBox[0][0])/self.cellSize)
            j = int((y-self.bBox[0][1])/self.cellSize)
            k = int((z-self.bBox[0][2])/self.cellSize)
            if (i,j,k) in self.cells:
                self.cells[(i,j,k)].append(self.data[r])
            else:
                self.cells[(i,j,k)]=[self.data[r],]

        if verbose:
            avg_hist={} # histogramm
            for ijk in self.cells.keys():
                if len(self.cells[ijk]) in avg_hist:
                    avg_hist[len(self.cells[ijk])]+=1
                else:
                    avg_hist[len(self.cells[ijk])]=1
            print (f"Cell size {self.cellSize}: " + \
                   f"{len(self.cells.keys())} active cells in " + \
                   f"{nx,ny,nz} = {nx*ny*nz} grid ")
            for i in range(1,999):
                if i in avg_hist:
                    print ("%2d" % i, "%4d" % avg_hist[i], '#' * int(80/avg_hist[1]*avg_hist[i]))

    # ~def __init__(self, data: np.ndarray, cellSize=None, sparse=True)

    def cellCount(self, minN=0):
        """
        method for xyzGrid()

        - reverse ijk to np.array() of shape (N, 4)
          with [,0],[,1],[,2] = x,y,z and [,3] = number of points in cell
        """
        
        self.p_grid = []
        
        for (i,j,k) in self.cells.keys():

            x = (i+0.5)*self.cellSize + self.bBox[0][0]
            y = (j+0.5)*self.cellSize + self.bBox[0][1]
            z = (k+0.5)*self.cellSize + self.bBox[0][2]
                
            self.p_grid.append(np.array([x,y,z,len(self.cells[(i,j,k)])]))

        self.p_grid = np.vstack(self.p_grid)
        self.p_grid = self.p_grid[self.p_grid[:,3] > minN]

    # ~def cellCount()

    def __str__(self):
        return f"gridData object {self.data.shape} with cellsize {self.cellSize}"