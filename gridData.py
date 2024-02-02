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

from re import S
import numpy as np
from scipy.spatial import KDTree

from xyzData import *
import config

# ---------------------------------------------------------------------------
# class gridData()
# ---------------------------------------------------------------------------

class gridData:

	def __init__(self, data, cellSize=None, sparse=True, bBox=None):
		"""
		constructor for gridData()
		"""

		self.cells = {}             # dict {(i,j,k): data}
		self.shape = (0,0,0)        # (nx,ny,nz)
		self.sparse = sparse        # True (default) if not every i,j,k defined
		self.cellSize = cellSize    # edge length of grid cell - (None):auto
		self.order = 2.0            # target points per cell for auto == order^3
									# examples 2:8, 2.75:20.8, 3.7:50
		self.data = data            # np.ndarray, 

		# source data - columns 0,1,2 must be x,y,z
		# -----------------------------------------
		if isinstance(data,np.ndarray):
			self.data = data        
			self.index = False
		elif isinstance(data,xyzData):
			self.data = data.current
			self.index = data.index
		else:
			print ('Error: requires np.ndarray or xyzData')
			raise TypeError
		if (self.data.ndim!=2) or (self.data.shape[1]<3):
			print ('Error: requires ndim==2 and shape(N, 3 or more')
			raise TypeError

		# determine bounding box
		# ----------------------
		if bBox==None:
			self.bBox = array3D_BBox(self.data)
		else:
			self.bBox = bBox

		# auto cell size algorithm
		# ------------------------
		if self.cellSize==None:

			rLen = pow(pow((self.bBox[1][0]-self.bBox[0][0]),2) + \
					   pow((self.bBox[1][1]-self.bBox[0][1]),2) + \
					   pow((self.bBox[1][2]-self.bBox[0][2]),2), 0.5) / pow(3,.5)

			nXYZ = pow(len(self.data),1./3.) # avg data points per axis in rLen, total is nXYZ^3
			self.cellSize = round(self.order/(nXYZ/rLen),1)

		# grid dimensions
		# ---------------
		(nx,ny,nz) = (int( (self.bBox[1][0]-self.bBox[0][0]) / self.cellSize ) + 1, \
					  int( (self.bBox[1][1]-self.bBox[0][1]) / self.cellSize ) + 1, \
					  int( (self.bBox[1][2]-self.bBox[0][2]) / self.cellSize ) + 1)
		self.shape = (nx,ny,nz)

		# create empty cells if not sparse
		# --------------------------------
		if not(self.sparse):
			self.cells = {(i,j,k): [] for i in range(nx) for j in range(ny) for k in range(nz)}

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

		if config.verbose:
			print(self)

			avg_hist={} # histogramm
			for ijk in self.cells.keys():
				if len(self.cells[ijk]) in avg_hist:
					avg_hist[len(self.cells[ijk])]+=1
				else:
					avg_hist[len(self.cells[ijk])]=1
			for i in range(1,13):
				if i in avg_hist:
					print ("%2d" % i, "%4d" % avg_hist[i], '#' * int(80/avg_hist[1]*avg_hist[i]))

	# ~def __init__(self, data: np.ndarray, cellSize=None, sparse=True)

	def cellCount(self, minN=0):
		"""
		method to reverse ijk to np.array() of shape (N, 4)
		with [,0],[,1],[,2] = x,y,z and [,3] = number of points in cell
		"""
		
		self.gridXYZ = []
		
		for (i,j,k) in self.cells.keys():

			x = (i+0.5)*self.cellSize + self.bBox[0][0]
			y = (j+0.5)*self.cellSize + self.bBox[0][1]
			z = (k+0.5)*self.cellSize + self.bBox[0][2]
				
			self.gridXYZ.append(np.array([x,y,z,len(self.cells[(i,j,k)])]))

		self.gridXYZ = np.vstack(self.gridXYZ)

		return self.gridXYZ[self.gridXYZ[:,3] > minN]

	# ~def cellCount()

	def fillGrid(self):
		"""
		method to fill sparse grid with empty cells
		""" 

		(nx,ny,nz) = self.shape

		for i in range(nx):
			for j in range(ny):
				for k in range(nz):
					if (i,j,k) not in self.cells:
						self.cells[(i,j,k)]=[]  

		self.sparse = False        
		self.cellCount() # update self.gridXYZ
		
	# ~def fillGrid()
		
	def __str__(self):
		(nx,ny,nz) = self.shape
		return f"gridData object of {self.data.shape} with cellsize {self.cellSize}: " + \
			   f"{len(self.cells.keys())} active cells in {nx,ny,nz} = {nx*ny*nz} grid "