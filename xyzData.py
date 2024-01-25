"""
xyzData.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
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
from scipy.spatial import KDTree
import time
import csv

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
# class xyzData()
# ---------------------------------------------------------------------------

class xyzData:

    def __init__(self, fileName=None):
        """
        constructor for xyzData()
        """

        self.fileName = fileName

        self.pData=[]   # raw data set, np.array()
        self.current=[] # current data set (filtered), np.array()
        self.bBox=[]    # current data bounding box

        self.index = {'x':0,'y':1,'z':2} # index storing self.pData
        
        self.exclude = \
            ['id','date','location residual','apparent stress', 'static stress drop' ,'dynamic stress drop' , \
             's wave frequency' ,'p wave energy' ,'s wave energy' ,'s:p energy ratio' , \
             'total radiated energy' ,'p(outlier)'] \
             + list(self.index.keys()) 

        self.csvCol={}  # origin column in csv file
        self.maxCol=0   # max column index saved in pData - shape is (N,maxCol+1)
        
        if self.fileName==None:
            pass
        else:
            self.read(self.fileName)

    # ~def __init__(self, fileName=None)

    def __str__(self):
        return f"{self.fileName}, {len(self.pData)} Lines, {len(self.current)} current Points"

    def read(self, fileName):
        """
        method to read xyzData
        """

        self.fileName=fileName
        self.maxCol = 0
        
        t0 = time.time()
        with open(fileName, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"') # defaults
    
            i,k = 0,0 # i: line counter, k: invalid lines
            
            for row in csv_reader:
                if i==0: # headers
                    if verbose: print ("reading column headers ...")

                    # listed headers, save j in dictionary & update maxCol
                    j=0
                    for csv_col_head in row:
                        # current column is j
                        csv_col = csv_col_head.lower() 
                        if csv_col in self.index.keys():
                            self.csvCol[csv_col]=j
                            self.maxCol=max(self.maxCol,self.index[csv_col])
                        j+=1

                    # unlisted headers, locate next available slot
                    listed=[]
                    for j in self.index: listed.append(self.index[j])
                    nextIndex = list(range(len(row)))
                    for j in listed: nextIndex.remove(j)
                    nextIndex.sort()

                    j=0                    
                    for csv_col_head in row:
                        csv_col = csv_col_head.lower() 
                        if csv_col not in self.exclude:
                            self.csvCol[csv_col]=j
                            self.index[csv_col]=nextIndex[0] # next available slot
                            self.maxCol=max(self.maxCol,self.index[csv_col])
                            nextIndex.remove(nextIndex[0])
                        j+=1
                        
                    # avoid repeat lookup
                    xpos = self.csvCol['x']
                    ypos = self.csvCol['y']
                    zpos = self.csvCol['z']
                    
                elif i>0: # read data
                    valid = True

                    try:
                        # x,y,z are required for valid line
                        x=float(row[xpos])
                        y=float(row[ypos])
                        z=float(row[zpos])
                    except:
                        valid = False
                        k+=1
                    
                    # row template with np.nan for missing data
                    rowData=[x,y,z]+[np.nan]*(self.maxCol-2)
                    
                    for csv_col in self.csvCol:
                        if csv_col not in ['x','y','z','id','date']:
                            j = self.csvCol[csv_col] # column index
                            try:
                                rowData[self.index[csv_col]]=float(row[j])
                            except:
                                pass

                    if valid:
                        self.pData.append(np.array(rowData))
                i+=1

        if verbose:
            print (f"time: {time.time()-t0} seconds")
            print (f"{i} Lines, {self.maxCol+1} columns")
            print (f"invalid data in {k} lines")
            
        self.pData = np.vstack(self.pData)

        self.current = self.pData
            
    # ~read(self, fileName)

    def filterIPR(self, p_IPR):
        """
        method to filter outliers using interpercentile range p_IPR = (0,50)
        """
        l0 = len(self.current)
        self.current = array3D_IPR(self.pData, p_IPR)
        self.bBox = array3D_BBox(self.current)
        
        if verbose:
            print (f"filterIPR ({p_IPR}%-{100-p_IPR}%) removed {l0-len(self.current)} lines")
        
        return self.current
        
    # ~filterIPR(self, p_IPR)

    def filterNaN(self, col):
        """
        method to filter on a column containing NaN
        
        arguments:
        -col integer: index / string: key for self.index[]
        """

        if isinstance(col,str):
            try:
                colStr=col
                col=self.index[colStr]
            except:
                return
            
        l0 = len(self.current)
        self.current = self.current[np.isfinite(self.current[:,col])]
 
        if verbose: print (f"filterNaN '{colStr}' [{col}] removed {l0-len(self.current)} lines")

        self.bBox = array3D_BBox(self.current)
        return self.current
        
    # ~filterNaN(self):

    def filterBBox(self, bBox, offset=0.0):
        """
        method to filter on a bounding box
        """
        ((x0,y0,z0),(x1,y1,z1)) = (bBox[0][0]-offset,bBox[0][1]-offset,bBox[0][2]-offset), \
                                  (bBox[1][0]+offset,bBox[1][1]+offset,bBox[1][2]+offset) 
        
        l0 = len(self.current)
        self.current = self.current[( \
            (self.current[:,0] > x0) & (self.current[:,0] < x1) & \
            (self.current[:,1] > y0) & (self.current[:,1] < y1) & \
            (self.current[:,2] > z0) & (self.current[:,2] < z1) )]

        if verbose: print (f"filterBBox {bBox} offset {offset} removed {l0-len(self.current)} lines")
            
        self.bBox = array3D_BBox(self.current)
        return self.current
        
    # ~filterBBox(self):

    def extractArrayN4(self, col):
        """
        method to extract np.array of shape (N, 4)
        
        arguments:
        -col integer: index / string: key for self.index[]
        """

        if isinstance(col,str):
            try:
                colStr=col
                col=self.index[colStr]
            except:
                print(f"extractArrayN4: '{col}' not found in self.index")
                return
            
        return np.hstack((self.current[:,0:3],self.current[:,col].reshape(-1,1)))
    
    # ~extractArrayN4(self, col)

    def mapData(self, source, newIndex='mapData-1', overwrite=True, maxDist=False, fill=np.nan):
        """
        method to map data from source to self using kdTree
        - one column from np.array of shape (N, 4) into newIndex
        - all columns from xyzData class object with new indices from source.index
        """
        
        # source data
        # -----------
        if isinstance(source,np.ndarray):
            sourceData = source
        elif isinstance(source,xyzData):
            sourceData = source.current
        else:
            return
        
        # target data
        # -----------
        targetData = self.current
        
        # kdTree
        # ------
        kdtree=KDTree(sourceData[:,0:3])
        dist,points=kdtree.query(targetData[:,0:3],1) # for ,2: points[i] becomes list

        if isinstance(source,np.ndarray): # map one column only

            if newIndex in self.index.keys():
                targetCol = self.index[newIndex]
            else:
                targetCol = self.maxCol+1
                self.index[newIndex]=targetCol
                overwrite=True
            if overwrite:
                # create new column with zeros
                targetData = np.hstack([targetData,np.zeros([len(targetData),1])])
    
                for i in range(len(points)):
                    if (maxDist is False) or (dist[i]<=maxDist):
                        targetData[i,targetCol]=sourceData[points[i],3]
                    else:
                        targetData[i,targetCol]=fill

        elif isinstance(source,xyzData): # map all columns

            # create new columns with zeros     
            for col in source.index.keys():
                if col not in self.index.keys():
                    self.maxCol+=1
                    self.index[col]=self.maxCol
                    targetData = np.hstack([targetData,np.zeros([len(targetData),1])])
            # map data
            for i in range(len(points)):    
                for col in source.index.keys():
                    if col not in ['x','y','z']:
                        if (col not in self.index.keys()) or overwrite:
                            if (maxDist is False) or (dist[i]<=maxDist):
                                targetData[i,self.index[col]]=sourceData[points[i],source.index[col]]
                            else:
                                targetData[i,self.index[col]]=fill  
        else:
            return

        self.current = targetData

    # ~def mapData(self, source, newIndex='mapData-1')