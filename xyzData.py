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

import time
from datetime import datetime
import csv
import config

import numpy as np
import pandas as pd

from scipy.spatial import KDTree
from stressUtils import *

# ---------------------------------------------------------------------------
# functions
# ---------------------------------------------------------------------------

def array3D_BBox(a):
    """
    bounding box of np.array() of shape (N, 3 or more) with a[,0],[,1],[,2] = x,y,z
    
    arguments:
    -a: np.array()
    
    returns: tuple ((x0,y0,z0),(x1,y1,z1))
    """
    
    if a.shape[0] > 0:
        return ((min(a[:,0]), min(a[:,1]), min(a[:,2])), (max(a[:,0]), max(a[:,1]), max(a[:,2])))
    else:
        return False

# ~def array3D_BBox(a)

def array3D_IPR(a, p_IPR=25.0):
    """
    filter in interpercentile range of np.array() of shape (N, 3 or more)
    
    arguments:
    -a: np.array()
    -p_IPR: float, interpercentile range (0,50)
    
    returns: np.array() of a.shape
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

# ~def array3D_IPR(a, p_IPR=25.0)

# ---------------------------------------------------------------------------
# class xyzData()
# ---------------------------------------------------------------------------

class xyzData:

    def __init__(self, fileName=None):
        """
        constructor for xyzData()
        class for np.array(dtype=float) with required columns x,y,z in [0,1,2]
        all string columns from csv file need interpretation methods, i.e. date
        
        arguments:
        -fileName:      (optional) argument fileName will invoke read() method
        """

        self.fileName = fileName    # self explanatory
        self.pData = np.array([])   # raw data set, np.array()
        self.current = np.array([]) # current data set (filtered), np.array()
        
                                    # current data bounding box ... nerd alert
        self.bBox = ((-np.inf,-np.inf,-np.inf),(np.inf,np.inf,np.inf))
        self.maxCol = 0             # max column index saved in pData - shape is (N,maxCol+1)
        self.dateCol = False
        self.stepCol = False
        
        if not(self.fileName==None):
            self.read(self.fileName)

    # ~def __init__(self, fileName=None)

    def __str__(self):
        return f"{self.fileName}, {len(self.pData)} Lines, {len(self.current)} current Points"

    def read(self, fileName, xyzCol=['x','y','z'], dateCol='Date', dateFormat='iso-8601', \
             limit=None, rmCRLF=True):
        """
        method to read xyzData
        
        arguments:
        -fileName:      string, csv file name
        -xyzCol:        list of strings, column names for x,y,z
        -dateCol:       string, column name for date
        -limit:         integer, limit number of lines
        -rmCRLF:        boolean, remove CRLF from csv file
        """

        self.fileName = fileName
        pDataArray = []   # raw data set, np.array()

        # argument xyzCol can be string or list (catch copy paste)
        if isinstance(xyzCol,str): xyzCol=xyzCol.split(',')
        elif isinstance(xyzCol,list) and len(xyzCol)==1: xyzCol=xyzCol[0].split(',')
        elif isinstance(xyzCol,list) and len(xyzCol)==3:pass
        else:
            print (f"xyzCol {xyzCol} not valid")
            raise ValueError

        # index for listed entries in pData with mandatory x,y,z in [0,1,2]
        self.index = {xyzCol[0]:0,xyzCol[1]:1,xyzCol[2]:2}
        if config.verbose:
            print ("xyz columns are: ", "x:", xyzCol[0], "y:", xyzCol[1], "z:", xyzCol[2])
        
        # optional excluded columns from csv file
        self.exclude = \
            ['id','date','location residual','apparent stress', 'static stress drop' ,'dynamic stress drop' , \
             's wave frequency' ,'p wave energy' ,'s wave energy' ,'s:p energy ratio' , \
             'total radiated energy' ,'p(outlier)'] \
             + list(self.index.keys()) 

        self.dateCol = dateCol      # date column name
        self.xyzCol = xyzCol        # x,y,z column names
        self.csvCol = {}            # origin column in csv file
        self.maxCol = 0             # max column index saved in pData - shape is (N,maxCol+1)

        t0 = time.time()
        with open(fileName, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"') # defaults
    
            i,k = 0,0 # i: line counter, k: invalid lines
            
            for row in csv_reader:
                if i==0: # headers
                    if config.verbose: print ("reading column headers ...")

                    # listed headers, save currentCol in dictionary & update maxCol
                    currentCol=0
                    for csv_col_head in row:
                        if rmCRLF:
                            # remove leading spaces, replace CRLF with space, remove double spaces
                            csv_col = csv_col_head.lstrip().replace('\r\n', ' ').replace('  ', ' ')
                        else:
                            # remove leading spaces
                            csv_col = csv_col_head.lstrip()

                        # lower case headers for 'x','y','z'
                        if csv_col in ['X','Y','Z']: csv_col=csv_col.lower()

                        if csv_col in self.index.keys():
                            self.csvCol[csv_col]=currentCol
                            self.maxCol=max(self.maxCol,self.index[csv_col])

                            if config.verbose:
                                print (csv_col, "listed in column", self.csvCol[csv_col], \
                                       "stored in", self.index[csv_col] )
                        currentCol+=1

                    # unlisted headers, locate available slots
                    listed=[]
                    for j in self.index: listed.append(self.index[j])
                    nextIndex = list(range(len(row)))
                    for j in listed: nextIndex.remove(j)
                    nextIndex.sort()

                    # unlisted headers, save currentCol in dictionary & update maxCol
                    currentCol=0                    
                    for csv_col_head in row:
                        if rmCRLF:
                            # remove leading spaces, replace CRLF with space, remove double spaces
                            csv_col = csv_col_head.lstrip().replace('\r\n', ' ').replace('  ', ' ')
                        else:
                            # remove leading spaces
                            csv_col = csv_col_head.lstrip()
                        
                        # lower case headers for 'x','y','z'
                        if csv_col in ['X','Y','Z']: csv_col=csv_col.lower()

                        if csv_col not in self.exclude:
                            self.csvCol[csv_col]=currentCol
                            self.index[csv_col]=nextIndex[0] # next available slot
                            self.maxCol=max(self.maxCol,self.index[csv_col])
                            nextIndex.remove(nextIndex[0])

                            if config.verbose:
                                print (csv_col, "unlisted in column", self.csvCol[csv_col], \
                                       "stored in", self.index[csv_col] )
                        currentCol+=1
                        
                    # avoid repeat lookup
                    xpos = self.csvCol[xyzCol[0]]
                    ypos = self.csvCol[xyzCol[1]]
                    zpos = self.csvCol[xyzCol[2]]
                    
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
                        if csv_col not in (xyzCol+['id',dateCol]):
                            colIndex = self.csvCol[csv_col]
                            try:
                                rowData[self.index[csv_col]]=float(row[colIndex])
                            except:
                                pass
                        elif csv_col == dateCol:
                            colIndex = self.csvCol[csv_col]
                            try:
                                if dateFormat == 'iso-8601':
                                    # use np.datetime64 directly
                                    rowData[self.index[csv_col]]=np.datetime64(row[colIndex].replace("/","-"))
                                else:
                                    # use datetime.strptime and dateFormat, replace / with -
                                    rowTime = datetime.strptime(row[colIndex].replace("/","-"), dateFormat)
                                    rowData[self.index[csv_col]] = np.datetime64(rowTime)
                            except:
                                # if dateCol is defined and not valid, line is invalid
                                valid = False
                                pass

                    if valid:
                        pDataArray.append(np.array(rowData))

                i+=1
                if not(limit==None):
                    if i>=(limit):
                        break

        if config.verbose:
            print (f"time: {time.time()-t0} seconds")
            print (f"{i} Lines, {self.maxCol+1} columns")
            print (f"invalid data in {k} lines")
            
        self.pData = np.vstack(pDataArray)

        self.current = self.pData
        self.bBox = array3D_BBox(self.current)
        
    # ~read(self, fileName)

    def reset(self):
        """
        method to reset self.current to self.pData
        """

        if len(self.pData)>0:
            self.current = self.pData
            self.bBox = array3D_BBox(self.current)
            self.maxCol = self.pData.shape[1]-1
        else:
            print (f'no data')
            raise ValueError
        
    # ~reset(self)

    def coordShift(self, dx=0.0, dy=0.0, dz=0.0):
        """
        method to shift x,y,z in self.current
        
        arguments:
        -dx, dy, dz: float, shift x,y,z
        """

        self.current[:,0]+=dx
        self.current[:,1]+=dy
        self.current[:,2]+=dz
        
        self.bBox = array3D_BBox(self.current)

    # ~coordShift(self, dx=0.0, dy=0.0, dz=0.0)

    def fromMetric(self, zShift=0.0, xScale=0.3048, yScale=0.3048, zScale=0.3048):
        """
        method to convert z to metric
        
        arguments:
        -xScale, yScale, zScale, zShift: float
        """
        
        # shift first, then scale
        self.current[:,2]=self.current[:,2]+zShift

        self.current[:,0]*=xScale
        self.current[:,1]*=yScale
        self.current[:,2]*=zScale

        self.bBox = array3D_BBox(self.current)

    def fromArray(self, data, index=False):
        """
        method to create xyzData from np.array()
        """

        self.xyzCol = ['x','y','z']
        self.index = {'x':0,'y':1,'z':2}

        if isinstance(data, np.ndarray):
            self.maxCol = data.shape[1]-1
            self.pData = data
            self.current = self.pData
            self.bBox = array3D_BBox(self.current)
            
        if index and isinstance(index, dict):
            self.index.update(index)
            self.maxCol = max(self.index.values())
        else:
            for i in range(3,self.pData.shape[1]):
                self.index[f'col{i}']=i+1
                
    # ~fromArray(self, data, index=False)

    def filterIPR(self, p_IPR):
        """
        method to filter outliers using interpercentile range p_IPR = (0,50)
        
        arguments:
        -p_IPR: float, interpercentile range (0,50)
        """
        
        l0 = len(self.current)
        self.current = array3D_IPR(self.pData, p_IPR)
        self.bBox = array3D_BBox(self.current)
        
        if config.verbose:
            print (f"filterIPR ({p_IPR}%-{100-p_IPR}%) removed {l0-len(self.current)} lines")
        
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
                print (f'{colStr} not in source index')
                raise TypeError
        
        l0 = len(self.current)
        self.current = self.current[np.invert(self.current[:,col]!=self.current[:,col])]
        
        if config.verbose: print (f"filterNaN [{colStr}] removed {l0-len(self.current)} lines")
        
        self.bBox = array3D_BBox(self.current)
            
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

        if config.verbose: print (f"filterBBox {bBox} offset {offset} removed {l0-len(self.current)} lines")
            
        self.bBox = array3D_BBox(self.current)
        return self.current
        
    # ~filterBBox(self):

    def filterCol(self, col='x', minVal=False, maxVal=False):
        """
        method to filter on a column value inside minVal,maxVal
        
        arguments:
        -col:       integer, index / string: key for self.index[]
        -minVal:    float, minimum value
        -maxVal:    float, maximum value
        """

        if isinstance(col,str):
            try:
                colStr=col
                col=self.index[colStr]
            except:
                print(f"filterCol: '{colStr}' not found in self.index")
                return
            
        l0 = len(self.current)
        
        if not(minVal is False):
            self.current = self.current[(self.current[:,col]>=minVal)]
        if not(maxVal is False):
            self.current = self.current[(self.current[:,col]<=maxVal)]
        
        colStr=col # only for output
        if config.verbose: print (f"filterCol '{colStr}',{minVal},{maxVal} removed {l0-len(self.current)} lines, " + \
                                  f"now {len(self.current)}")

        if len(self.current)==0:
            print (f"filterCol '{col}' removed all lines")
            raise ValueError

        self.bBox = array3D_BBox(self.current)
    
    # ~filterCol(self, col='x', minVal=False, maxVal=False)

    def extractArrayN4(self, col):
        """
        method to extract np.array() of shape (N, 4)
        
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

    def extractStress(self, indices='default', xyz=False):
        """
        method to extract np.array() of shape (N, 6) or (N, 9) if xyz=True
        
        arguments:
        -indices:       stress column notation
        -xyz boolean:   True if xyz stress components are in array[0:2]
        """
        
        if indices=='default' or indices=='xyz':
            col=['sxx','syy', 'szz', 'sxy','sxz', 'syz']
        if indices=='123':
            col=['s11','s22', 's33', 's12','s13', 's23']
        
        for i in range(6):
            try:
                col[i]=self.index[col[i]]
            except:
                print(f"extractStress: '{col[i]}' not found in self.index")
                return

        if xyz:
            return np.hstack((self.current[:,0:3], \
                   np.hstack([self.current[:,col[i]].reshape(-1, 1) for i in range(6)])))
        else:
            return np.hstack([self.current[:,col[i]].reshape(-1, 1) for i in range(6)])

    # ~extractStress()

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