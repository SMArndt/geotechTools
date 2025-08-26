"""
xyzData.py - Copyright 2025 S.M.Arndt, Cavroc Pty Ltd
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

from calendar import c
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

def array3D_BBox_ChatGPT(points):
    """
    Calculate the bounding box of a 3D array.

    Parameters:
    - points: NumPy array of points.

    Returns:
    - limits: Tuple of min and max limits for x, y, z axes.
    """
    min_vals = points.min(axis=0)
    max_vals = points.max(axis=0)
    return (min_vals, max_vals)

def bBoxInside(bBoxOuter, bBoxInner):
    """
    Returns True if bBoxInner is completely inside bBoxOuter.

    Parameters:
    - bBoxOuter: tuple of ((xmin, ymin, zmin), (xmax, ymax, zmax))
    - bBoxInner: tuple of ((xmin, ymin, zmin), (xmax, ymax, zmax))
    """
    for i in range(3):
        if not (bBoxOuter[0][i] <= bBoxInner[0][i] and bBoxInner[1][i] <= bBoxOuter[1][i]):
            return False
    return True

def roundBBox(bbox, digits=1):
    """
    Rounds all coordinates in a 3D bounding box to the specified number of decimal places.
    
    Parameters:
        bbox (tuple): A tuple of two 3D coordinate tuples, e.g., ((x1, y1, z1), (x2, y2, z2))
        digits (int): Number of decimal places to round to (default is 1)
    
    Returns:
        tuple: A new bbox with rounded float values.
    """
    return tuple(tuple(round(coord, digits) for coord in point) for point in bbox)

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

def indicesMap(indices='xyz'):

    mapStressStrain = {
        'xyz': ['sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz'],
        '123': ['s11', 's22', 's33', 's12', 's13', 's23'],
        'stress-Flac3D' : ['Sxx(MPa)','Syy(MPa)','Szz(MPa)','Sxy(MPa)','Sxz(MPa)','Syz(MPa)'],
        'strain-Flac3D' : ['Exx(-)','Eyy(-)','Ezz(-)','Exy(-)','Exz(-)','Eyz(-)'],
        'stress-xyz': ['sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz'],
        'stress-123': ['s11', 's22', 's33', 's12', 's13', 's23'],
        'STRESS-xyz': ['SXX', 'SYY', 'SZZ', 'SXY', 'SXZ', 'SYZ'],
        'STRESS-123': ['S11', 'S22', 'S33', 'S12', 'S13', 'S23'],
        'strain-xyz': ['exx', 'eyy', 'ezz', 'exy', 'exz', 'eyz'],
        'strain-123': ['e11', 'e22', 'e33', 'e12', 'e13', 'e23'],
        'STRAIN-xyz': ['EXX', 'EYY', 'EZZ', 'EXY', 'EXZ', 'EYZ'],
        'STRAIN-123': ['E11', 'E22', 'E33', 'E12', 'E13', 'E23'],
        'PBT' : ['P-Axis Plunge (\u00b0)', 'P-Axis Trend (\u00b0)', 'B-Axis Plunge (\u00b0)', 'B-Axis Trend (\u00b0)',  \
                 'T-Axis Plunge (\u00b0)', 'T-Axis Trend (\u00b0)']
    }

    if indices not in mapStressStrain.keys():
        print(f"'{indices}' no association found for indices")
        raise ValueError

    return mapStressStrain.get(indices, [])

# ~def indicesMap(indices='xyz')

def calcRvalue(row):
    """
    calculate R value for row

    arguments:
    -row: list of the three principals, i.e (P,B,T)

    returns: R value (float)
    """

    # list absolute values of the three principals, i.e (P,B,T)
    abs_values = [abs(x) for x in row]

    # maximum and minimum absolute value
    max_abs = max(abs_values)
    min_abs = min(abs_values)

    # compare the middle value with the maximum absolute value
    
    if abs_values[1] == max_abs:
        mid_abs = max(abs_values[0], abs_values[2])
    elif abs_values[1] == min_abs:
        mid_abs = min(abs_values[0], abs_values[2])
    else:
        mid_abs = abs_values[1]
    
    # divide difference between the max and mid by max and min values
    return (max_abs - mid_abs) / (max_abs - min_abs)

# ~def calcRvalue(row)

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
        
        Note: Common dateFormats are '%d-%m-%Y' or '%Y-%m-%d %H:%M'
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

                    # remove leading spaces, fix encoding issues, remove UTF-8 BOM
                    cleanRow = []
                    for csv_col in row:

                        # remove leading spaces
                        csv_col = csv_col.lstrip()
                        # fix encoding issues with degree symbol often found in mXrap files and
                        csv_col = csv_col.replace('Â°', '°')
                        # remove 'ï»¿' from UTF-8 BOM (Byte Order Mark) gone rogue
                        csv_col = csv_col.replace('ï»¿', '')
                        if rmCRLF:
                            # replace CRLF with space, remove double spaces
                            csv_col = csv_col.replace('\r\n', ' ').replace('  ', ' ')

                        cleanRow.append(csv_col)

                    row = cleanRow

                    # listed headers, save currentCol in dictionary & update maxCol
                    currentCol=0
                    for csv_col in row:

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
                    for csv_col in row:

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
                    
                    if valid:

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
                                    k+=1
                                    pass

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

    def fromImperial(self, zShift=0.0, xScale=0.3048, yScale=0.3048, zScale=0.3048, site=False):
        """
        Convert coordinates from imperial to metric. Z is shifted then scaled.
    
        Arguments:
        - zShift: float, shift applied to Z before scaling
        - xScale, yScale, zScale: float, scale factors (default 0.3048 for [ft]] to [m])
        - site: string, optional site name, site specific transformations
        """

        match site:
            case _:  # no transformations
                self.current[:, 2] = self.current[:, 2] + zShift

        # Apply scaling to all axes
        self.current[:, 0] *= xScale
        self.current[:, 1] *= yScale
        self.current[:, 2] *= zScale

        # Update bounding box
        self.bBox = array3D_BBox(self.current)

    def fromMetric(self, zShift=0.0, xScale=0.3048, yScale=0.3048, zScale=0.3048):
        """
        method to convert z to metric
        
        arguments:
        -xScale, yScale, zScale, zShift: float
        """
        
        # shift first, then scale

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
    
    def filterTimeRange(self, startDT='1970-01-01 00:00:00.000', endDT='2100-01-01 00:00:00.000'):
        """
        method to filter on a datetime range (not date range!)
        
        arguments:
        -startDT: string, format 'YYYY-MM-DD HH:MM:SS.sss'
        -endDT:   string, format 'YYYY-MM-DD HH:MM:SS.sss'
        """
        
        # note if startDT or endDT are in the format 'YYYY-MM-DD' datetime becomes a date
        # object, not a datetime, and comparison (<,>,=) fails
        
        if self.dateCol:
            dateIndex = self.index[self.dateCol]

            l0 = len(self.current)
            self.current = self.current[(self.current[:,dateIndex]>=np.datetime64(startDT)) & \
                                        (self.current[:,dateIndex]<=np.datetime64(endDT))]
        else:
            print (f"no date column")
            raise ValueError

        if config.verbose: print (f"filterDateRange {startDT} to {endDT} removed {l0-len(self.current)} lines")

        self.bBox = array3D_BBox(self.current)
    
    # ~filterTimeRange(self, startDT, endDT)

    def filterBBox(self, bBox, offset=0.0):
        """
        method to filter on a bounding box
        
        arguments:
        -bBox: tuple ((x0,y0,z0),(x1,y1,z1))
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
        -col:       integer, index / string: key for self.index[]
        
        returns: np.array()
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

    def extractStress(self, indices='xyz', xyz=False):
        """
        deprecated, use extractTensor() instead
        """
        print("--- extractStress: deprecated method, use extractTensor() instead!")

        return self.extractTensor(indices, xyz)

    # ~extractStress()

    def extractTensor(self, indices='xyz', xyz=False):
        """
        method to extract np.array() of shape (N, 6) or (N, 9) if xyz=True
        
        arguments:
        -indices:   string, stress column notation, 'xyz' or '123'
        -xyz:       boolean, True if xyz stress components are in array[0:2]
        
        returns: np.array()
        """

        if indices == 'PBT':
            print ("Warning: extractTensor 'PBT' is not a tensor, returns trend/plunge (6 components)")

        colNames = indicesMap(indices)
        col = [None] * 6

        for i in range(6):
            try:
                col[i]=self.index[colNames[i]]
            except:
                print(f"extractStress: '{col[i]}' not found in self.index")
                return

        if xyz:
            return np.hstack((self.current[:,0:3], \
                   np.hstack([self.current[:,col[i]].reshape(-1, 1) for i in range(6)])))
        else:
            return np.hstack([self.current[:,col[i]].reshape(-1, 1) for i in range(6)])

    # ~extractTensor()

    def extractPBT(self, colNames=False):
        """
        method to extract np.array() of shape (N, 9) for columns P,B,T
        
        arguments:
        -colNames
        
        returns: np.array()
        """

        # default column names for P,B,T
        if not colNames:
            colNames = ['P-Axis scale', 'P-Axis Trend (o)', 'P-Axis Plunge (o)', \
                        'B-Axis scale', 'B-Axis Trend (o)', 'B-Axis Plunge (o)', \
                        'T-Axis scale', 'T-Axis Trend (o)', 'T-Axis Plunge (o)']
        col = [None] * 9

        for i in range(9):
            try:
                print(colNames[i])
                col[i]=self.index[colNames[i]]
            except:
                print(f"extractPBT: '{col[i]}' not found in self.index")
                return

        return np.hstack([self.current[:,col[i]].reshape(-1, 1) for i in range(9)])

    # ~extractPBT()

    def addPrincipals(self, indices, prefix, dipDir=True, normals=True):
        """
        add principal values from tensor (six components) to self.current

        optional:
        - add dip, dipDir for each principal component to self.current
        - add normals for each principal component to self.current
        """

        # extractTensor returns np.array() of shape (N, 6)
        tArray = self.extractTensor(indices=indices, xyz=False)

        prinList, dipsList, vectList = [], [], [] # containers (lists) of length tArray

        for i in range(len(tArray)):
            # get full tensor
            T=unpackStress(tArray[i])

            # prinList[] - principal stresses from (eigenvalues, eigenvectors)
            sp, ev = getPrincipalStress(T)
            prinList.append(sp)

            # dipsList[] - plunge and trend for each principal stress
            dip_dipDir = getStressOrientation(ev)
            dipsList.append(dip_dipDir)

            if normals:
                # get normal vectors for each dip and direction
                vectList.append(getNormals(dip_dipDir))

        # ---------------------------------------------------------------------------
        # update self.current: principals, optional orientation(Dip/DipDir) & normals
        # ---------------------------------------------------------------------------

        # prefix + 1,2,3,Dip1,DipDir1,Dip2,DipDir2,Dip3,DipDir3,N1x,N1y,N1z,N2x,N2y,N2z,N3x,N3y,N3z

        if prefix:
            i = self.maxCol
            indexNew = {prefix+'1' : i+1, prefix+'2' : i+2, prefix+'3' : i+3}

            self.index.update(indexNew)
            self.maxCol+=3
            self.current = np.hstack([self.current, np.array(prinList)])

        if prefix and dipDir:
            i = self.maxCol
            indexNew.update({prefix+'Dip1' : i+1, prefix+'DipDir1' : i+2, \
                             prefix+'Dip2' : i+3, prefix+'DipDir2' : i+4, \
                             prefix+'Dip3' : i+5, prefix+'DipDir3' : i+6})

            self.index.update(indexNew)
            self.maxCol+=6
            self.current = np.hstack([self.current, np.array(dipsList).reshape((len(tArray),-1))])

        if prefix and normals:
            i = self.maxCol
            indexNew.update({prefix+'N1x' : i+1, prefix+'N1y' : i+2, prefix+'N1z' : i+3, \
                             prefix+'N2x' : i+4, prefix+'N2y' : i+5, prefix+'N2z' : i+6, \
                             prefix+'N3x' : i+7, prefix+'N3y' : i+8, prefix+'N3z' : i+9})

            self.index.update(indexNew)
            self.maxCol+=9
            self.current = np.hstack([self.current, np.array(vectList).reshape((len(tArray),-1))])

        # ~def addPrincipals(self, indices, prefix, dipDir=True, normals=True)

    def addNormals(self, prefix, indices='PBT'):
        """
        add normals for each plunge and trend in indices [trend1,plunge1,...]
        (uses extractTensor to get 6 components although not a tensor)

        arguments:
        -prefix: string, prefix for new columns
        -indices: string indicesMap() for column names
        """

        tArray = self.extractTensor(indices=indices, xyz=False)
        # dipsList = [[plunge1,trend1],...]
        dipsList = tArray.reshape(len(tArray), 3, 2)
        print (dipsList.shape)
        vectList = []

        for dip_dipDir in dipsList:

            # get normal vectors for each dip and direction
            vectList.append(getNormals(dip_dipDir))

        # ------------------------------------------------------------------------------
        # update self.current with normals, prefix + N1x,N1y,N1z,N2x,N2y,N2z,N3x,N3y,N3z

        i = self.maxCol
        indexNew = {prefix+'N1x' : i+1, prefix+'N1y' : i+2, prefix+'N1z' : i+3, \
                    prefix+'N2x' : i+4, prefix+'N2y' : i+5, prefix+'N2z' : i+6, \
                    prefix+'N3x' : i+7, prefix+'N3y' : i+8, prefix+'N3z' : i+9}

        self.index.update(indexNew)
        self.maxCol+=9
        self.current = np.hstack([self.current, np.array(vectList).reshape((len(tArray),-1))])

        # ~def addNormals(self, indices, prefix)

    def mapData(self, source, newIndex='mapData-1', overwrite=True, maxDist=False, \
                selectStep = False, fill=np.nan):
        """
        method to map data from source to self using kdTree
        - one column from np.array of shape (N, 4) into newIndex
        - all columns from xyzData class object with new indices from source.index
        - if selectStep is set, only map data for selectStep with stepCol in source
        """

        # source data
        # -----------
        if isinstance(source,np.ndarray):
            sourceData = source
            print (f"sourceData array: {sourceData.shape}")
        elif isinstance(source,xyzData):
            sourceData = source.current
            if self.stepCol:
                stepCol = self.index[self.stepCol]
            print (f"sourceData xyzData: {sourceData.shape}")
        else:
            return
        
        # target data
        # -----------
        targetData = self.current
        
        # kdTree
        # ------
        t0 = time.time()

        kdtree=KDTree(sourceData[:,0:3])
        dist,points=kdtree.query(targetData[:,0:3],1) # for ,2: points[i] becomes list

        if config.verbose:
            print (f"kdTree time mapping {source}: {time.time()-t0} seconds")

        if isinstance(source,np.ndarray): # map one column only

            if newIndex in self.index.keys():
                targetCol = self.index[newIndex]
            else:
                targetCol = self.maxCol+1
                self.index[newIndex]=targetCol
                overwrite=True
            if overwrite:
                # create new column with np.nan
                targetData = np.hstack([targetData,np.full([len(targetData),1],np.nan)])
    
                for i in range(len(points)):
                    if (maxDist is False) or (dist[i]<=maxDist):
                        targetData[i,targetCol]=sourceData[points[i],3]
                    else:
                        targetData[i,targetCol]=fill
                        
        # issue: check if overwrite works correctly

        elif isinstance(source,xyzData): # map all columns

            selfKeys = self.index.keys() # list of keys in self.index before mapping

            # create new columns with np.nan for each column in source
            for col in source.index.keys():
                # check if column exists in self, otherwise add new column with np.nan, 
                # skip x,y,z if naming is different to source (i.e. 'Easting','Northing','Elevation')
                if (col not in self.index.keys()) and (col not in source.xyzCol):
                    self.maxCol+=1
                    self.index[col]=self.maxCol

                    #adding new column with np.nan
                    targetData = np.hstack([targetData,np.full([len(targetData),1],np.nan)])

            # map data
            for i in range(len(points)):

                # selective fill based on step
                if not(selectStep) or (targetData[i,stepCol]==selectStep):

                    # each column from source data to be mapped to self
                    for col in source.index.keys():

                        # skip x,y,z columns in source.xyzCol
                        if col not in source.xyzCol:

                            # don't overwrite if column exists in self unless flag is set
                            if (col not in selfKeys) or overwrite:
                                if (maxDist is False) or (dist[i]<=maxDist):
                                    targetData[i,self.index[col]]=sourceData[points[i],source.index[col]]
                                else:
                                    targetData[i,self.index[col]]=fill  

        else:
            return

        if config.verbose:
            print (f"mapData time: {time.time()-t0} seconds")

        self.current = targetData
        self.bBox = array3D_BBox(self.current)

    # ~def mapData(self, source, newIndex='mapData-1')
    
    def associateSteps(self, seqFile=None, headers=['Start','Step'], newIndex='Step', overwrite=True):
        """
        method to associate time steps to self.current
        
        arguments:
        -seqFile: string, csv file name
        """

        endOfTime = np.datetime64('2100-01-01 00:00:00.000')
        optN = 0 # mark range to cycle_list
        
        if seqFile:
            self.seqFile = seqFile
        else:
            print (f"no sequence file")
            return
        
        seqData = pd.read_csv(seqFile)
        seqData.dropna(how='all', inplace=True)
        
        # test if self.current has date column
        
        if self.dateCol:
            dateIndex = self.index[self.dateCol]
        else:
            print (f"no date column")
            return
        if not((headers[0] in seqData) and (headers[1] in seqData)):
            print (f"no Start or Step column in {seqFile}")
            return
        
        self.filterNaN(self.dateCol)
        
        # create new column for steps in self.current if required
    
        if newIndex in self.index.keys():
            if overwrite:
                targetCol = self.index[newIndex]
            else:
                print (f"{newIndex} exists and overwrite==False")
                return
        else:
            self.maxCol += 1
            targetCol = self.maxCol
            self.index[newIndex]=targetCol
            self.current = np.hstack([self.current,np.zeros([len(self.current),1])])

        self.stepCol = newIndex

        # create list of dates plus 'endOfTime' to avoid repeat lookup
        
        dateList = np.append(np.array(list(seqData[headers[0]])[:],dtype='datetime64'),endOfTime)
        
        for n in range(len(seqData)):

            if (dateList[n] >= dateList[n+1]):
                print (f"Start dates in {seqFile} not in ascending order (line {n+3})")
                raise ValueError

        for n in range(len(seqData)-1):

            if seqData[headers[1]][n] >= seqData[headers[1]][n+1]:
                print (f"Steps in {seqFile} not in ascending order (line {n+3})")
                raise ValueError
            if seqData[headers[1]][n+1]-seqData[headers[1]][n]!=1:
                print (f"Warning: Step {n+2} missing in {seqFile} (line {n+3})")

        # loop over dates in self.current and associate steps from seqData
        t0 = time.time()

        for j in range(len(self.current[:,dateIndex])):

            if config.verbose:
                if (j % config.nHeartbeat == 0):
                    print (f"associateSteps: processing line {j} time {time.time()-t0} optN {optN}")

            # cyclic shift list of dates in seqData to use most recent first
            # --------------------------------------------------------------
            optRange = list(range(len(seqData))[optN:]) + list(range(len(seqData))[:optN])

            for n in optRange:

                if (dateList[n] <= self.current[j,dateIndex] < dateList[n+1]):
                    self.current[j,targetCol]=seqData[headers[1]][n]
                    optN = n
                    break
                
        # store current data in pData to keep Step column on reset

        self.pData = self.current

    # ~def associateSteps(self, seqFile=None, headers=['Start','Step'], newIndex='Step', overwrite=True)

    def uniqueSteps(self):
        """
        method to return unique steps in self.current
        
        returns: np.array()
        """

        if self.stepCol:
            stepCol = self.index[self.stepCol]
        else:
            print (f"no step column")
            return

        return np.unique(self.current[:,stepCol])

    def sortedArray(self, col='x', reverse=False):
        """
        method to sort self.current based on col
            
        arguments:
        -col:       integer, index / string: key for self.index[]
        -reverse:   boolean, True for descending order
        """
        # sort based on column index or string key
        # this code 100% written by GitHub Co-pilot
        # ----------------------------------------

        if isinstance(col,str):
            try:
                colStr=col
                col=self.index[colStr]
            except:
                print(f"sortedArray: '{colStr}' not found in self.index")
                return
        elif isinstance(col,int):
            if col>self.maxCol:
                print(f"sortedArray: col {col} not in self.current")
                return
        else:
            print(f"sortedArray: no valid col argument")
            return

        self.current = self.current[self.current[:,col].argsort()]
        if reverse:
            self.current = self.current[::-1]
                
        self.bBox = array3D_BBox(self.current)

        # ~def sortedArray(self, col='x', reverse=False)

    def addRValue(self, indices, newColumn='R-Value'):
        """
        method to add R value to self.current - see calcRvalue()

        arguments:
        -indices: list of three column names (i.e. 'P/B/T-Axis scale' or 'S1/2/3')
        """

        rColumn = np.zeros((self.current.shape[0], 1))

        for i in range(self.current.shape[0]):

            row = [self.current[i][self.index[x]] for x in indices ]
            rColumn[i] = calcRvalue(row)

        self.current = np.hstack([self.current, rColumn])

        indexNew = {newColumn : self.maxCol+1}
        self.index.update(indexNew)
        self.maxCol+=1
