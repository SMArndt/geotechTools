"""
stlGeom.py
"""

# ---------------------------------------------------------------------------
# imports
# ---------------------------------------------------------------------------

import numpy as np
import stl
import os
import time

# ---------------------------------------------------------------------------
# global variables
# ---------------------------------------------------------------------------

verbose = False

# ---------------------------------------------------------------------------
# functions
# ---------------------------------------------------------------------------

def array3D_BBox(a):
    """
    bounding box of np.array() of shape (N, 3 or more) with a[,0],[,1],[,2] = x,y,z
    
    returns tuple ((x0,y0,z0),(x1,y1,z1))
    """
    
    return ((min(a[:,0]), min(a[:,1]), min(a[:,2])), (max(a[:,0]), max(a[:,1]), max(a[:,2])))

# ---------------------------------------------------------------------------
# class stlGeom()
# ---------------------------------------------------------------------------

class stlGeom:

    def __init__(self, fileName=None):
        """
        constructor for stlGeom()
        """

        self.fileName = fileName
        self.pathName = None
        self.stlMesh = None
        self.bBox = None
        
        if self.fileName==None:
            pass
        else:
            self.read(self.fileName)

    # ~def __init__(self, fileName=None)

    def __str__(self):
        if self.stlMesh:
            if self.fileName:
                return f"File {self.fileName} {self.bBox} {len(self.stlMesh)} triangles"
            elif self.pathName:
                return f"Path {self.pathName} {self.bBox} {len(self.stlMesh)} triangles"
        else:
            return f"Path {self.pathName} {self.bBox} no Mesh"

    def read(self, fileName):
        """
        method to read stlGeom
        """

        self.fileName=fileName

        t0 = time.time()
        self.stlMesh = stl.mesh.Mesh.from_file(fileName)
        
        if verbose:
            print (f"time: {time.time()-t0} seconds")
            
        self.minMax()
            
    # ~read(self, fileName)

    def readPath(self, pathName, recursive=True, combine=False):
        """
        walk path, combine stl files and determine bounding box
        """

        self.fileName=None
        self.stlMesh = None
        self.bBox = None
        self.pathName=pathName
        oldData=[]

        t0 = time.time()

        for osFolder, osSubfolders, osFilenames in os.walk(pathName):
            if verbose:
                print(f"Reading {osFolder}")
            for osFilename in osFilenames:
                if combine and self.bBox:
                    oldData=self.stlMesh.data # save current data

                if os.path.splitext(osFilename)[1]=='.stl':
                    try:
                        self.stlMesh=stl.mesh.Mesh.from_file(os.path.join(osFolder,osFilename))
                        if verbose: print(f"... success reading {osFilename}")

                        # update bounding box
                        if self.bBox:
                            oldBox = self.bBox
                            self.minMax()
                            self.bBox = array3D_BBox(np.vstack((oldBox,self.bBox)))
                        else:
                            self.minMax()

                        # combine meshes
                        if combine and len(oldData) and self.bBox:
                            self.stlMesh = stl.mesh.Mesh(np.concatenate( \
                                [oldData,self.stlMesh.data]))
                    except:
                        print(f"..... error reading {osFilename}")
                        pass

            if not recursive:
                break

        if not combine:
            self.stlMesh = None
        
        if verbose:
            print (f"time: {time.time()-t0} seconds")
            print (self.bBox)

    # ~readPath(self, fileName)

    def minMax(self):
        """
        method to return bounding box of stl mesh using attributes min_ and max_
        """
        
        self.bBox = tuple(map(tuple,(self.stlMesh.min_, self.stlMesh.max_)))
        return self.bBox

    # ~minMax(self)

    def write(self, fileName, mode='ASCII'):
        """
        write stlGeom to file
        """

        if self.stlMesh:
            self.stlMesh.save(fileName, mode=stl.Mode.ASCII)
            if verbose:
                print (f"saved: {fileName}")
        else:
            print (f"no mesh to save")
            
    # ~write(self, fileName)