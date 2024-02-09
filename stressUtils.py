"""
stressUtils.py - Copyright 2024 S.M.Arndt, Cavroc Pty Ltd
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
import math

def rot_x(theta):
    """
    3D rotation matrix theta, positive (right hand rule) around x axis
    """

    c,s = math.cos(theta), math.sin(theta)

    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def rot_y(theta):
    """
    3D rotation matrix theta, positive (right hand rule) around y axis
    """

    c,s = math.cos(theta), math.sin(theta)

    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def rot_z(theta):
    """
    3D rotation matrix theta, positive (right hand rule) around z axis
    """

    c,s = math.cos(theta), math.sin(theta)

    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

def unpackStress(s):
    """
    returns np.array(3,3) from stress s = [Sxx,Syy,Szz,Sxy,Sxz,Syz]
    """
    return np.array([[s[0],s[3],s[4]],[s[3],s[1],s[5]],[s[4],s[5],s[2]]])

def packStress(T):
    """
    returns stress s = [Sxx,Syy, Szz, Sxy,Sxz,Syz] from np.array(3,3)
    """
    return [T[0,0],T[1,1],T[2,2],T[0,1],T[0,2],T[1,2]]

def getCartesianStress(Principals):
    """
    Principals are a list of Stress value, dip and trend for S1, S2, S3
    
    returns Cartesian stress tensor as np.array(3,3)
    """
    
    Cartesian = np.zeros((3,3))                      
    for(value,dip,trend) in Principals:
    
        T=np.array([[0,0,0],[0,value,0],[0,0,0]])

        # rotation matrix R = Rz(-trend) * Rx(-dip)
        R = np.dot(rot_z(math.radians(-trend)), rot_x(math.radians(-dip)))
        # rotated tensor R * T * R.T
        T = np.dot(np.dot(R, T), R.T)
    
        Cartesian += T
    
    return Cartesian

def getPrincipalStress(Cartesian):
    
    eigenvalues, eigenvectors = np.linalg.eig(Cartesian)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    return eigenvalues, eigenvectors

def getStressOrientation(eigenvectors):
    """
    returns plunge and trend for each principal stress
    """
    
    orientations=[]
    for v in eigenvectors.T:
        # calculate plunge and trend
        sign = np.sign(v[2])
        plunge = math.degrees(math.asin(v[2]*sign))
        trend = (math.degrees(math.atan2(v[0]*sign,v[1]*sign))+180)%360

        orientations.append([plunge,trend])
    
    return orientations

def getPlaneOrientation(normal):
    """
    returns dip and direction for normal vector
    """
    
    # calculate dip and direction
    sign = np.sign(normal[2])
    dip = math.degrees(math.acos(normal[2]*sign))
    ddir = (math.degrees(math.atan2(normal[0]*sign,normal[1]*sign))+360)%360

    return (dip,ddir)