# transform: Functions for manipulation of Cartesian coordinates
#
# Copyright (C) 2023 Yang Wang
#
# Author:
#     Yang Wang, 2023
#     ( yangwang@yzu.edu.cn; yangwang2008@gmail.com )
#
# --------------------------------------------------------------------
#  This file is part of GenInfi.
#   
#  GenInfi is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#   
#  GenInfi is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#   
#  You should have received a copy of the GNU General Public License
#  along with GenInfi.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------
#
#  Created on March, 2023
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#    Mar 26, 2023:
#    -  Completed all basic functions in this source file
#

import numpy as np

def cart2pol( x, y ):
    rho = np.sqrt( x*x + y*y )
    phi = np.arctan2( y, x )
    return (rho, phi)

def pol2cart( rho, phi ):
    x = rho * np.cos( phi )
    y = rho * np.sin( phi )
    return (x, y)

# Centerize
def centerize( xyz ):

    NAt = xyz.shape[0]
    cc = xyz.mean( axis=0 )
    for i in range( NAt ):
        xyz[i,:] = xyz[i,:] - cc

    return xyz

# Translation operation on a set of coordinates
def  moveCoord( xyz, tv ):
# tv <vector>: translational vector

    xyz1 = np.zeros( xyz.shape )
    for k in range( len(tv) ):
        xyz1[ :, k ] = xyz[ :, k ] + tv[k]

    return xyz1

# Rotation operation on a set of coordinates
def rotateCoord( xyz, tv, ang ):
#   tv <vector>: directional vector of the rotation axis
#   ang <number>: rotation angle in degree
    PREC = np.finfo(float).eps*1E4

    xyz = np.array( xyz )
    tv /= np.linalg.norm( np.array(tv) )
    ha = ang*np.pi/360 # half of rotation angle

    # case for C_2 rotation
    if abs( ang - 180 ) <= PREC:
        xyz1 = -xyz + 2*(xyz@tv)*tv
        return xyz1

    xyz1 = np.zeros( xyz.shape )
    Q = xyz - np.outer (xyz@tv, tv ) 
    if xyz.ndim == 1:
        N = 1
    else:
        N = xyz.shape[0]
    for k in range( N ):
        if Q.ndim == 1:
            Qk = Q
        else:
            Qk = Q[k,:]
        if np.linalg.norm( Qk ) == 0:
            xyz1[k,:] = xyz[k,:]
            continue

        u = np.cross( Qk, tv ) / ( np.tan(ha) + PREC )
        v = u - Qk
        v /= np.linalg.norm(v)
        if xyz.ndim == 1:
            xyz1 = xyz + 2*np.sin(ha)*np.linalg.norm(Qk) * v
        else:
            xyz1[k,:] = xyz[k,:] + \
                2*np.sin(ha)*np.linalg.norm(Qk) * v

    return xyz1 


