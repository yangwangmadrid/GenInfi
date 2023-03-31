# clarene: Create coordinates of clarene <h1,h2,h3,h4,h5,h6>
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
from equiangular_hexagon import *
from transform import *
from write_xyz import *
from addH import *

pi = np.pi
cos = np.cos
sin = np.sin
    

def clarene( h, ifAddH=False, ifWrite=False ):
# h: 6 positive intergers indicating the sides of the equiangular hexagon
# which can be generated from h_arr = equiangular_hexagon( L )

    if len(h) != 6:
        print( 'h must be an array containing 6 positive integers' )
        exit(1)

    for i in range(6):
        if h[i] % 2 != 0 or h[i] <= 0:
            print( 'h[%i] must be an even positive integer' % i )
            exit(1)


    #========== Constants ==========
    a = 1.42 # (in Angstrom) C--C bond length in benzene
    s = np.sqrt(3)/2

    # I. Create the central points of all hexagons:
    Nhexa = sum( h )
    print( '%i hexagons in clarene <%i,%i,%i,%i,%i,%i>' % \
            (Nhexa, h[0],h[1],h[2],h[3],h[4],h[5]) )
    NAt = 4 * Nhexa
    print( '%i C atoms in clarene <%i,%i,%i,%i,%i,%i>' % \
            (NAt, h[0],h[1],h[2],h[3],h[4],h[5]) )


    xy_h = np.zeros( (Nhexa, 2) ) # 2D coordinates of hexagon centers
    xy = np.zeros( (NAt, 2) )

    # Only work out 2D coordinates, (x,y), for the molecule is planar (z = 0)
    # 1. Get all centers of the 6 corner hexagons:
    xy_h[0:6,:] = 3/2*a * equiangular_hexagon_coord( h )
    xy_h[0:6,:] = centerize( xy_h[0:6,:] )
    Ncorner = 6
    counter = Ncorner

    # 2. Generate all intermediate points along each side of superhexagon:
    for k in range( 6 ):
        # Side k: all intermediate points between corners k---(k+1):
        if k < 5:
            k1 = k + 1
        else:
            k1 = 0
        nsextet = h[k] // 2
        dv = ( xy_h[k1,:] - xy_h[k,:] ) / nsextet
        for i in range( nsextet-1 ):
            xy_h[ counter + i, : ] = xy_h[ k, : ] + dv * (i+1)
        counter += nsextet-1 

    Nsex = Nhexa // 2
    assert( counter == Nsex )


    # Add non-sextet hexagons:
    # -- First, get linkage matrix between all sextet rings:
    lm_h = np.zeros( (Nsex, Nsex),dtype='i' )
    for j in range( Nsex ):
        for k in range( j+1, Nsex ):
            r = np.linalg.norm( xy_h[j,:] - xy_h[k,:] ) / (np.sqrt(3)*a)
            if r < 2 + 1E-3:
                lm_h[j,k] = 1
                lm_h[k,j] = 1
    assert( np.array_equal( np.sum(lm_h,axis=0), 2*np.ones((1,Nsex))[0] ) )
    # -- Then, get left-neighbor-sextet of each sextet ring:
    v_shift = np.zeros( (Nsex, 2) ) # Shift vector from each sextet to nonsextet
    cc = xy_h[0:Nsex,:].mean( axis=0 ) # Center of the molecule
    for j in range( Nsex ):
        ix_nb = np.where( lm_h[j,:] == 1 )[0]
        v1 = xy_h[ ix_nb[0], : ] - xy_h[ j, : ]
        v2 = xy_h[ ix_nb[1], : ] - xy_h[ j, : ]
        v0 = xy_h[ j, : ] - cc
        vn1 = np.cross( np.append( v1, 0. ), np.append( v0, 0. ) )
        vn2 = np.cross( np.append( v2, 0. ), np.append( v0, 0. ) )
        if vn1[2] < 0 and vn2[2] > 0:
            # Rotate by 30 degree (counter-clockwise):
            v_shift3 = rotateCoord( np.append( v1, 0. ), [ 0, 0, 1 ], 30 )
            v_shift[j,:] = v_shift3[ 0:2 ] / np.sqrt(3)
        elif vn1[2] > 0 and vn2[2] < 0:
            # Rotate by 30 degree (counter-clockwise):
            v_shift3 = rotateCoord( np.append( v2, 0. ), [ 0, 0, 1 ], 30 )
            v_shift[j,:] = v_shift3[ 0:2 ] / np.sqrt(3)
        else:
            print( 'ERROR: Impossible' )
            exit(1)
    # -- Last, add nonsextet rings:
    for j in range( Nsex ):
        xy_h[ Nsex + j, : ] = xy_h[ j, : ] + v_shift[ j, : ]


    # II. Create the corner points of each hexagon:
    # 1. Generate the first set of six radial vectors centered at the central
    #    point of each hexagon, which is already generated above:
    t = np.arange( 6 ) * 2*pi / 6
    v_n = np.column_stack( (cos(t)*a, sin(t)*a) )

    iC = 0
    for ih in range( Nhexa ):
        # Translation by radial vectors:
        for i in range( 6 ):
            xy_tmp = xy_h[ ih, : ] + v_n[ i, : ]
            # Rule out already generated C atoms:
            nC = iC
            flag = True
            for k in range( nC ):
                if np.linalg.norm( xy_tmp - xy[ k, : ] ) < 1E-3:
                    flag = False # Found duplicated
                    break
            if flag: # Not duplicated
                xy[ iC, : ] = xy_tmp
                iC += 1

    # Centerize:
    xyz = np.zeros( (NAt, 3) )
    xyz[ :, 0:2 ] = centerize( xy )

    # Add hydrogen atoms:
    if ifAddH:
        [ xyz, elem ] = addH( xyz )
    else:
        elem = []
        for i in range( NAt ):
            elem.append( 'C' )


    # Write to file:
    if ifWrite:
        fmt = 'xyz'
        outp = 'clr_R%i-%i_%i_%i_%i_%i_%i.%s' % \
                (Nhexa, h[0],h[1],h[2],h[3],h[4],h[5], fmt)
        write_xyz( outp, xyz, elem, outp )

    return xyz



def enum_clarene( Nring ):
    hh_arr = equiangular_hexagon_even( Nring )
    Niso = hh_arr.shape[0]
    print( '%i isomers of [%i]clarenes' % (Niso, Nring) )

    for i in range( Niso ):
        print()
        print( 'Generating isomer %i ...' % (i+1) )
        clarene( hh_arr[i,:], True, True )

