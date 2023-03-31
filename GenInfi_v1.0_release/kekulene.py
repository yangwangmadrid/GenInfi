# kekulene: Create coordinates of generalized kekulene [h1,h2,h3,h4,h5,h6]
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
    

def generalized_kekulene( h, ifAddH=False, ifWrite=False ):
# h: 6 positive intergers indicating the sides of the equiangular hexagon
# which can be generated from h_arr = equiangular_hexagon( L )

    if len(h) != 6:
        print( 'h must be an array containing 6 positive integers' )
        exit(1)

    for i in range(6):
        if h[i] % 1 != 0 or h[i] <= 0:
            print( 'h[%i] must be a positive integer' % i )
            exit(1)


    #========== Constants ==========
    a = 1.42 # (in Angstrom) C--C bond length in benzene
    s = np.sqrt(3)/2

    # I. Create the central points of all hexagons:
    Nhexa = sum( h )
    print( '%i hexagons in kekulene [%i,%i,%i,%i,%i,%i]' % \
            (Nhexa, h[0],h[1],h[2],h[3],h[4],h[5]) )
    NAt = 4 * Nhexa
    print( '%i C atoms in kekulene [%i,%i,%i,%i,%i,%i]' % \
            (NAt, h[0],h[1],h[2],h[3],h[4],h[5]) )


    xy_h = np.zeros( (Nhexa, 2) ) # 2D coordinates of hexagon centers
    xy = np.zeros( (NAt, 2) )

    # Only work out 2D coordinates, (x,y), for the molecule is planar (z = 0)
    # 1. Get all centers of the 6 corner hexagons:
    xy_h[0:6,:] = np.sqrt(3)*a * equiangular_hexagon_coord( h )
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
        dv = ( xy_h[k1,:] - xy_h[k,:] ) / h[k]
        for i in range( h[k]-1 ):
            xy_h[ counter + i, : ] = xy_h[ k, : ] + dv * (i+1)
        counter += h[k]-1 

    assert( counter == Nhexa )


    # II. Create the corner points of each hexagon:
    # 1. Generate the first set of six radial vectors centered at the central
    #    point of each hexagon, which is already generated above:
    t = np.arange( 0.5, 6.0 ) * 2*pi / 6
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
    NC = NAt # Number of C atoms

    # Add hydrogen atoms:
    if ifAddH:
        [ xyz, elem ] = addH( xyz )
    else:
        elem = []
        for i in range( NAt ):
            elem.append( 'C' )
    NH = len(elem) - NC # Number of C atoms
    if ifAddH and NC != NH*2:
        print( 'ERROR: Unable to construct the generalized kekulene '\
                '[%i,%i,%i,%i,%i,%i]' % (h[0],h[1],h[2],h[3],h[4],h[5]) )
        print( '  Possible reasons:')
        print( '    - The given side lengths form an unequiangular hexagon' )
        print( '    - Fused rings are found in this generalized kekulene' )
        print( 'ABORTED' )
        return []


    # Write to file:
    if ifWrite:
        fmt = 'xyz'
        outp = 'kek_R%i-%i_%i_%i_%i_%i_%i.%s' % \
                (Nhexa, h[0],h[1],h[2],h[3],h[4],h[5], fmt)
        write_xyz( outp, xyz, elem, outp )

    return xyz



def enum_generalized_kekulene( Nring ):
    hh_arr = equiangular_hexagon( Nring )
    Niso = hh_arr.shape[0]
    print( '%i isomers of [%i]kekulenes' % (Niso, Nring) )

    for i in range( Niso ):
        print()
        print( 'Generating isomer %i ...' % (i+1) )
        generalized_kekulene( hh_arr[i,:], True, True )

