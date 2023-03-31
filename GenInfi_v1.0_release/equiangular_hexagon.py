# equiangular_hexagon: Functions related to equiangular hexagons
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
import math

def mod_double( x, n ):
    PREC = np.finfo(float).eps*1E4

    x_int = math.floor( x + PREC )
    t = math.floor( x_int / n )
    r = x - n*t
    if r < PREC:
        r = 0
    return r


def equiangular_hexagon( L, ifVerbatim=True ):
    a0_min = math.ceil( L/2 )
    a0_max = math.ceil( 2/3*L ) - 1
    #print( a0_min, a0_max )

    if ifVerbatim:
        print( 'Equiangular hexagon with perimeter of %i:' % L )
        print( ' %2s %2s %2s %2s %2s %2s' % \
              ( 'H1', 'H2', 'H3', 'H4', 'H5', 'H6' ) )
        print( '-------------------' )

    h_arr = []
    counter = 0
    for a0 in range( a0_min, a0_max+1):
        a1 = L - a0

        b = ( 2*a1 - a0 ) / 3
        c = ( 2*a0 - a1 ) / 3
        #print( b, c )

        for y in np.arange( -(b-1), (b-1)/2 + 0.1, 1/3 ):
            for x in np.arange( 0, (b+y-1) + 0.1, 1/3 ):
                #print( 'y = %.4f x = %.4f' % ( y, x ) )

                h1 = b - 2*y
                h2 = c - x + y
                h3 = b + x + y
                h4 = c - 2*y
                h5 = b - x + y
                h6 = c + x + y
            
                # Rule out cases containing any noninteger number:
                if( mod_double(h1,1) != 0 or mod_double(h2,1) != 0 \
                    or mod_double(h3,1) != 0 or mod_double(h4,1) != 0 \
                    or mod_double(h5,1) != 0 or mod_double(h6,1) != 0 ):
                    continue
                else:
                    # Convert h1, h2, ..., h6 to integers:
                    h1 = round( h1 )
                    h2 = round( h2 )
                    h3 = round( h3 )
                    h4 = round( h4 )
                    h5 = round( h5 )
                    h6 = round( h6 )

                assert( h1 + h2 + h3 + h4 + h5 + h6 == L )


                # Reorder so as to make the minimum sequence:
                hh = np.zeros( (12, 6), dtype='i' )
                hh[0,:] = [ h1, h2, h3, h4, h5, h6 ]
                hh[1,:] = np.flipud( hh[0,:] )
                for k in range( 1, 6 ):
                    for j in range( 0, 5 ):
                        hh[2*k,j] = hh[2*k-2,j+1]
                    hh[2*k,5] = hh[2*k-2,0]
                    hh[2*k+1,:] = np.flipud( hh[2*k,:] )

                hh_min = hh[0,:]
                for k in range( 1, 12 ):
                    for j in range( 0, 6 ):
                        if hh[k,j] < hh_min[j]:
                            hh_min = hh[k,:]
                            break
                        elif hh[k,j] > hh_min[j]:
                            break

                #print( '%5i ' % (counter+1), end='' )
                #for k in range( len( hh_min ) ):
                #    print( ' %2i' % hh_min[k], end='' )
                #print()

                h_arr.append( hh_min )
                counter += 1

    h_arr = np.unique( h_arr, axis=0 )
    h_arr = np.flipud( h_arr )

    if ifVerbatim:
        for i in range( len( h_arr ) ):
            for k in range( len( h_arr[i] ) ):
                print( ' %2i' % h_arr[i][k], end='' )
            print()
        print( '-------------------' )
    
        print( 'Totally %i enumerated possibilities' % counter )
        print( 'Totally %i UNIQUE possibilities' % len( h_arr ) )

    return h_arr


def equiangular_hexagon_even( L, ifVerbatim=True ):
    if L % 1 != 0 or L < 6:
        print( 'ERROR: L must be an even positive integer (>=6) for '\
                'the perimeter of the hexagon' )
    if L % 2 != 0:
        print( 'ERROR: L must be a multiple of two for '\
                'the perimeter of the hexagon' )

    a0_min = math.ceil( L/2 )
    a0_max = math.ceil( 2/3*L ) - 1
    #print( a0_min, a0_max )

    if ifVerbatim:
        print( 'Equiangular hexagon with perimeter of %i:' % L )
        print( ' %2s %2s %2s %2s %2s %2s' % \
              ( 'H1', 'H2', 'H3', 'H4', 'H5', 'H6' ) )
        print( '-------------------' )

    h_arr = []
    counter = 0
    for a0 in range( a0_min, a0_max+1):
        a1 = L - a0

        b = ( 2*a1 - a0 ) / 3
        c = ( 2*a0 - a1 ) / 3
        #print( b, c )

        for y in np.arange( -(b-1), (b-1)/2 + 0.1, 1/3 ):
            for x in np.arange( 0, (b+y-1) + 0.1, 1/3 ):
                #print( 'y = %.4f x = %.4f' % ( y, x ) )

                h1 = b - 2*y
                h2 = c - x + y
                h3 = b + x + y
                h4 = c - 2*y
                h5 = b - x + y
                h6 = c + x + y
            
                # Rule out cases containing any odd number:
                if( mod_double(h1,2) != 0 or mod_double(h2,2) != 0 \
                    or mod_double(h3,2) != 0 or mod_double(h4,2) != 0 \
                    or mod_double(h5,2) != 0 or mod_double(h6,2) != 0 ):
                    continue
                else:
                    # Convert h1, h2, ..., h6 to integers:
                    h1 = round( h1 )
                    h2 = round( h2 )
                    h3 = round( h3 )
                    h4 = round( h4 )
                    h5 = round( h5 )
                    h6 = round( h6 )

                assert( h1 + h2 + h3 + h4 + h5 + h6 == L )


                # Reorder so as to make the minimum sequence:
                hh = np.zeros( (12, 6), dtype='i' )
                hh[0,:] = [ h1, h2, h3, h4, h5, h6 ]
                hh[1,:] = np.flipud( hh[0,:] )
                for k in range( 1, 6 ):
                    for j in range( 0, 5 ):
                        hh[2*k,j] = hh[2*k-2,j+1]
                    hh[2*k,5] = hh[2*k-2,0]
                    hh[2*k+1,:] = np.flipud( hh[2*k,:] )

                hh_min = hh[0,:]
                for k in range( 1, 12 ):
                    for j in range( 0, 6 ):
                        if hh[k,j] < hh_min[j]:
                            hh_min = hh[k,:]
                            break
                        elif hh[k,j] > hh_min[j]:
                            break

                #print( '%5i ' % (counter+1), end='' )
                #for k in range( len( hh_min ) ):
                #    print( ' %2i' % hh_min[k], end='' )
                #print()

                h_arr.append( hh_min )
                counter += 1

    h_arr = np.unique( h_arr, axis=0 )
    h_arr = np.flipud( h_arr )

    if ifVerbatim:
        for i in range( len( h_arr ) ):
            for k in range( len( h_arr[i] ) ):
                print( ' %2i' % h_arr[i][k], end='' )
            print()
        print( '-------------------' )
    
        print( 'Totally %i enumerated possibilities' % counter )
        print( 'Totally %i UNIQUE possibilities' % len( h_arr ) )

    return h_arr


def equiangular_hexagon_coord( h ):
# h: 6 positive intergers indicating the sides of the equiangular hexagon
# which can be generated from h_arr = equiangular_hexagon( L )
    s = np.sqrt(3)/2

    xy = np.zeros( (6, 2) )
    xy[1,:] = [ h[0], 0 ]
    xy[2,:] = [ h[0] + h[1]/2, -s*h[1] ]
    xy[3,:] = [ h[0] + h[1]/2 - h[2]/2, -s*(h[1]+h[2]) ]
    xy[4,:] = [ h[0] + h[1]/2 - h[2]/2 - h[3], -s*(h[1]+h[2]) ]
    xy[5,:] = [ h[0] + h[1]/2 - h[2]/2 - h[3] - h[4]/2, -s*(h[1]+h[2]-h[4]) ]

    return xy


def compSeq( hh, kk ):
    N = len(hh)
    b = 0
    for j in range( N ):
        if hh[j] == kk[j]:
            continue
        else:
            b = hh[j] - kk[j]
            return b

    return b


def minSeq( hh ):
    hh = np.array(hh)
    min_seq = hh
    N = len( hh )
    for j in range( 1, N ):
        seq = hh[ np.concatenate( (range(j,N), range(j)) ) ]
        if compSeq( seq, min_seq ) < 0:
            min_seq = seq

    return min_seq


def hexa_chirality( hh ):
    # Minimal sequence of hh:
    min_seq1 = minSeq( hh )
    min_seq2 = minSeq( np.flipud( hh ) )

    b = compSeq( min_seq1, min_seq2 )
    if b == 0:
        chirality = 'A'
    elif b < 0:
        chirality = 'R'
    else:
        chirality = 'L'

    return chirality


# Get all distinct sequence from a given cyclic array, hh[]
# NOTE: the chirality of hh[] is FIXED during the enumeration!
def all_uniq_cyclicSeq( hh ):
    hh = np.array( hh )
    N = len( hh )
    hh_all = [ hh ]
    for j in range( 1, N ):
        hh_all.append( hh[ np.concatenate( (range(j,N), range(j)) ) ] )

    hh_uniq = np.unique( hh_all, axis=0 )
    return hh_uniq


def assert_equiangular_hexagon( hh ):
    if hh[0] + hh[1] != hh[3] + hh[4]:
        return False
    if hh[1] + hh[2] != hh[4] + hh[5]:
        return False
    if hh[2] + hh[3] != hh[5] + hh[0]:
        return False

    return True
