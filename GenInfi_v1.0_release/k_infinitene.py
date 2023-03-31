# k_infinitene: Create coordinates of generalized kekulene-based infinitene
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
#
########################################################################
# MORE DETAILS:
#
# Create coordinates of generalized kekulene-based infinitene:
# i.e., K-infinitene [h1,h2,h3,h4,h5,h6|k1,k2,k3,k4,k5,k6](d)
# subject to the following conditions:
# i) h1,...,h6, k1,...,k6 >=2 and are all POSITIVE INTEGER numbers;
# ii) (h1+...+h6) <= (k1+...+k6)
#     and when (h1+...+h6) == (k1+...+k6), let h6 <= k1
# iii)   -k1 <= d <= h6 and d must be an INTEGER number
# NOTE: One should take into account the chirality of the two kekulenes:
#    (1) Let's define the canonical notation of a kekulene [h1,...,h6]:
#        i.e., [h1,...,h6] is the smallest sequence
#          We can compare two sequences, (h1,...,h6) and (h1',...,h6'):
#            If h1<=h1' AND h2<=h2', ..., AND h6<=h6'
#            Then (h1,...,h6) <= (h1',...,h6')
#    (2) We can define the chirality of a kekulene [h1,...,h6]:
#        If the canonical notation of (h1,...,h6) is EQUAL TO the canonical
#        notation of its reverse fliplr(h1,...,h6), then [h1,...,h6] is
#        ACHIRAL (A);
#        Otherwise, [h1,...,h6] is CHIRAL:
#          If (h1,...,h6) < fliplr(h1,...,h6),
#          Then [h1,...,h6] is RIGHT (R);
#          If (h1,...,h6) > fliplr(h1,...,h6),
#          Then [h1,...,h6] is LEFT (L).
#     (3) When both [h1,...,h6] and [k1,...,k6] are chiral, we may have
#         4 different K-infinitenes:
#             (R,R), (L,L), (R,L), (L,R)
#         When [h1,...,h6] is chiral and [k1,...,k6] is achiral, we may
#         have 2 different K-infinitenes:
#             (R,A), (L,A)
#         When both [h1,...,h6] and [k1,...,k6] are achiral, we have only
#         1 K-infinitene:
#             (A,A)

import numpy as np
from kekulene import *
from write_xyz import *
from addH import *
from hmo import *
import glob

def zcoord_kinfi_top( xyz_kek1, RG, hexa_cc, ixOrder, zm, h6, k1, d ):
    # Cutting position to which all left rings (including ix_cut) on the 
    # contact side have z = -zm,  and all rings right to ix_cut on the contact 
    # side have z = zm (for kekulene [h1,...,h6])
    if d <= h6 - k1:
        ix_cut = d + k1 # Cut at east(E) edge
        print( 'Case I (d = %i; d <= h6 - k1)' % d )
    else:
        ix_cut = h6 # Cut at northeast(NE) edge
        print( 'Case II (d = %i; d > h6 - k1)' % d )

    if ix_cut < 0:
        print( 'ERROR: ix_cut (%i) must be a nonnegative integer' % ix_cut )
        exit(1)

    NAt = xyz_kek1.shape[0]

    # Ring indices of the bottom side:
    ix_bottom = ixOrder[ range( h6+1 ) ]

    # Determine z-coordinates of hexagon centers:
    NR = hexa_cc.shape[0]
    z_hexa_cc = np.zeros( (NR, 1) )
    # 1) For bottom side:
    Nhexa_bottom = len( ix_bottom )
    Nhexa_remain = NR - Nhexa_bottom
    zstep = ( 2 * zm )/( Nhexa_remain + 1 )
    z_hexa_cc[ range(ix_cut+1) ] = -abs( zm )
    z_hexa_cc[ range(ix_cut+1,Nhexa_bottom) ] = abs( zm )
    # 2) All other sides:
    #zval = np.arange( zm - zstep, -zm, -zstep )
    zval = np.zeros( (NR-Nhexa_bottom, 1) )
    for i in range( NR-Nhexa_bottom ):
        zval[i] =  zm - zstep + i*( -zstep )
    #z_hexa_cc[ range(Nhexa_bottom,NR) ] = zval.reshape( (len(zval),1) )
    z_hexa_cc[ range(Nhexa_bottom,NR) ] = zval

    
    # Determine z-coordinates of ATOMS:
    z_atoms = 999999 * np.ones( (NAt, ) )
    for j in range( NR ):
        rg = RG[ ixOrder[j] ]
        for k in range( len(rg) ):
            iat = rg[k]
            if z_atoms[iat] == 999999:
                z_atoms[iat] = z_hexa_cc[j]
            else: # For coronenes
                w = 0.5 # Empirical proportion
                z_atoms[iat] = w*z_hexa_cc[j] + (1-w)*z_atoms[iat]
    # Correct again for ix_cut ring in case of coronenes:
    rg = RG[ ixOrder[ix_cut] ]
    for k in range( len(rg) ):
        iat = rg[k]
        z_atoms[iat] = z_hexa_cc[ix_cut]

    # Update xyz_kek1:
    xyz_kek1[ :, 2 ] = z_atoms

    return xyz_kek1


def zcoord_kinfi_bottom( xyz_kek2, RG, hexa_cc, ixOrder, zm, h6, k1, d ):
    # Cutting position to which all left rings (including ix_cut) on the 
    # contact side have z = -zm,  and all rings right to ix_cut on the contact 
    # side have z = zm (for kekulene [k1,...,k6])
    if d <= h6 - k1:
        ix_cut = k1 # Cut at southeast(SE) edge
        print( 'Case I (d = %i; d <= h6 - k1)' % d )
    else:
        ix_cut = h6 - d # Cut at east(E) edge
        print( 'Case II (d = %i; d > h6 - k1)' % d )

    if ix_cut < 0:
        print( 'ERROR: ix_cut (%i) must be a nonnegative integer' % ix_cut )
        exit(1)

    NAt = xyz_kek2.shape[0]

    # Determine z-coordinates of hexagon centers:
    NR = hexa_cc.shape[0]
    z_hexa_cc = np.zeros( (NR, 1) )
    # 1) For top side:
    ix_top = range( ix_cut+1 )
    z_hexa_cc[ ix_top ] = abs( zm )
    # 2) All other sides:
    ix_rest = np.flipud( range(ix_cut+1,NR) )
    zstep = ( 2 * zm ) / ( len(ix_rest) + 1 )
    zval = np.zeros( (len(ix_rest), 1 ) )
    for i in range( len(ix_rest) ):
        zval[i] =  zm - zstep + i*( -zstep )
    #zval = np.arange( zm - zstep, -zm, -zstep )
    #z_hexa_cc[ ix_rest ] = zval.reshape( (len(zval),1) )
    z_hexa_cc[ ix_rest ] = zval

    # Determine z-coordinates of ATOMS:
    z_atoms = 999999 * np.ones( (NAt, ) )
    for j in range( NR ):
        rg = RG[ ixOrder[j] ]
        for k in range( len(rg) ):
            iat = rg[k]
            if z_atoms[iat] == 999999:
                z_atoms[iat] = z_hexa_cc[j]
            else: # For coronenes
                w = 0.5 # Empirical proportion
                z_atoms[iat] = w*z_hexa_cc[j] + (1-w)*z_atoms[iat]
    # Correct again for ix_cut ring in case of coronenes:
    rg = RG[ ixOrder[ix_cut] ]
    for k in range( len(rg) ):
        iat = rg[k]
        z_atoms[iat] = z_hexa_cc[ix_cut]

    # Update xyz_kek2:
    xyz_kek2[ :, 2 ] = z_atoms

    return xyz_kek2

def reorderRings_kek_top( coord, RG ):
    lm_ring = ringLinkage( RG )
    NR = len( RG )
    # Verify that each ring has only 2 neighbors so as to avoid dead loop
    # below:
    if len( np.where( np.sum( lm_ring, axis=0 ) > 2 )[0] ):
        print( np.where( np.sum( lm_ring, axis=0 ) > 2 )[0] )
        print( 'ERROR: Wrong connectivity between rings of cyclopolyarenes:' )
        print( np.sum( lm_ring, axis=0 ) )
        print( 'Rings are incorrectly identified as:' )
        for i in range( len(RG) ):
            print( '%2i:  ' % (i+1), end='' )
            for at in RG[i]:
                print( ' %i' % (at+1), end='' )
            print()
        print( 'ABORTED' )
        exit(1)

    hexa_cc = np.zeros( (NR, 3) )
    for j in range( NR ):
        hexa_cc[ j, : ] = coord[ RG[j], : ].mean( axis=0 )

    # Find the lowest rings:
    ix1 = np.where( abs( hexa_cc[:,1] - min(hexa_cc[:,1]) ) < 1E-3 )[0]

    # Find the lowest-leftmost ring:
    ix2 = ix1[ np.where( abs( hexa_cc[ix1,0] - min(hexa_cc[ix1,0]) ) \
            < 1E-3 )[0] ]

    # The first ring:
    ixOrder = np.array( ix2 )
    # Determine the second ring:
    ixRemain = np.setdiff1d( ix1, ixOrder )
    for j in ixRemain:
        if lm_ring[ ixOrder[0], j ]:
            ixOrder = np.append( ixOrder, j )
            break
    # The rest of the rings:
    ixRemain = np.setdiff1d( np.arange(NR), np.asarray(ixOrder) )
    while len(ixOrder) < NR:
        for j in ixRemain:
            if lm_ring[ ixOrder[len(ixOrder)-1], j ]:
                ixOrder = np.append( ixOrder, j )
                break
        ixRemain = np.setdiff1d( ixRemain, ixOrder )

    return ixOrder


def reorderRings_kek_bottom( coord, RG ):
    lm_ring = ringLinkage( RG )
    NR = len( RG )
    # Verify that each ring has only 2 neighbors so as to avoid dead loop
    # below:
    if len( np.where( np.sum( lm_ring, axis=0 ) > 2 )[0] ):
        print( np.where( np.sum( lm_ring, axis=0 ) > 2 )[0] )
        print( 'ERROR: Wrong connectivity between rings of cyclopolyarenes:' )
        print( np.sum( lm_ring, axis=0 ) )
        print( 'Rings are incorrectly identified as:' )
        for i in range( len(RG) ):
            print( '%2i:  ' % (i+1), end='' )
            for at in RG[i]:
                print( ' %i' % (at+1), end='' )
            print()
        print( 'ABORTED' )
        exit(1)

    hexa_cc = np.zeros( (NR, 3) )
    for j in range( NR ):
        hexa_cc[ j, : ] = coord[ RG[j], : ].mean( axis=0 )

    # Find the topmost rings:
    ix1 = np.where( abs( hexa_cc[:,1] - max(hexa_cc[:,1]) ) < 1E-3 )[0]

    # Find the leftmost-topmost ring:
    ix2 = ix1[ np.where( abs( hexa_cc[ix1,0] - min(hexa_cc[ix1,0]) ) \
            < 1E-3 )[0] ]

    # The first ring:
    ixOrder = np.array( ix2 )
    # Determine the second ring:
    ixRemain = np.setdiff1d( ix1, ixOrder )
    for j in ixRemain:
        if lm_ring[ ixOrder[0], j ]:
            ixOrder = np.append( ixOrder, j )
            break
    # The rest of the rings:
    ixRemain = np.setdiff1d( np.arange(NR), np.asarray(ixOrder) )
    while len(ixOrder) < NR:
        for j in ixRemain:
            if lm_ring[ ixOrder[len(ixOrder)-1], j ]:
                ixOrder = np.append( ixOrder, j )
                break
        ixRemain = np.setdiff1d( ixRemain, ixOrder )

    return ixOrder


def hexagon_centers_noclar( coord, RG, topOrBottom ):
    NR = len( RG )

    # Reorder rings according to the linkage between rings:
    if topOrBottom == 'top':
        ixOrder = reorderRings_kek_top( coord, RG )
    elif topOrBottom == 'bottom':
        ixOrder = reorderRings_kek_bottom( coord, RG )
    else:
        print( 'ERROR: topOrBottom should take the value of \'top\' ' \
              'or \'bottom\'' )
        exit(1)

    # Hexagon centers:
    hexa_cc = np.zeros( (NR, 3) )
    for j in range( NR ):
        hexa_cc[ j, : ] = coord[ RG[j], : ].mean( axis=0 )

    # Ordered hexagon centers according to ixOrder:
    ordhexa_cc = hexa_cc[ ixOrder, : ]

    return (hexa_cc, ordhexa_cc, ixOrder)



def k_infinitene( hh, kk, d, ifWrite=False ):
#  hh: an array containing h1, ..., h6 for kekulene [h1,h2,h3,h4,h5,h6]
#  kk: an array containing k1, ..., k6 for kekulene [k1,k2,k3,k4,k5,k6]
#  If ifWrite is TRUE, then write coordinate to an external *.xyz file
    #====================== Check validity of inputs ======================
    # Range of d:
    if d < -kk[0] or d > hh[5]:
        print( 'ERROR: d must lie in range %i <= d <= %i (d = %i)' % \
                (-kk[0], hh[5], d) )
        exit(1)
    #======================================================================

    #============================= Constants ==============================
    a = 1.42 # (in Angstrom) C--C bond length in benzene
    zm = 1.60 # A fixed distance between upper kekulene [h1,...h6] and 
              # upper kekulene [k1,...,k6]
    #======================================================================

    #======================================================================
    # For convenience of programming, we impose the constraint that
    # h6 must be greater than or equal to k1
    # So, we need to rename hh and kk
    hh0 = hh # Backup
    kk0 = kk # Backup
    d0 = d # Backup
    if_hk_switched = False
    if hh[5] < kk[0]:
        # Note:
        # (1) Make kk0[0] as the contact side
        hh = np.array( [ kk0[1], kk0[2], kk0[3], kk0[4], kk0[5], kk0[0] ], \
                dtype='i' )
        # Note:
        # (1) Make hh0[5] as the contact side
        # (2) Clockwise-->counterclockwise to maintain the chirality of infi.
        kk = np.array( [ hh0[5], hh0[4], hh0[3], hh0[2], hh0[1], hh0[0] ], \
                dtype='i' )
        d = -d0
        if_hk_switched = True
    #======================================================================

    # Initialize:
    Nhexa = sum( hh ) + sum( kk )
    print( 'K-infinitene [%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i):' % \
            (hh0[0],hh0[1],hh0[2],hh0[3],hh0[4],hh0[5], \
            kk0[0],kk0[1],kk0[2],kk0[3],kk0[4],kk0[5], d0) )
    print( '%i hexagons' % Nhexa )
    NAt = 4 * Nhexa
    print( '%i atoms' % NAt )
    xyz = np.zeros( (NAt, 3) )

    # I. Generate coordinates of two kekulenes, [h1,...,h6] and [k1,...,k6]:
    # (1) Reorder hh[] so that h6 becomes the contact (i.e. BOTTOM) side:
    h6 = hh[5] # backup hh[5]
    k1 = kk[0]
    hh = np.array( [ hh[2], hh[3], hh[4], hh[5], hh[0], hh[1] ], dtype='i' )
    xyz_kek1 = generalized_kekulene( hh, False, False )
    NAt1 = xyz_kek1.shape[0]
    #print( xyz_kek1 )


    # Determine all hexagons, sextets and their centers:
    # Ring list:
#    RG1 = getRings_from_coord( xyz_kek1 )
    RG1 = getRings_from_coord_corrected( xyz_kek1 )
    # Check validity of kekulene hh:
    if len( RG1 ) == 0:
        print( '!!! WARNING: wrong structure of kekulene [%i,%i,%i,%i,%i,%i]' \
                % hh )
        print( '!!! Aborted' )
        xyz = []
        return
    [ hexa_cc1, ordhexa_cc1, ixOrder1 ] = \
            hexagon_centers_noclar( xyz_kek1, RG1, 'top' )

    #   (2) No need to reorder kk[], as k1 is already the contact (TOP) side:
    xyz_kek2 = generalized_kekulene( kk, False, False )
    # Determine all hexagons, sextets and their centers:
    #RG2 = getRings_from_coord( xyz_kek2 )
    RG2 = getRings_from_coord_corrected( xyz_kek2 )
    # Check validity of kekulene kk:
    if len( RG2 ) == 0:
        print( '!!! WARNING: wrong structure of kekulene [%i,%i,%i,%i,%i,%i]' \
                % kk )
        print( '!!! Aborted' )
        xyz = []
        return
    [ hexa_cc2, ordhexa_cc2, ixOrder2 ] = \
            hexagon_centers_noclar( xyz_kek2, RG2, 'bottom' )

    # II. Generate (x,y)-coord. of joined infinitene [h1,...,h6|k1,...,k6](d):
    # Facts:
    #   For any given kekulene [h1,...,h6], all rings form a superhexagon,
    #   and the actual length of side h_i of the superhexagon is h_i*sqrt(3)*a
    #
    # 1. Translate kekulene [h1,...,h6] so that the h6-side coincides with the
    #    x-axis meanwhile the leftmost ring of h6-side is at the origin:
    xyz[ range(NAt1), : ] = moveCoord( xyz_kek1, -ordhexa_cc1[0,:] )

    # 2. Translate kekulene [k1,...,k6] so that the k1-side coincides with the
    #    x-axis meanwhile the leftmost ring of k1-side is at the origin:
    xyz[ range(NAt1,NAt), : ] = moveCoord( xyz_kek2, -ordhexa_cc2[0,:] )
    hexa_cc2 = moveCoord( hexa_cc2, -ordhexa_cc2[0,:] )

    # 3. Translate kekulene [k1,...,k6] along x-direction by the variable d:
    shift_x = d * np.sqrt(3) *a
    xyz[ range(NAt1,NAt), : ] = \
            moveCoord( xyz[ range(NAt1,NAt), : ], [ shift_x, 0., 0. ] )


    # III. Generate z-coordinates of the joined infinitene:
    #          [h1,h2,h3,h4,h5,h6|k1,k2,k3,k4,k5,k6](d):
    # 1. Kekulene [h1,h2,h3,h4,h5,h6]:
    xyz_kek1 = xyz[ range(NAt1), : ]
    xyz_kek1 = zcoord_kinfi_top( xyz_kek1, RG1, hexa_cc1, ixOrder1, \
            zm, h6, k1, d )
    xyz[ range(NAt1), : ] = xyz_kek1

    # 2. Kekulene [k1,k2,k3,k4,k5,k6]:
    xyz_kek2 = xyz[ range(NAt1,NAt), : ]
    xyz_kek2 = zcoord_kinfi_bottom( xyz_kek2, RG2, hexa_cc2, ixOrder2, \
            zm, h6, k1, d )
    xyz[ range(NAt1,NAt), : ] = xyz_kek2

    # Get rings of the K-infinitene:
    #RG = getRings_from_coord_corrected( xyz )
    lm = adjMat( xyz ) # Adjacent matrix
    RG = rings_kinfi( lm, NAt1, RG1, RG2, ixOrder1, ixOrder2, h6, k1, d )
    
    # IV. In case that the two kekulenes have been switched previously for
    # convenience of programming, now we switch them back:
    if if_hk_switched:
        print( 'h<->k switched and mirror structure is to be generated' )
        # Reflection about x-axis:
        for j in range( NAt ):
            xyz[ j, 2 ] = -xyz[ j, 2 ]

    # HMO total energy and eigenvalues:
    hmoSol = hmo( xyz ) # NOTE: xyz[] is so far only for CARBON atoms

    # Add hydrogen atoms:
    #[ xyz, elem ] = addH( xyz )
    [ xyz, elem ] = addH_corr( xyz, RG )

    # Check if number of H atoms is correct:
    NC = NAt
    NH = xyz.shape[0] - NAt
    print( '%i H atoms in K-infinitene '\
            '[%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i)' % \
            (NH, hh0[0],hh0[1],hh0[2],hh0[3],hh0[4],hh0[5], \
            kk0[0],kk0[1],kk0[2],kk0[3],kk0[4],kk0[5], d) )
    # Exceptional case::
    if NH != NC // 2:
        print( '!!! WARNING: Wrong number of H atoms found in K-infinitene '\
                '[%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i)' % \
                (hh0[0],hh0[1],hh0[2],hh0[3],hh0[4],hh0[5], \
                kk0[0],kk0[1],kk0[2],kk0[3],kk0[4],kk0[5], d) )
        print( '!!! %i H atoms are expected' % (NC // 2) )
        print( 'This infinitene structure should be adjusted manually' )

    # Get Huckel eigenvalues:
    hmoSol = hmo( xyz[range(NC),:] ) # xyz[] updated after calling addH_corr()
    print( '\nHuckel pi energy = %.6f |beta|' % hmoSol.Etot )

    # Output xyz to external *.gjf file:
    if ifWrite:
        fmt = 'xyz'
        outp = 'kinf_R%i-%i_%i_%i_%i_%i_%i-%i_%i_%i_%i_%i_%i_d%i.%s' % \
                (Nhexa,hh0[0],hh0[1],hh0[2],hh0[3],hh0[4],hh0[5], \
                kk0[0],kk0[1],kk0[2],kk0[3],kk0[4],kk0[5], d0, fmt)
        write_xyz( outp, xyz, elem, outp )

    return (xyz, elem, hmoSol.E)


# Enumerate infinitenes from two generalized kekulenes by assigning cyclic
# sequences of them, hh and kk, with all possible values of shift d.
#
# NOTE: For the input sequences hh and kk, the CONTACT sides are FIXED,
#       i.e., hh(6)|kk(1) OR kk(1)|hh(6) are fixed, while other indices
#       hh(i) and kk(j) can flip CYCLICALLY. Thus, the script will 
#       automatically enumerate all possible combinations of chiralities of
#       hh and kk.
#       Hence, this script is different from (more specific than) the
#       script enum_hk_all_k_infinitene().
#
def enum_hk_k_infinitene( hh, kk ):
    print( '='*60 )

    # I. Canonical notation of K-infinitene:
    #    [h1,h2,h3,h4,h5,h6|k1,k2,k3,k4,k5,k6](d)
    #    subject to the following conditions:
    #   i) h1,...,h6, k1,...,k6 >=1 and are all integer numbers;
    #  ii) (h1+...+h6) <= (k1+...+k6)
    #     and when (h1+...+h6) == (k1+...+k6), let h6 <= k1
    # iii) -k1 <= d <= h6 and d must be an integer number
    h_sum = sum( hh )
    k_sum = sum( kk )
    if h_sum > k_sum: # Swap hh <--> kk
        hh0 = hh
        kk0 = kk
        # (1) Make kk0(1) as the contact side:
        hh = [ kk0[1], kk0[2], kk0[3], kk0[4], kk0[5], kk0[0] ]
        # (1) Make hh0[5] as the contact side
        # (2) Clockwise-->counterclockwise to maintain the chirality of infi.
        kk = [ hh0[5], hh0[4], hh0[3], hh0[2], hh0[1], hh0[0] ]

    # II. Chirality of hh and kk:
    h_chir = hexa_chirality( hh )
    hh_arr = []
    if h_chir == 'R':
        hh_arr.append( hh ) # R
        hh_arr.append( [hh[4], hh[3], hh[2], hh[1], hh[0], hh[5]] ) # L
    elif h_chir == 'L':
        hh_arr.append( [hh[4], hh[3], hh[2], hh[1], hh[0], hh[5]] ) # R
        hh_arr.append( hh ) # L
    elif h_chir == 'A':
        hh_arr.append( hh ) # A

    k_chir = hexa_chirality( kk )
    kk_arr = []
    if k_chir == 'R':
        kk_arr.append( kk ) # R
        kk_arr.append( [kk[0], kk[5], kk[4], kk[3], kk[2], kk[1]] ) # L
    elif k_chir == 'L':
        kk_arr.append( [kk[0], kk[5], kk[4], kk[3], kk[2], kk[1]] ) # R
        kk_arr.append( kk ) # L
    elif k_chir == 'A':
        kk_arr.append( kk ) # A


    # III. Get all combinations of chirality of [hh|kk]:
    NMol = 0
    xyz_all = []
    HMOeig_all = []
    hkd_all = []
    for i1 in range( len(hh_arr) ):
        hh_ = hh_arr[ i1 ]
        hh_0 = hh_
        for i2 in range( len(kk_arr) ):
            kk_ = kk_arr[ i2 ]
            kk_0 = kk_
            if h_sum == k_sum and hh_0[5] > kk_0[0]: # Swap hh <--> kk
                # Make kk_0[0] as the contact side:
                hh_ = [ kk_0[1], kk_0[2], kk_0[3], kk_0[4], kk_0[5], kk_0[0] ]
                # Make hh0(6) as the contact side:
                kk_ = [ hh_0[5], hh_0[0], hh_0[1], hh_0[2], hh_0[3], hh_0[4] ]

            # Generate the K-infinitene:
            for d in range( -kk_[0], hh_[5]+1 ):
                print( '-'*50 )
                print( 'Generating', end='' )
                print( ' [%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i) ...' % \
                        (hh_[0],hh_[1],hh_[2],hh_[3],hh_[4],hh_[5], \
                        kk_[0],kk_[1],kk_[2],kk_[3],kk_[4],kk_[5], d) )
                NMol += 1
                [ xyz, elem, HMOeig ] = k_infinitene( hh_, kk_, d, False )
                xyz_all.append( xyz )
                HMOeig_all.append( HMOeig )
                hkd_all.append( [hh_[0],hh_[1],hh_[2],hh_[3],hh_[4],hh_[5], \
                        kk_[0],kk_[1],kk_[2],kk_[3],kk_[4],kk_[5], d ] )
                print( '-'*50 )

    print( '\n%i original structures generated' % NMol )
    HMOeig_all = np.array( np.round( HMOeig_all, 8 ) )
    HMOeig_all, ix = np.unique( HMOeig_all, axis=0, return_index=True )
    NMol_uniq = HMOeig_all.shape[0]
    print( '%i unique structures obtained' % NMol_uniq )
    # Update hkd_all and xyz_all:
    hkd_all = np.array( hkd_all, dtype='i' )[ ix, : ]
    xyz_uniq = []
    for i in range( NMol_uniq ):
        xyz_uniq.append( xyz_all[ ix[i] ] )
    print( '='*50 )

    return (hkd_all, xyz_uniq, elem, HMOeig_all)



# Enumerate infinitenes from two generalized kekulenes by assigning cyclic
# sequences of them, hh and kk, with all possible values of shift d.
#
# NOTE: The input sequences hh and kk are cyclically invariant and this
# script will automatically enumerate all nonequivalent [h1,h2,h3,h4,h5,h6]
# and [k1,k2,k3,k4,k5,k6] combinations and also take into account
# chirality.
# Hence, this script is more GENERAL than enum_hk_k_infinitene().
#
def enum_hk_all_k_infinitene( hh, kk ):
    hh_uniq = all_uniq_cyclicSeq( hh )
    kk_uniq = all_uniq_cyclicSeq( kk )

    hkd_all = []
    xyz_all = []
    HMOeig_all = []
    NMol = 0
    for i1 in range( hh_uniq.shape[0] ):
        hh_ = hh_uniq[ i1, : ]
        for i2 in range( kk_uniq.shape[0] ):
            kk_ = kk_uniq[ i2, : ]
            hkd_arr, xyz_arr, elem, HMOeig_arr = \
                    enum_hk_k_infinitene( hh_, kk_ )
            for j in range( len( xyz_arr ) ):
                NMol += 1
                hkd_all.append( hkd_arr[j,:] )
                xyz_all.append( xyz_arr[j] )

           # print( HMOeig_arr.shape )
           # np.append( HMOeig_all, [HMOeig_arr], axis=0 )
            if i1 == 0 and i2 == 0:
                HMOeig_all = HMOeig_arr
            else:
                HMOeig_all = np.vstack( (HMOeig_all, HMOeig_arr) )

    print( '\nA total of %i original structures enumerated' % NMol )
    HMOeig_all = np.array( np.round( HMOeig_all, 8 ) )
    HMOeig_all, ix = np.unique( HMOeig_all, axis=0, return_index=True )
    NMol_uniq = HMOeig_all.shape[0]
    print( 'A total of %i unique structures obtained' % NMol_uniq )

    # Write coordinates to xyz files:
    for i in range( NMol_uniq ):
        hkd = hkd_all[ ix[i] ]
        xyz = xyz_all[ ix[i] ]
        NR = len(xyz) // 6
        Nocc = 2*NR
        Ehmo = 2*np.sum( HMOeig_all[ i, range(Nocc) ], axis=0 )
        fmt = 'xyz'
        outp = 'kinf_R%i-%i_%i_%i_%i_%i_%i-%i_%i_%i_%i_%i_%i_d%i.%s' % \
                (NR,hkd[0],hkd[1],hkd[2],hkd[3],hkd[4],hkd[5], \
                hkd[6],hkd[7],hkd[8],hkd[9],hkd[10],hkd[11], hkd[12], fmt)
        title = 'K-infinitene [%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i) '\
                'EHMO: %.8f' % (hkd[0],hkd[1],hkd[2],hkd[3],hkd[4],hkd[5], \
                hkd[6],hkd[7],hkd[8],hkd[9],hkd[10],hkd[11], hkd[12], Ehmo)

        write_xyz( outp, xyz, elem, title )

    return (NMol, NMol_uniq)


# Enumerate all possible infinitenes of a given ring size, Nring, based on 
# two generalized kekulenes, taking into account of all possible shifts, d,
# and all combinations of charalities of the kekulenes.
#
def enum_Nring_all_k_infinitene( Nring ):
    # Enumerate all possible combinations of ring sizes of the two kekulenes:
    L1max = Nring // 2
    L_all = []
    NL_all = 0
    for L1 in range( 6, L1max+1 ):
        NL_all += 1
        L2 = Nring - L1
        L_all.append( [L1, L2 ] )

    # Enumerate all possible equiangular_hexagon for each of L_all[]:
    #counter_all = 0
    #counter_all_uniq = 0
    counter1 = 0
    for j in range( NL_all ):
        print( '#'*60 )
        L1 = L_all[j][0]
        L2 = L_all[j][1]
        counter1 += 1
        print( '(%i) RING SIZES OF KEKULENES: %2i + %2i' % (counter1, L1, L2) )
        print( '#'*60 )
        hh_arr = equiangular_hexagon( L1, False )
        kk_arr = equiangular_hexagon( L2, False )
    
        counter2 = 0
        for i1 in range( hh_arr.shape[0] ):
            hh = hh_arr[ i1, : ]
            for i2 in range( kk_arr.shape[0] ):
                kk = kk_arr[ i2, : ]
                counter2 += 1
                print( '#### (%i.%i):' % (counter1, counter2) )
                print( '####    kekulene 1: [%2i,%2i,%2i,%2i,%2i,%2i]' % \
                        (hh[0], hh[1], hh[2], hh[3], hh[4], hh[5]) )
                print( '####    kekulene 2: [%2i,%2i,%2i,%2i,%2i,%2i]' % \
                        (kk[0], kk[1], kk[2], kk[3], kk[4], kk[5]) )
                NMol, NMol_uniq = enum_hk_all_k_infinitene( hh, kk )
                #counter_all += NMol
                #counter_all_uniq += NMol_uniq
        print( '#'*60 )

    print( '\nAll [%i]K-infinitene structures have been successfully'
          ' generated' % Nring )
    Niso = len( glob.glob( 'kinf_R%i-*_*-*_*_d*.xyz' % Nring ) )
    print( 'A total of %i unique isomers (excluding enantiomers)' % Niso )
    print( '\nALL DONE' )
