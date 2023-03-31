# topo: Topology and graph calculations
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
#  Updates:
#    March, 2023:
#    -  Added some functions for GenInfi
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Feb 4, 2021:
#    - Added function readRings()
#
#    Dec 27, 2020:
#    - Added functions to get connectivity and neighbors of rings
#
#    Dec 26, 2020:
#    - Added functions to determine the rings
#
#    May 31, 2020:
#    - Fixed a bug in function xyzNonH() for choosing non-hydrogen atoms
#

import numpy as np

# Only choose non-hydrogen atoms from a set of xyz coordinates:
#
#  Return:
#           xyz1: Cartesian coordinates of non-hydrogen atoms
#             at: List of indices of atoms ( starting from 0 )
#
def xyzNonH( elem, xyz ):
    at = []
    iat = -1
    for e in elem:
        iat += 1
        ecap = e.upper()
        if ecap == '1' or ecap == 'H':
            continue
        at.append( iat )
    xyz1 = xyz[ at, : ]
    return xyz1, at
# enddef xyzNonH()


# Adjacency matrix from a set of xyz coordinates:
#   xyz: a ndarray storing a set of Cartesian coordinates
def adjMat( xyz ):

    MAX_BONDLEN = 1.80

    N = xyz.shape[0] # Number of atoms
    A = np.zeros( ( N, N ), dtype='i' )

    for i in range( N ):
        for j in range( i+1, N ):
            r = np.linalg.norm( xyz[i,:] - xyz[j,:] )
            if r <= MAX_BONDLEN:
                A[i,j] = 1
                A[j,i] = 1

    return A
# enddef adjMat()


# Get bond list from an adjacent matrix
#
#    A: Adjacent matrix
#   at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of the two atoms
def bonds( A, at=[] ):

    N = len( A )
    B = []
    for i in range( N ):
        for j in range( i+1, N ):
            if A[i,j] == 1:
                B.append( np.array([ i, j ]) )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return B

    B0 = []
    for bnd in B:
        B0.append( np.array([ at[bnd[0]], at[bnd[1]] ]) )
    return B0

# enddef bonds()


# Get neighbors of each atom
#
#    A: Adjacent matrix
#   at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of all neighbors of 
#   each atom
def neighbors( A, at=[] ):
    N = len( A )
    nblist = []
    for i in range( N ):
        ix = [ j for j, a in enumerate(A[i,:]) if a > 0 ] 
        nblist.append( ix )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return nblist

    nblist0 = []
    for nb in nblist:
        x = [ at[a] for a in nb ]
        nblist0.append( x )
    return nblist0
# enddef neighbors()


# Remove EXACTLY the bonds specified by list bnd[] from the whole bond list B[]:
# NOTE: The two atomic indices in bnd[] are in ascending order, and so are the
# two atomic indices in each element of list B[]
def removeExactBonds( B, bnd ):

    Bnew = []
    for bd in B:
        if np.array_equal( bd, bnd ):
            continue
        Bnew.append( bd )

    return Bnew
# enddef removeExactBonds()


# Remove all bonds INVOLVING atoms in bnd[] from bond list:
# NOTE: Note the difference between removeBonds() and removeExactBonds()
def removeBonds( B, bnd ):
 
    # Number of bonds in the list:
    NB = len( B )
 
    Bnew = []
    for i in range(NB):
        if B[i][0] == bnd[0] or B[i][1] == bnd[0] or \
           B[i][0] == bnd[1] or B[i][1] == bnd[1]:
            continue
        Bnew.append( B[i] )
 
    return Bnew
# enddef removeBonds()


# Remove all bonds INVOLVING given atoms:
def removeBonds_by_Atoms( B, Atoms ):
 
    # Number of bonds in the list:
    NB = len( B )
    # Number of given atoms:
    if not isinstance( Atoms, list ):
        Atoms = [ Atoms ]
    NA = len( Atoms )
 
    Bnew = []
    for i in range( NB ):
        flag = False
        for j in range( NA ):
            # Found the bond involving atom j:
            if B[i][0] == Atoms[j] or B[i][1] == Atoms[j]:
                flag = True
                break
        if flag:
            continue
        Bnew.append( B[i] )
 
    return Bnew
# enddef removeBonds_by_Atoms()


# Remove all involved atoms in bnd[] from atom list:
def removeAtoms_from_Bond( A, bnd ):
 
    # Number of atoms in the list:
    NAt = len( A )
 
    Anew = []
    for i in range(NAt):
        if A[i] == bnd[0] or A[i] == bnd[1]:
            continue
        Anew.append( A[i] )
 
    return Anew
# enddef removeAtoms_from_Bonds()


# Check if a ring is monocyclic
#
#  rg: a ring specified by a list of indices of atoms (starting from 0)
#  rglist: List of rings to store the results of determination; For initial
#          call, set it an empty list
def isMonocyclic( rg, nblist ):
    nrg = len( rg )

    for j in range( nrg ):
        iat = rg[j]
        if j == 0:
            iat_1 = rg[nrg-1]
            iat1 = rg[1]
        elif j == nrg-1:
            iat_1 = rg[j-1]
            iat1 = rg[0]
        else:
            iat_1 = rg[j-1]
            iat1 = rg[j+1]
        nb = nblist[ iat ]
    
        for a in nb:
            if a == iat_1 or a == iat1:
                continue
            # The rest of neighbors of iat:
            if a in rg:
                return False

    return True
# enddef isMonocyclic()


# Canonical numbering of a ring. The rule is as follows:
#   - 1. Keep the first two numbers as small as possible.
#   - 2. Being clockwise or counterclockwise does not matter.
def canonicalRing( rg ):
    nrg = len( rg )
    ix = rg.index( min(rg) )
    if ix == 0:
        at_1 = rg[ nrg-1 ]
    else:
        at_1 = rg[ ix-1 ]

    if ix == nrg-1:
        at1 = rg[0]
    else:
        at1 = rg[ ix+1 ]

    num = list( range( nrg ) )
    loop = num[ ix+1: ] + num[ 0 : ix ]
    if at_1 < at1: # Reverse
        loop.reverse()

    loop = [ ix ] + loop

    return [ rg[i] for i in loop ]
# enddef canonicalRing()


# Get rings from a given atom
#
#  rglist: List of rings to store the results of determination; For initial
#          call, set it an empty list
#  rg: Current ring; For initial call, rg = [iat], where iat is the index
#      of the given atom
#  iat: Index of the given atom, which must be from 0, 1, 2, ..., NAt
#       NOTE: It may not the actual atomic index in the whole molecule that
#       contains hydrogens; it is the appearance index among non-H atoms
#  nblist: List, each entry indicating indices (starting from 0) of all 
#          neighbors of each atom
#
# Return:
#   List, each entry indicating indices (starting from 0) of each ring
def rings_from_an_atom( rglist, rg, nblist, iat ):

    # Number of atoms in current rg[]:
    nrg = len( rg )
    nb = nblist[ iat ]

    for at_next in nb:
        # Check if at_next is already in rg[]:
        #ix = [ i for i, a in enumerate(rg) if a == at_next ]
        #ix = rg.index( at_next )
        rg_curr = rg.copy()

        # Avoid large rings:
        if len( rg_curr ) > 7:
            continue

        if at_next in rg:
            #print( at_next, rg )
            #input()
            ix = rg.index( at_next )
            # Backward atom:
            if ix == nrg-2:
                continue
            # Ring complete:
            elif rg.index( at_next ) == 0:
                if isMonocyclic( rg_curr, nblist ):
                    rg_curr = canonicalRing( rg_curr )
                    # Check if rg is a new ring:
                    ifNew = True
                    for rg_k in rglist:
                        if rg_curr == rg_k:
                            ifNew = False
                            break
                    if ifNew:
                        rglist.append( rg_curr )
                        # print( rglist )

            continue
        else:
            # Ring growth:
            rg_curr.append( at_next )

        rings_from_an_atom( rglist, rg_curr, nblist, at_next )

    return rglist
# enddef rings_from_an_atom()


# Get ring list of a molecular graph
#
#  nblist: List, each entry indicating indices (starting from 0) of all 
#          neighbors of each atom
#  at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of atoms in each ring
def rings( nblist, at=[] ):
    NAt = len( nblist )

    rglist = []
    for iat in range( NAt ):
        rglist = rings_from_an_atom( rglist, [ iat ], nblist, iat )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return rglist

    rglist0 = []
    for rg in rglist:
        x = [ at[a] for a in rg ]
        rglist0.append( x )
    return rglist0
# enddef rings()


# Determine if two rings are adjacent
#
def ifAdjacetRings( rg1, rg2 ):
    nrg1 = len( rg1 )
    nrg2 = len( rg2 )
    for i in range( nrg1 ):
        if i < nrg1-1:
            bnd1 = [ rg1[i], rg1[i+1] ]
        else:
            bnd1 = [ rg1[i], rg1[0] ]
        # Reorder the two atom labels in bnd1:
        if bnd1[1] < bnd1[0]:
            bnd1 = [ bnd1[1], bnd1[0] ]

        for j in range( nrg2 ):
            if j < nrg2-1:
                bnd2 = [ rg2[j], rg2[j+1] ]
            else:
                bnd2 = [ rg2[j], rg2[0] ]
            # Reorder the two atom labels in bnd2:
            if bnd2[1] < bnd2[0]:
                bnd2 = [ bnd2[1], bnd2[0] ]

            if bnd1 == bnd2:
                return True

    return False
# enddef ifAdjacetRings()


# Get neighbors of each ring
#
#   rglist: List of indices of atoms (starting from 0) of each ring
#
# Return:
#   List, each entry indicating indices (starting from 0) of all neighbors of 
#   each ring
def ringNeighbors( rglist ):
    N = len( rglist )

    # Get connectivity between rings:
    A = np.zeros( ( N, N ), dtype='i' )
    for i in range( N ):
        for j in range( i+1, N ):
            A[ i, j ] = ifAdjacetRings( rglist[i], rglist[j] )
            A[ j, i ] = A[ i, j ]
    #print( A )

    nblist = []
    for i in range( N ):
        ix = [ j for j, a in enumerate(A[i,:]) if a > 0 ] 
        nblist.append( ix )
    return nblist
# enddef RingNeighbors()


# Read rings from *.rings file and return them to RG[]
#   ringsFileName: Name of *.rings file where atomic indices of rings are stored
def readRings( ringsFileName ):
    with open( ringsFileName ) as f:
        RG = [ [int(x) for x in line.split()] for line in f ]

    return RG
# enddef readRings()


def getRings_from_coord( xyz ):
    # Adjacent matrix:
    A = adjMat( xyz )
    # List of neighbors:
    nblist = neighbors( A )
    # Ring list:
    return rings( nblist )

# A corrected version of getRings_from_coord() for special cases like
# kekulene[1, 1, 2, 1, 1, 2]:
def getRings_from_coord_corrected( xyz ):
    RG0 = getRings_from_coord( xyz )
    lm = adjMat( xyz )
    CN = np.sum( lm, axis=0 ) # Coordination number of C atoms

    RG = []
    ix_ringdel0 = []
    for j in range( len(RG0) ):
        rg = RG0[j]
        ifProRing = True
        for k in range( len(rg) ):
            if CN[ rg[k] ] == 2:
                ifProRing = False
                break
            elif CN[ rg[k] ] != 3:
                print( 'ERROR: impossible coordination number (%i) for C#%i' %\
                        CN[ rg[k] ], rg[k] )
        if not ifProRing:
            RG.append( rg )
        else:
            ix_ringdel0.append( j )

    # Recover misremoved, valid ring:
    lm_ring0 = ringLinkage( RG0 )
    lm_ring = ringLinkage( RG )
    CN_ring = np.sum( lm_ring, axis=0 )
    ix_lowcn = np.where( CN_ring < 2 )[0]
    for j in range( len(ix_lowcn) ):
        for k in range( len(RG0) ):
            if np.array_equal( np.sort(RG[ ix_lowcn[j] ]), np.sort(RG0[k]) ):
                ix_lowcn0[j] = k

    for j in range( len(ix_ringdel0) ):
        conn_num = 0
        for k in range( len(ix_lowcn) ):
            if lm_ring0[ ix_ringdel0[j], ix_lowcn0[k] ]:
                conn_num += 1
        if conn_num == 2:
            RG.append( RG0[ ix_ringdel0[j] ] )
        elif conn_num > 2:
            print( 'ERROR: impossible connectivity number of ring: %i' % \
                    RG0[ ix_ringdel0[j] ] )
#    if xyz[0,0] > 1:
#        print( RG )
#        input()

    return RG


def ringLinkage( RG ):
    NR = len( RG )
    lm_ring = np.zeros( (NR, NR), dtype='i' )
    for j in range( NR ):
        rg1 = RG[j]
        for k in range( j+1, NR ):
            rg2 = RG[k]
            if np.intersect1d( rg1, rg2 ).size == 2:
                lm_ring[j,k] = 1
                lm_ring[k,j] = 1

    return lm_ring


def rings_kinfi( lm, NAt1, RG1, RG2, ixOrder1, ixOrder2, h6, k1, d ):
    # Cutting position to which all left rings (including ix_cut) on the 
    # contact side have z = -zm,  and all rings right to ix_cut on the contact 
    # side have z = zm (for kekulene [h1,...,h6])
    if d <= h6 - k1:
        ix_cut1 = d + k1 # Cut at east(E) edge
        ix_cut2 = k1 # Cut at southeast(SE) edge
        #print( 'Case I (d = %i; d <= h6 - k1)' % d )
    else:
        ix_cut1 = h6 # Cut at northeast(NE) edge
        ix_cut2 = h6 - d # Cut at east(E) edge
        #print( 'Case II (d = %i; d > h6 - k1)' % d )

    if ix_cut1 < 0:
        print( 'ERROR: ix_cut1 (%i) must be a nonnegative integer' % ix_cut1 )
        exit(1)
    if ix_cut2 < 0:
        print( 'ERROR: ix_cut2 (%i) must be a nonnegative integer' % ix_cut2 )
        exit(1)

    NR1 = len( RG1 )
    NR2 = len( RG2 )
    for i in range( NR2 ):
        for k in range( len(RG2[i]) ):
            RG2[i][k] += NAt1

    ix_cutRing1 = ixOrder1[ix_cut1]
    ix_cutnextRing1 = ixOrder1[ix_cut1+1]
    cutBond1 = np.intersect1d( RG1[ix_cutRing1], RG1[ix_cutnextRing1] )
    ix_cutRing2 = ixOrder2[ix_cut2]
    ix_cutnextRing2 = ixOrder2[ix_cut2+1]
    cutBond2 = np.intersect1d( RG2[ix_cutRing2], RG2[ix_cutnextRing2] )

    RG = []
    # For rings in top kekulene 1:
    for j in range( NR1 ):
        if j == ix_cut1+1:
            rg = RG1[ ixOrder1[j] ]
            ix = np.where( (rg == cutBond1[0]) | (rg == cutBond1[1]) )[0]
            at_remain = np.setdiff1d( rg, cutBond1 )
            CN_at_remain = np.zeros( len(at_remain), dtype='i' )
            k = 0
            for at1 in at_remain:
                for at2 in np.setdiff1d( at_remain, at1 ):
                    if lm[ at1, at2 ]:
                        CN_at_remain[k] += 1
                k += 1
            at_link = at_remain[ np.where( CN_at_remain == 1 ) ]
            if len( at_link ) != 2:
                print( 'ERROR: Something wrong with linking atoms' )
                exit(1)
            # Neighbors of at_link[0]:
            nbs_1 = np.where( lm[ at_link[0], : ] == 1 )[0]
            nb_1 = nbs_1[ nbs_1 > NAt1-1 ]
            # Neighbors of at_link[1]:
            nbs_2 = np.where( lm[ at_link[1], : ] == 1 )[0]
            nb_2 = nbs_2[ nbs_2 > NAt1-1 ]
            RG.append( np.append( at_remain, [nb_1, nb_2] ).tolist() )
        else:
            RG.append( RG1[ ixOrder1[j] ] )

    # For rings in top kekulene 2:
    for j in range( NR2 ):
        if j == ix_cut2+1:
            rg = RG2[ ixOrder2[j] ]
            ix = np.where( (rg == cutBond2[0]) | (rg == cutBond2[1]) )[0]
            at_remain = np.setdiff1d( rg, cutBond2 )
            CN_at_remain = np.zeros( len(at_remain), dtype='i' )
            k = 0
            for at1 in at_remain:
                for at2 in np.setdiff1d( at_remain, at1 ):
                    if lm[ at1, at2 ]:
                        CN_at_remain[k] += 1
                k += 1
            at_link = at_remain[ np.where( CN_at_remain == 1 ) ]
            if len( at_link ) != 2:
                print( 'ERROR: Something wrong with linking atoms' )
                exit(1)
            # Neighbors of at_link[0]:
            nbs_1 = np.where( lm[ at_link[0], : ] == 1 )[0]
            nb_1 = nbs_1[ nbs_1 < NAt1 ]
            # Neighbors of at_link[1]:
            nbs_2 = np.where( lm[ at_link[1], : ] == 1 )[0]
            nb_2 = nbs_2[ nbs_2 < NAt1 ]
            RG.append( np.append( at_remain, [nb_1, nb_2] ).tolist() )
        else:
            RG.append( RG2[ ixOrder2[j] ] )

    return RG

