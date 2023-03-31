# addH: Functions to automatically add hydrogen atoms to carbon backbones
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
from topo import *

def addH( xyz ):
    MAX_ITER = 100
    d_CH = 1.09 # in Angstrom; C--H bond length
    MIN_HH_DIST = 1.5 # in Angstrom

    NC = xyz.shape[0] # Number of C atoms
    lm = adjMat( xyz ) # Adjacent matrix
    CN = np.sum( lm, axis=0 ) # Coordination number of C atoms
    ixC2 = np.where( CN == 2)[0] # All 2-coordinated C atoms
    NH = len( ixC2 ) # Number of H atoms
    # assert( NC == NH*2 )

    # Determine coordinates of H atoms
    xyzH = np.zeros( (NH, 3) )
    for j in range( NH ):
        iC = ixC2[j]
        iC_nb = np.where( lm[iC,:] == 1 )[0]
        v1 = xyz[ iC_nb[0], : ] - xyz[ iC, : ]
        v1 /= np.linalg.norm( v1 )
        v2 = xyz[ iC_nb[1], : ] - xyz[ iC, : ]
        v2 /= np.linalg.norm( v2 )
        vH = -( v1 + v2 )
        vH /= np.linalg.norm( vH )
        xyzH[ j, : ] = xyz[ iC, : ] + vH*d_CH


    # Adjust H atoms that are too close to each other:
    for j in range( NH ):
        for k in range( j+1, NH ):
            v = xyzH[ k, : ] - xyzH[ j, : ]
            r = np.linalg.norm( v )
            if r < MIN_HH_DIST:
                print( 'H atoms %i--%i too close to each other: %.3f Angstrom' \
                        % (NC+j+1, NC+k+1, r) )
            
                iter = 1
                while iter <= MAX_ITER:
                    # Adjust H positions according to H--H distance:
                    u = v / r
                    xyzH[ j, : ] = xyzH[ j, : ] - u*0.2
                    xyzH[ k, : ] = xyzH[ k, : ] + u*0.2
                
                    # Adjust H-j and H-k positions according to C--H distance:
                    for iH in [ j, k ]:
                        for iC in range( NC ):
                            w = xyzH[ iH, : ] - xyz[ iC, : ]
                            rw = np.linalg.norm( w )
                            if rw <= d_CH*1.2:
                                break
                        u = w / rw
                        xyzH[ iH, : ] = xyzH[ iH, : ] + u*( d_CH - rw )
                
                    # Update H--H distance:
                    r = np.linalg.norm( xyzH[ k, : ] - xyzH[ j, : ] )
                    if r >= MIN_HH_DIST:
                        print( 'H%i--H%i distance converged to %.3f Angstrom '
                              'after %i iterations' % (NC+j, NC+k, r, iter) )
                        break
                    iter += 1
                    if iter == MAX_ITER:
                        print( 'WARNING: H%i--H%i distance NOT converged: ' % \
                            (NC+j, NC+k) )
                        print( '%.3f Angstrom', r )

    xyz = np.row_stack( (xyz, xyzH) )

    elem = []
    for i in range( NC ):
        elem.append( 'C' )
    for i in range( NC, NC+NH ):
        elem.append( 'H' )

    return (xyz, elem)


# A corrected version of addH() for cases with fused rings 
# (e.g., K-infi [1,1,1,1,1,1|1,1,4,1,1,4](-1))
def addH_corr( xyz, RG_corr ):
    MAX_ITER = 100
    d_CH = 1.09 # in Angstrom; C--H bond length
    MIN_HH_DIST = 1.5 # in Angstrom

    NC = xyz.shape[0] # Number of C atoms
    lm = adjMat( xyz ) # Adjacent matrix
    CN = sum( lm ) # Coordination number of C atoms
    ixC2 = np.where( CN == 2)[0] # All 2-coordinated C atoms
    NH = len( ixC2 ) # Number of H atoms
    # assert( NC == NH*2 )

    #======================================================================
    # Bug fix for fused-ring k-infinitenes:
    RG = getRings_from_coord( xyz )
    #RG = getRings_from_coord_corrected( xyz )
    lm_ring = ringLinkage( RG )
    CN_ring = sum( lm_ring )
    max_CN_ring = max( CN_ring )
    if max_CN_ring > 2:
        ix_proRing = []
        for j in range( len( RG ) ):
            flag = False
            for k in range( len( RG_corr ) ):
                if np.array_equal( np.sort(RG[j]), np.sort(RG_corr[k]) ):
                    flag = True
                    break
            if not flag:
                ix_proRing.append(j)
        #print( ix_proRing ); input()

        # Find all neighboring rings of the problematic rings:
        ix_nbRing = []
        for j in range( len( ix_proRing ) ):
            ix_rgnb = np.where( lm_ring[ ix_proRing[j], : ] == 1 )[0]
            ix_nbRing += ix_rgnb.tolist()
        ix_nbRing = np.unique( ix_nbRing, axis=0 ) # uniq
        # Exclude proRing from nbRing:
        ix_nbRing = np.setdiff1d( ix_nbRing, ix_proRing ) 
        #for k in range( len( ix_nbRing ) ):
        #    print( RG[ ix_nbRing[k] ] )
        #input()

        list_ix_proC = [] # List of indices of problematic C atoms
        for j in range( len( ix_proRing ) ):
            proRing = RG[ ix_proRing[j] ]
            # Find problematic C pairs that should not be bonded:
            for iat1 in range( len( proRing ) ):                
                if iat1 < len( proRing )-1:
                    iat2 = iat1 + 1
                else:
                    iat2 = 0
                at1 = proRing[ iat1 ]
                at2 = proRing[ iat2 ]
                #print(at1,at2); input()
                ifFoundBond = False
                for k in range( len( ix_nbRing ) ):
                    nbRing = RG[ ix_nbRing[k] ]
                    for jat1 in range( len( nbRing ) ):
                        if jat1 < len( nbRing )-1:
                            jat2 = jat1 + 1
                        else:
                            jat2 = 0
                        nb_at1 = nbRing[ jat1 ]
                        nb_at2 = nbRing[ jat2 ]
                        #print(nb_at1,nb_at2); input()
                        if (at1==nb_at1 and at2==nb_at2) or \
                                (at1==nb_at2 and at2==nb_at1):
                            ifFoundBond = True
                            break
                    if ifFoundBond:
                        break
                if not ifFoundBond:
                    list_ix_proC.append( [ at1, at2 ] )
        #print( list_ix_proC ), input()

        for j in range( len(list_ix_proC) ):
            ix_proC = list_ix_proC[ j ]
            if len( ix_proC ) != 2:
                print( 'ERROR: something wrong with ix_proC = %i' % ix_proC )
                print( 'ABORTED!!!' )

            # Enlarge the distance between the two atoms in ix_proC:
            Rmax = 2.00 + 0.01 # Target distance
            ix1 = ix_proC[0]
            ix2 = ix_proC[1]
            print( 'Enlarging C%i--C%i distance ...' % (ix1, ix2) )
            v12 = xyz[ ix2, : ] - xyz[ ix1, : ]
            dv = ( Rmax - np.linalg.norm(v12) )/2 * v12 / np.linalg.norm(v12)
            xyz[ ix1, : ] = xyz[ ix1, : ] - dv
            xyz[ ix2, : ] = xyz[ ix2, : ] + dv

            # All C atoms where H atoms are added:
            ixC2 = np.append( ixC2, ix_proC ) 


        ixC2 = np.unique( ixC2, axis=0 )

        lm = adjMat( xyz ) # UPDATE adjacent matrix
        CN = sum( lm ) # UPDATE coordination number of C atoms
        ixC2_new = np.where( CN == 2 )[0] # UPDATE all 2-coordinated C atoms
        #print( np.sort(ixC2) ), print( np.sort(ixC2_new) ), input()
        assert( np.array_equal( np.sort(ixC2), np.sort(ixC2_new) ) )     
        NH = len( ixC2 ) # UPDATE the number of H atoms
    #======================================================================


    # Determine coordinates of H atoms
    xyzH = np.zeros( (NH, 3) )
    for j in range( NH ):
        iC = ixC2[j]
        iC_nb = np.where( lm[iC,:] == 1 )[0]
        v1 = xyz[ iC_nb[0], : ] - xyz[ iC, : ]
        v1 /= np.linalg.norm( v1 )
        v2 = xyz[ iC_nb[1], : ] - xyz[ iC, : ]
        v2 /= np.linalg.norm( v2 )
        vH = -( v1 + v2 )
        vH /= np.linalg.norm( vH )
        xyzH[ j, : ] = xyz[ iC, : ] + vH*d_CH


    # Adjust H atoms that are too close to each other:
    for j in range( NH ):
        for k in range( j+1, NH ):
            v = xyzH[ k, : ] - xyzH[ j, : ]
            r = np.linalg.norm( v )
            if r < MIN_HH_DIST:
                print( 'H atoms %i--%i too close to each other: %.3f Angstrom' \
                        % (NC+j+1, NC+k+1, r) )
            
                iter = 1
                while iter <= MAX_ITER:
                    # Adjust H positions according to H--H distance:
                    u = v / r
                    xyzH[ j, : ] = xyzH[ j, : ] - u*0.1
                    xyzH[ k, : ] = xyzH[ k, : ] + u*0.1
                
                    # Adjust H-j and H-k positions according to C--H distance:
                    for iH in [ j, k ]:
                        for iC in range( NC ):
                            w = xyzH[ iH, : ] - xyz[ iC, : ]
                            rw = np.linalg.norm( w )
                            if rw <= d_CH*1.2:
                                break
                        u = w / rw
                        xyzH[ iH, : ] = xyzH[ iH, : ] + u*( d_CH - rw )
                
                    # Update H--H distance:
                    r = np.linalg.norm( xyzH[ k, : ] - xyzH[ j, : ] )
                    if r >= MIN_HH_DIST:
                        print( 'H%i--H%i distance converged to %.3f Angstrom '
                              'after %i iterations' % (NC+j, NC+k, r, iter) )
                        break
                    iter += 1
                    if iter == MAX_ITER:
                        print( 'WARNING: H%i--H%i distance NOT converged: ' % \
                            (NC+j, NC+k) )
                        print( '%.3f Angstrom', r )

    xyz = np.row_stack( (xyz, xyzH) )

    elem = []
    for i in range( NC ):
        elem.append( 'C' )
    for i in range( NC, NC+NH ):
        elem.append( 'H' )


    return (xyz, elem)

