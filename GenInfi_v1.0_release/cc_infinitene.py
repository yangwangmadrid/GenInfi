# cc_infinitene: Create coordinates of generalized coronene-clarene infinitene
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
# Created on July 14, 2022
# By Yang WANG (yangwang@yzu.edu.cn)
#
#
########################################################################
# MORE DETAILS:
#
# Create coordinates of generalized coronene-clarene infinitene:
# i.e., CC-infinitene <coronene|k1,k2,k3,k4,k5,k6>(d)
# subject to the following conditions:
# i) k1,...,k6 >=2 and are all EVEN numbers;
# iii)   -k1 <= d <= 0 and d must be a EVEN number
# NOTE: One should take into account the chirality of the clarene:
#    (1) Let define the canonical notation of a clarene <k1,...,k6>:
#        i.e., (k1,...,k6) is the smallest sequence
#          We can compare two sequences, (k1,...,k6) and (k1',...,k6'):
#            If k1<=k1' AND k2<=k2', ..., AND k6<=k6'
#            Then (k1,...,k6) <= (k1',...,k6')
#    (2) We can define the chirality of a clarene <k1,...,k6>:
#        If the canonical notation of (k1,...,k6) is EQUAL TO the canonical
#        notation of its reverse fliplr(k1,...,k6), then <k1,...,k6> is
#        ACHIRAL (A); 
#        Otherwise, <k1,...,k6> is CHIRAL:
#          If (k1,...,k6) < fliplr(k1,...,k6),
#          Then <k1,...,k6> is RIGHT (R);
#          If (k1,...,k6) > fliplr(k1,...,k6),
#          Then <k1,...,k6> is LEFT (L).
#     (3) Since coronene has an A chirality, when <k1,...,k6> is chiral, 
#         we may have 2 different CC-infinitenes:
#             (A,R), (A,L)
#         When <k1,...,k6> is achiral, we only have a unique CC-infinitene:
#             (A,A)

from c_infinitene import *

# For construction of coronene-clarene infinitenes:
def zcoord_ccinfi_coronene( xyz_clr1, RG, hexa_cc, ixOrder, zm, k1, d ):
    # Cutting position to which all left rings (including ix_cut) on the 
    # contact side have z = -zm,  and all rings right to ix_cut on the 
    # contact side have z = zm (for clarene <h1,...,h6>)
    if d == -k1: # special case
        ix_cut = 0
    else:
        ix_cut = 1

    NAt = xyz_clr1.shape[0]

    # Ring indices of the bottom side:
    ix_bottom = ixOrder[ 0 ]

    # Determine z-coordinates of hexagon centers:
    NR = hexa_cc.shape[0]
    z_hexa_cc = np.zeros( (NR, 1) )

    # 1) For bottom side:
    if ix_cut == 0:
        Nhexa_bottom = 2
    elif ix_cut == 1:
        Nhexa_bottom = 3
    else:
        print( 'ERROR: ix_cut (%i) must be either 0 or 1' % ix_cut )

    Nhexa_remain = NR - Nhexa_bottom
    zstep = ( 2 * zm )/( Nhexa_remain + 1 )
    z_hexa_cc[ range(ix_cut+1) ] = -abs( zm )
    z_hexa_cc[ range(ix_cut+1,Nhexa_bottom) ] = abs( zm )

    # 2) All other sides:
    zval = np.zeros( (NR-Nhexa_bottom, 1) )
    for i in range( NR-Nhexa_bottom ):
        zval[i] =  zm - zstep + i*( -zstep )
        z_hexa_cc[ range(Nhexa_bottom,NR) ] = zval

    # Determine z-coordinates of ATOMS:
    z_atoms = 999999 * np.ones( (NAt, ) )
    for j in range( NR ):
        rg = RG[ ixOrder[j] ]
        for k in range( len(rg) ):
            iat = rg[k]
            if z_atoms[iat] == 999999:
                z_atoms[iat] = z_hexa_cc[j]

            else:
                # Exclude the special case of kek2 on the top side when it 
                # contains an odd number of hexagons on the top side:
                if abs( z_atoms[iat] - abs(zm) ) < 1E-3 and \
                        abs( z_hexa_cc[j] + abs(zm) ) < 1E-3:
                        pass
                elif abs( z_atoms[iat] + abs(zm) ) < 1E-3 and \
                    abs( z_hexa_cc[j] - abs(zm) ) < 1E-3:
                        pass
                elif abs( z_atoms[iat] - abs(zm) ) < 1E-3 and \
                    abs( z_hexa_cc[j] + abs(zm) - abs(zstep) ) < 1E-3:
                        pass
                elif abs( z_atoms[iat] + abs(zm) ) < 1E-3 and \
                    abs( z_hexa_cc[j] - abs(zm) + abs(zstep) ) < 1E-3:
                        pass
                else:
                    z_atoms[iat] = ( z_hexa_cc[j] + z_atoms[iat] ) / 2
    

    # Update xyz_clr1:
    xyz_clr1[ :, 2 ] = z_atoms

    return xyz_clr1


def zcoord_ccinfi_bottom( xyz_clr2, RG, hexa_cc, ixOrder, zm, k1, d ):
    # Cutting position to which all left rings (including ix_cut) on the 
    # contact side have z = -zm,  and all rings right to ix_cut on the 
    # contact side have z = zm (for clarene <k1,...,k6>)
    if d == -k1: # special case
        ix_cut = k1
        print( 'Case I (d = %i; d == -k1)' % d )
    else: # special cases
        ix_cut = 1 - d
        print( 'Case II (d = %i; -k1 < d <= 0)' % d )

    NAt = xyz_clr2.shape[0]

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
    z_hexa_cc[ ix_rest ] = zval

    # Determine z-coordinates of ATOMS:
    z_atoms = 999999 * np.ones( (NAt, ) )
    for j in range( NR ):
        rg = RG[ ixOrder[j] ]
        for k in range( len(rg) ):
            iat = rg[k]
            if z_atoms[iat] == 999999:
                z_atoms[iat] = z_hexa_cc[j]

    # Update xyz_clr2:
    xyz_clr2[ :, 2 ] = z_atoms

    return xyz_clr2


# Create coordinates of coronene:
def coronene():
    pi = np.pi
    cos = np.cos
    sin = np.sin

    #========== Constants ==========
    a = 1.42; # (in Angstrom) C--C bond length in benzene
    s = np.sqrt(3)/2

    # I. Create the central points of all hexagons:
    xyz = np.zeros( (24, 3) ) # 24 C atoms in coronene

    # Only work out 2D coordinates, (x,y), for the molecule is planar (z = 0)
    # 1. Get all centers of the 6 corner hexagons:
    t = pi/6 + np.arange(6)*pi/3
    # 2D coord. of hexagon centers
    xy_h1 = 2*s*a * np.array( [ cos(t) ] ).reshape((6,1))
    xy_h2 = 2*s*a * np.array( [ sin(t) ] ).reshape((6,1))
    xy_h = np.column_stack( (xy_h1,xy_h2) )

    # II. Create the corner points of each hexagon:
    # 1. Generate the first set of six radial vectors centered at the central
    # point of each hexagon, which is already generated above:
    t = np.arange(6) * 2*pi / 6
    v_n1 = a * np.array( [ cos(t) ] ).reshape((6,1))
    v_n2 = a * np.array( [ sin(t) ] ).reshape((6,1))
    v_n = np.column_stack( (v_n1,v_n2) )

    iC = 0
    for ih in range(6):
        # Translation by radial vectors:
        for j in range(6):
            xy_tmp = xy_h[ ih, : ] + v_n[ j, : ]
            # Rule out already generated C atoms:
            nC = iC
            flag = True
            for k in range( nC+1 ):
                if np.linalg.norm( xy_tmp - xyz[ k, 0:2 ] ) < 1E-3:
                    flag = False # Found duplicated
                    break
            if flag: # Not duplicated
                xyz[ iC, 0:2 ] = xy_tmp
                iC += 1

    # Rings:
    RG1 = getRings_from_coord( xyz )
    # Reorder RG1[], because when forming infinitene, 
    # the algorithm requires that coronene be placed on top of the clarene:
    RG = [ [] ]*6
    for j in range(7):
        vr = xyz[ RG1[j], : ].mean( axis=0 )
        # Skip the central ring:
        if np.linalg.norm( vr ) < 1E-6:
            continue
        # Get radial angle:
        _, theta = cart2pol( vr[0], vr[1] )
        theta *= 180/pi
        irg = round( (theta+30)/60 )+1
        if irg == -1:
            irg = 5
        RG[ irg ] = RG1[j]
    ixOrder = range(6)

    # Hexagon centers:
    hexa_cc = np.zeros( (6, 3) )
    # Sextet centers:
    sex_cc = np.zeros( (3, 3) )
    for j in range(6):
        hexa_cc[ j, : ] = xyz[ RG[j], : ].mean( axis=0 )
        if j % 2 == 0:
            sex_cc[ j//2, : ] = hexa_cc[ j, : ]

    return (xyz, hexa_cc, sex_cc, RG, ixOrder)



def cc_infinitene( kk, d, ifWrite=False ):
#  kk: an array containing k1, ..., k6 for clarene <k1,k2,k3,k4,k5,k6>
#  If ifWrite is TRUE, then write coordinate to an external *.xyz file
    #====================== Check validity of inputs ======================
    if d % 2 == 1:
        print( 'ERROR: d must be an even number (d = %i)' % d )
    # Range of d:
    if d < -kk[0] or d > 0:
        print( 'ERROR: d must lie in range %i <= d <= %i (d = %i)' % \
                (-kk[0], 0, d) )
        exit(1)
    #======================================================================

    #============================= Constants ==============================
    a = 1.42 # (in Angstrom) C--C bond length in benzene
    zm = 1.60 # A fixed distance between upper clarene <h1,...h6>
              # and upper clarene <k1,...,k6>
    #======================================================================

    # Initialize:
    Nhexa = 6 + sum( kk ) # Coronene has 6 rings
    print( 'CC-infinitene <coronene|%i,%i,%i,%i,%i,%i>(%i):' % \
            (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
    print( '%i hexagons' % Nhexa )
    NAt = 4 * Nhexa
    print( '%i atoms' % NAt )
    xyz = np.zeros( (NAt, 3) )

    # I. Generate coordinates of coronene and the clarene <k1,...,k6>:
    # (1) Coronene:
    # Coord. and all hexagons, sextets and their centers:
    k1 = kk[0]
    [ xyz_clr1, hexa_cc1, sex_cc1, RG1, ixOrder1 ] = coronene()
    NAt1 = xyz_clr1.shape[0]
    #print( xyz_clr1 )

    # (2) Clarene: No need to reorder kk[], 
    #              as k1 is already the contact (TOP) side:
    xyz_clr2 = clarene( kk, False, False )
    # Determine all hexagons, sextets and their centers:
    RG2 = getRings_from_coord_corrected( xyz_clr2 )
    # Check validity of clarene kk:
    if len( RG2 ) == 0:
        print( '!!! WARNING: wrong structure of clarene <%i,%i,%i,%i,%i,%i>' \
                % kk )
        print( '!!! Aborted' )
        xyz = []
        return
    [ hexa_cc2, sex_cc2, ixOrder2 ] = \
            hexagon_centers( xyz_clr2, RG2, 'bottom' )

    # II. Generate (x,y)-coord. of joined infinitene <coronene|k1,...,k6>(d):
    # Facts:
    #   For coronene, all three Clar rings form a supertriangle,
    #   and the actual length of side h_i of the supertriangle is 3*a
    #
    # 1. Translate coronene so that the bottom sextet ring coincides with the 
    #    x-axis and meanwhile is also placed at the origin:
    xyz[ range(NAt1), : ] = moveCoord( xyz_clr1, -sex_cc1[0,:] )

    # 2. Translate clarene <k1,...,k6> so that the k1-side coincides with the
    # x-axis meanwhile the leftmost sextet ring of h6-side is at the origin: 
    xyz[ range(NAt1,NAt), : ] = moveCoord( xyz_clr2, -sex_cc2[0,:] )
    sex_cc2 = moveCoord( sex_cc2, -sex_cc2[0,:] )

    # 3. Translate clarene <k1,...,k6> along x-direction by the variable d:
    shift_x = d * 3/2 *a
    xyz[ range(NAt1,NAt), : ] = \
            moveCoord( xyz[ range(NAt1,NAt), : ], [ shift_x, 0., 0. ] )

    # III. Generate z-coordinates of the joined infinitene:
    #          <coronene|k1,k2,k3,k4,k5,k6>(d):
    # 1. Coronene:
    xyz_clr1 = xyz[ range(NAt1), : ]
    xyz_clr1 = zcoord_ccinfi_coronene( xyz_clr1, RG1, hexa_cc1, ixOrder1, \
            zm, k1, d )
    #print( xyz_clr1 ); input()
    xyz[ range(NAt1), : ] = xyz_clr1

    # 2. Clarene <k1,k2,k3,k4,k5,k6>:
    xyz_clr2 = xyz[ range(NAt1,NAt), : ]
    xyz_clr2 = zcoord_ccinfi_bottom( xyz_clr2, RG2, hexa_cc2, ixOrder2, \
            zm, k1, d )
    #print( xyz_clr2 ); input()
    xyz[ range(NAt1,NAt), : ] = xyz_clr2

   # HMO total energy and eigenvalues:
    hmoSol = hmo( xyz ) # NOTE: xyz[] is so far only for CARBON atoms
    print( '\nHuckel pi energy = %.6f |beta|' % hmoSol.Etot )

    # Add hydrogen atoms:
    [ xyz, elem ] = addH( xyz )

    # Check if number of H atoms is correct:
    NC = NAt
    NH = xyz.shape[0] - NAt
    print( '%i H atoms in CC-infinitene '\
            '<coronene|%i,%i,%i,%i,%i,%i>(%i)' % \
            (NH, kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
    # Exceptional case::
    if NH != NC // 2:
        print( '!!! WARNING: Wrong number of H atoms found in CC-infinitene '\
                '<coronene|%i,%i,%i,%i,%i,%i>(%i)' % \
                (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
        print( '!!! %i H atoms are expected' % (NC // 2) )
        print( 'This infinitene structure should be adjusted manually' )

    # Output xyz to external *.gjf file:
    if ifWrite:
        fmt = 'xyz'
        outp = 'ccinf_R%i-1_1_1_1_1_1-%i_%i_%i_%i_%i_%i_d%i.%s' % \
                (Nhexa, kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d, fmt)
        write_xyz( outp, xyz, elem, outp )

    return (xyz, elem, hmoSol.E)


# Enumerate coronene-clarene infinitenes by assigning cyclic sequence of 
# the clarene, kk, with all possible values of shift d.
#
# NOTE: For the input sequence, kk, the CONTACT side is FIXED,
#       i.e., coronene|kk(1) is fixed, while other indices kk(j) can flip 
#       CYCLICALLY. Thus, the script will automatically enumerate all 
#       possible chiralities of kk.
#       Hence, this script is different from (more specific than) the
#       script enum_k_all_gc_infinitene().
#
def enum_k_cc_infinitene( kk ):
    print( '='*60 )

    # I. Canonnical notation of CC-infinitene:
    # <coronene|k1,k2,k3,k4,k5,k6>(d)
    # subject to the following conditions:
    #  i) k1,...,k6 >=2 and are all EVEN numbers;
    # ii) -k1 <= d <= 0 and d must be a EVEN number
    k_sum = sum( kk )

    # II. Generate <coronene|kk>:
    NMol = 0
    xyz_all = []
    HMOeig_all = []
    kd_all = []
    # Generate the CC-infinitene:
    for d in range( -kk[0], 1, 2 ):
        print( '-'*50 )
        print( 'Generating', end='' )
        print( ' <coronene|%i,%i,%i,%i,%i,%i>(%i) ...' % \
                        (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
        NMol += 1
        [ xyz, elem, HMOeig ] = cc_infinitene( kk, d, False )
        xyz_all.append( xyz )
        HMOeig_all.append( HMOeig )
        kd_all.append( [ kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d ] )
        print( '-'*50 )

    print( '\n%i original structures generated' % NMol )
    HMOeig_all = np.array( np.round( HMOeig_all, 8 ) )
    HMOeig_all, ix = np.unique( HMOeig_all, axis=0, return_index=True )
    NMol_uniq = HMOeig_all.shape[0]
    print( '%i unique structures obtained' % NMol_uniq )
    # Update kd_all and xyz_all:
    kd_all = np.array( kd_all, dtype='i' )[ ix, : ]
    xyz_uniq = []
    for i in range( NMol_uniq ):
        xyz_uniq.append( xyz_all[ ix[i] ] )
    print( '='*50 )

    return (kd_all, xyz_uniq, elem, HMOeig_all)


# Enumerate coronene-clarene infinitenes by assigning cyclic sequence the 
# clarene, kk, with all possible values of shift d.
#
# NOTE: The input sequences kk are cyclically invariant and this script 
# will automatically enumerate all nonequivalent <k1,k2,k3,k4,k5,k6> 
# combinations.
# Hence, this script is more GENERAL than enum_k_cc_infinitene().
#
def enum_k_all_cc_infinitene( kk ):
    kk_uniq = all_uniq_cyclicSeq( kk )

    kd_all = []
    xyz_all = []
    HMOeig_all = []
    NMol = 0
    for i2 in range( kk_uniq.shape[0] ):
        kk_ = kk_uniq[ i2, : ]
#        print(kk_)
#        input()
        kd_arr, xyz_arr, elem, HMOeig_arr = enum_k_cc_infinitene( kk_ )
        for j in range( len( xyz_arr ) ):
            NMol += 1
            kd_all.append( kd_arr[j,:] )
            xyz_all.append( xyz_arr[j] )

        # print( HMOeig_arr.shape )
        # np.append( HMOeig_all, [HMOeig_arr], axis=0 )
        if i2 == 0:
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
        kd = kd_all[ ix[i] ]
        xyz = xyz_all[ ix[i] ]
        NR = len(xyz) // 6
        Nocc = 2*NR
        Ehmo = 2*np.sum( HMOeig_all[ i, range(Nocc) ], axis=0 )
        fmt = 'xyz'
        outp = 'ccinf_R%i-1_1_1_1_1_1-%i_%i_%i_%i_%i_%i_d%i.%s' % \
                (NR,kd[0],kd[1],kd[2],kd[3],kd[4],kd[5], kd[6], fmt)
        title = 'CC-infinitene <coronene|%i,%i,%i,%i,%i,%i>(%i) '\
                'EHMO: %.8f' % (kd[0],kd[1],kd[2],kd[3],kd[4],kd[5], \
                kd[6], Ehmo)

        write_xyz( outp, xyz, elem, title )

    return (NMol, NMol_uniq)


# Enumerate all possible coronene-clarene infinitenes of a given ring size,
# Nring, taking into account of all possible shifts, d.
#
def enum_Nring_all_cc_infinitene( Nring ):

    L2 = Nring - 6 # Ring sizes of the clarene

    # Enumerate all possible equiangular_hexagon_even for L2:
    kk_arr = equiangular_hexagon_even( L2, False )

    counter2 = 0
    for i2 in range( kk_arr.shape[0] ):
        kk = kk_arr[ i2, : ]
        counter2 += 1
        print( '#### (%i):' % counter2 )
        print( '####    coronene | clarene 2: <%2i,%2i,%2i,%2i,%2i,%2i>' % \
                (kk[0], kk[1], kk[2], kk[3], kk[4], kk[5]) )
        NMol, NMol_uniq = enum_k_all_cc_infinitene( kk )
        #counter_all += NMol
        #counter_all_uniq += NMol_uniq
    print( '#'*60 )

    print( '\nAll [%i]CC-infinitene structures have been successfully'
          ' generated' % Nring )
    Niso = len( glob.glob( 'ccinf_R%i-*_*-*_*_d*.xyz' % Nring ) )
    print( 'A total of %i unique isomers (excluding enantiomers)' % Niso )
    print( '\nALL DONE' )
