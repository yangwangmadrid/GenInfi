# GenInfi:
#   A handy program for generating structures of generalized infinitenes,
#   generalized kekulenes, as well as clarenes
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
#    -  Completed the first decent, stable version of codes
#

VERSION = '1.0 (Mar. 2023)'
 
import re
import sys
sys.stdout.flush()
import os.path
import numpy as np
from equiangular_hexagon import * 
from kekulene import *
from clarene import *
from k_infinitene import *
from c_infinitene import *
from cc_infinitene import *
from topo import *


def printUsage():
    print( 'Usage: python geninfi.py <kinf|cinf|ccinf|kek|clr> <parameters>' )
    print( 'Parameters:' )
    print( '    Nring' )
    print( ' or h1 h2 h3 h4 h5 h6 (for ccinf/kek/clr)' )
    print( ' or h1 h2 h3 h4 h5 h6 d (for ccinf)' )
    print( ' or h1 h2 h3 h4 h5 h6 k1 k2 k3 k4 k5 k6 (for kinf/cinf)' )
    print( ' or h1 h2 h3 h4 h5 h6 k1 k2 k3 k4 k5 k6 d (for kinf/cinf)' )

#==============================================================================
#  THE MAIN PROGRAM:
#==============================================================================
 
# Welcome message:
print( 'GenInfi version', VERSION )
print( '  -- A program for generating infinitenes, kekulenes and clarenes' )
print( 'Written by Yang WANG (yangwang@yzu.edu.cn)' )
print( 'Copyright (C) 2023 Yang Wang' )
print( '' )
 
if len( sys.argv ) < 3:
    printUsage()
    exit(1)
elif len( sys.argv ) == 3:
    jobType = 'enum_Nring'
elif len( sys.argv ) == 8:
    jobType = 'enum_k'
elif len( sys.argv ) == 9:
    jobType = 'enum_kd'
elif len( sys.argv ) == 14:
    jobType = 'enum_hk'
elif len( sys.argv ) == 15:
    jobType = 'enum_hkd'
else:
    print( 'ERROR: Incorrect number of arguments' )
    printUsage()
    exit(1)


#-------- Parse command-line arguments: --------
molType = sys.argv[1].lower()
if molType != 'kinf' and molType != 'cinf' and molType != 'ccinf' \
        and molType != 'kek' and molType != 'clr':
    print( 'ERROR: Type of molecules must be kinf, cinf, ccinf, kek, or clr' )
    printUsage()
    exit(1)

# Check compatibility of molType and jobType:
if jobType == 'enum_k':
    if molType == 'kinf' or molType == 'cinf':
        print( 'ERROR: Not enough parameters for generating %s' % molType )
        printUsage()
        exit(1)
    else:
        kk = [ int(x) for x in sys.argv[2:8] ]
        # Check validity of kk[]:
        if not assert_equiangular_hexagon( kk ):
            print( 'ERROR: The specified side lengths do not form '\
                    'an equiangular hexagon' )
            exit(1)
elif jobType == 'enum_kd':
    if molType == 'kinf' or molType == 'cinf':
        print( 'ERROR: Not enough parameters for generating %s' % molType )
        printUsage()
        exit(1)
    elif molType == 'kek' or molType == 'clr':
        print( 'ERROR: Too many parameters for generating %s' % molType )
        printUsage()
        exit(1)
    else:
        kk = [ int(x) for x in sys.argv[2:8] ]
        print( kk )
        if not assert_equiangular_hexagon( kk ):
            print( 'ERROR: The specified side lengths do not form '\
                    'an equiangular hexagon' )
            exit(1)
        d = int( sys.argv[-1] )

elif jobType == 'enum_hk' or jobType == 'enum_hkd':
    if molType == 'ccinf' or molType == 'kek' or molType == 'clr':
        print( 'ERROR: Too many parameters for generating %s' % molType )
        printUsage()
        exit(1)
    else:
        hh = [ int(x) for x in sys.argv[2:8] ]
        if not assert_equiangular_hexagon( hh ):
            print( 'ERROR: The first 6 specified side lengths do not form '\
                    'an equiangular hexagon' )
            exit(1)
        kk = [ int(x) for x in sys.argv[8:14] ]
        if not assert_equiangular_hexagon( kk ):
            print( 'ERROR: The last 6 specified side lengths do not form '\
                    'an equiangular hexagon' )
            exit(1)
        if jobType == 'enum_hkd':
            d = int( sys.argv[-1] )

elif jobType == 'enum_Nring':
    Nring = int( sys.argv[2] )
    if molType == 'kek' or molType == 'clr':
        if Nring < 6:
            print( 'ERROR: Number of rings (Nring) must be at least 6' )
            exit(1)
    else:
        if Nring < 12:
            print( 'ERROR: Number of rings (Nring) must be at least 12' )
            exit(1)
    if ( molType == 'cinf' or molType == 'ccinf' or molType == 'clr' ) \
            and Nring % 2 == 1:
        print( 'ERROR: Number of rings (Nring) must be an even number for %s' \
                % molType )
        exit(1)
#-----------------------------------------------


if molType == 'kek':
    print( '*'*60 )
    if jobType == 'enum_Nring':
        print( 'ENUMERATING ALL GENERALIZED KEKULENES WITH %i RINGS' % Nring )
        print()
        enum_generalized_kekulene( Nring )
    elif jobType == 'enum_k':
        print( 'GENERATING GENERALIZED KEKULENE [%i,%i,%i,%i,%i,%i]' \
                % (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5]) )
        print()
        generalized_kekulene( kk, True, True )
    print( '*'*60 )
elif molType == 'clr':
    print( '*'*60 )
    if jobType == 'enum_Nring':
        print( 'ENUMERATING ALL GENERALIZED CLARENES WITH %i RINGS' % Nring )
        print()
        enum_clarene( Nring )
    elif jobType == 'enum_k':
        # Check that all kk[] are even numbers:
        for k in kk:
            if k % 2 != 0:
                print( 'ERROR: All specified side lengths must be even '
                      'numbers for clarenes' )
                exit(1)
        print( 'GENERATING GENERALIZED CLARENE <%i,%i,%i,%i,%i,%i>' \
                % (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5]) )
        print()
        clarene( kk, True, True )
    print( '*'*60 )
elif molType == 'kinf':
    print( '*'*60 )
    if jobType == 'enum_Nring':
        print( 'ENUMERATING ALL K-INFINITENES WITH %i RINGS' % Nring )
        print()
        enum_Nring_all_k_infinitene( Nring )
    elif jobType == 'enum_hk' or jobType == 'enum_hkd':
        if jobType == 'enum_hk':
            print( 'GENERATING ALL K-INFINITENES WITH KEKULENES'\
                    '[%i,%i,%i,%i,%i,%i] AND [%i,%i,%i,%i,%i,%i]' \
                    % (hh[0],hh[1],hh[2],hh[3],hh[4],hh[5], \
                    kk[0],kk[1],kk[2],kk[3],kk[4],kk[5]) )
            print()
            enum_hk_all_k_infinitene( hh, kk )
        if jobType == 'enum_hkd':
            print( 'GENERATING K-INFINITENES '\
                    '[%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i](%i)' \
                    % (hh[0],hh[1],hh[2],hh[3],hh[4],hh[5], \
                    kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
            print()
            k_infinitene( hh, kk, d, True )
    print( '*'*60 )
elif molType == 'cinf':
    print( '*'*60 )
    if jobType == 'enum_Nring':
        print( 'ENUMERATING ALL C-INFINITENES WITH %i RINGS' % Nring )
        print()
        enum_Nring_all_c_infinitene( Nring )
    elif jobType == 'enum_hk' or jobType == 'enum_hkd':
        # Check that all hh[] and kk[] are even numbers:
        for k in hh + kk:
            if k % 2 != 0:
                print( 'ERROR: All specified side lengths must be even '
                      'numbers for C-infinitenes' )
                exit(1)
        if jobType == 'enum_hk':
            print( 'GENERATING ALL C-INFINITENES WITH CLARENES'\
                    '<%i,%i,%i,%i,%i,%i> AND <%i,%i,%i,%i,%i,%i>' \
                    % (hh[0],hh[1],hh[2],hh[3],hh[4],hh[5], \
                    kk[0],kk[1],kk[2],kk[3],kk[4],kk[5]) )
            print()
            enum_hk_all_c_infinitene( hh, kk )
        if jobType == 'enum_hkd':
            print( 'GENERATING C-INFINITENES '\
                    '<%i,%i,%i,%i,%i,%i|%i,%i,%i,%i,%i,%i>(%i)' \
                    % (hh[0],hh[1],hh[2],hh[3],hh[4],hh[5], \
                    kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
            print()
            c_infinitene( hh, kk, d, True )
    print( '*'*60 )
elif molType == 'ccinf':
    print( '*'*60 )
    if jobType == 'enum_Nring':
        print( 'ENUMERATING ALL C-INFINITENES WITH %i RINGS' % Nring )
        print()
        enum_Nring_all_cc_infinitene( Nring )
    elif jobType == 'enum_k' or jobType == 'enum_kd':
        # Check that all kk[] are even numbers:
        for k in kk:
            if k % 2 != 0:
                print( 'ERROR: All specified side lengths must be even '
                      'numbers for CC-infinitenes' )
                exit(1)
        if jobType == 'enum_k':
            print( 'GENERATING ALL CC-INFINITENES WITH CORONENE'\
                    'AND CLARENE <%i,%i,%i,%i,%i,%i>' \
                    % (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5]) )
            print()
            enum_k_all_cc_infinitene( kk )
        if jobType == 'enum_kd':
            print( 'GENERATING CC-INFINITENES '\
                    '<1,1,1,1,1,1|%i,%i,%i,%i,%i,%i>(%i)' \
                    % (kk[0],kk[1],kk[2],kk[3],kk[4],kk[5], d) )
            print()
            cc_infinitene( kk, d, True )
    print( '*'*60 )


print( 'ALL DONE' )
#==============================================================================
