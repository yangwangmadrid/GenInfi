# write_xyz: Write Cartesian coordinates of atoms to an external *.xyz file
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
def write_xyz( outp, coord, elem, title='' ):
    Natom = coord.shape[0] # Number of all atoms

    with open( outp, 'w' ) as writer:
        writer.write( '%3.0f' % Natom )
        writer.write( '\n%s\n' % title )
        for i in range( Natom ):
            writer.write( '%3s  %12.6f %12.6f %12.6f\n' % \
                    (elem[i], coord[i,0],coord[i,1],coord[i,2]) )

    print( 'Coordinates written to file %s' % outp )
