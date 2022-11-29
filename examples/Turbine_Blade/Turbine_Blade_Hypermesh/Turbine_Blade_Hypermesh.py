
# GMT - Copyright (C) 2022 Filip Zaoral, IT4Innovations,
#                          VSB-Technical University of Ostrava, Czech Republic

# This file is a part of GMT.

# See the LICENSE.txt file in the GMT root directory for license information.

# GMT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.

# GMT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Initialization:
from glob import glob
from sys import path, argv, platform
import os
from shutil import copy
from subprocess import run

# Pathfinding:
if platform == "win32":
    GMT = os.getenv('GMTPATH').lower()
    CWD = os.path.dirname(os.path.abspath(__file__)).lower()
else:
    GMT = os.getenv('GMTPATH')
    CWD = os.path.dirname(os.path.abspath(__file__))
os.path.dirname(argv[0])                                                       # Locates the current directory as the working directory.
if not GMT:
    raise Exception("Error: No environment variable named 'GMTPATH' was found")
path.append(GMT)                                                               # Locates the Gmsh Meshing Tool library directory.
if CWD == GMT:
    from lib import GMTImport
    if __name__ == '__main__':
        GMTImport.Import()
else:
    FileName = os.path.basename(os.path.abspath(__file__)).replace('.py','.gcf')
    
    # Copy the Gmsh Configuration File to the Gmsh Meshing Tool directory:
    copy(FileName,os.path.join(GMT,'GMT_Configuration_File.gcf'))
    
    # Write the working directory path to the Gmsh Configuration File:
    ModelName = glob('*.py')[0].replace('.py','')
    with open(os.path.join(GMT,'GMT_Configuration_File.gcf'),'a',             \
              encoding = 'utf8') as GCF:                                       # Opens the Gmsh Configuration File for writing.
        GCF.write('\n\nWorkingDirectoryPath       = ' + CWD + ' # Absolute pa'\
                  'th to the working directory.')                              # Writes the current directory path into the Gmsh Configuration File.
        GCF.write('\nName                       = ' + ModelName + ' # Model ' \
                  'Name (IGS: <Model Name>.<Volume Name>.igs, STL: <Model Nam'\
                  'e>.<Volume Name>.<Surface Name>.stl). [string]')            # Writes the model name into the Gmsh Configuration File.
    
    # Launch of the Gmsh Meshing Tool Launcher:
    if platform == "win32":
        # import GMTLaunch; GMTLaunch
        os.startfile(os.path.join(GMT,'GMTLaunch.py'))
    else:
        run("python3 " + os.path.join(GMT,'GMTLaunch.py %s'),shell=True)
