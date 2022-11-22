
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

from sys import argv, platform
if platform == "win32":
    from lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
else:
    import gmsh                                                                # Locates the Gmsh source library directory.

gmsh.initialize(argv)
gmsh.option.setNumber('Mesh.ColorCarousel',2)                                  # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
gmsh.option.setNumber('Mesh.SurfaceEdges',1)                                   # Display edges of surface mesh?
gmsh.option.setNumber('Mesh.SurfaceFaces',1)                                   # Display faces of surface mesh?
gmsh.option.setNumber('Mesh.VolumeEdges',0)                                    # Display edges of volume mesh?
gmsh.option.setNumber('Mesh.VolumeFaces',0)                                    # Display faces of volume mesh?
if ('-nopopup' not in argv):
    gmsh.fltk.run()
gmsh.finalize()
