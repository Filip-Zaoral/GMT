
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
