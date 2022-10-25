from sys import platform
if platform == "win32":
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
else:
    from Lib.gmsh.api import gmsh                                              # Locates the Gmsh source library directory.
import sys
from os import remove
import numpy as np
from math import pi

gmsh.initialize(sys.argv)
gmsh.logger.start()
s1 = [(2,9),(2,10),(2,11),(2,12),(2,13)]
ms = 0.05; nbl = 5; hbl = 25e-3; gbl = 1.5
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt")
gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                    # Removes duplicate triangles when importing STL files.
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
gmsh.merge("Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.mesh.classifySurfaces(pi/4,True,True,pi)                            # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.geo.synchronize()
gmsh.option.setNumber('Mesh.SurfaceEdges',1)
gmsh.option.setNumber('Mesh.SurfaceFaces',1)
gmsh.option.setNumber('Mesh.VolumeEdges',0)
gmsh.option.setNumber('Mesh.VolumeFaces',0)
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('General.NumThreads',8)                                  # Sets the maximum number of CPU threads to use for 2D/3D meshing.

shellNames = []
shellGroups = {}
ss = gmsh.model.getEntities(2)
ns = len(ss)
for i in range(ns):
    sName = gmsh.model.getEntityName(ss[i][0],ss[i][1])                        # Attempt to get entity labels read from the content of the IGES files.
    sName = sName.split('/')                                                   # Extracts the names of BC groups from entity labels.
    if (len(sName) in [0,1]) or ((len(sName) > 1) and  (len(sName[1]) == 0)):
        shellNames.append("A" + str(i + 1))                                    # Extracts names of surface BC groups from model part names.
    else:
        nName = shellNames.count(sName[1])
        if nName == 0:
            shellNames.append(sName[1])
        else:
            shellNames.append(sName[1] + "_" + str(i + 1))
for p in ss:
    sName = shellNames[ss.index(p)]                                            # Get entity labels read from the name of the STL files.
    gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")       # Prints names of all surface entities with successfuly identified labels/names of BC groups.
    if sName not in shellGroups:
        shellGroups[sName] = []
        shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
for sName,sTag in shellGroups.items():
    g = gmsh.model.addPhysicalGroup(2,sTag,shellNames.index(sName) + 1)        # Creates boundary surface groups.
    gmsh.model.setPhysicalName(2,g,sName)                                      # Assigns names to the boundary surface groups.

# generate curved mesh?
try:
    gmsh.model.mesh.generate(2)
except:
    pass
gmsh.model.geo.synchronize()

s1Tags = [i[1] for i in s1]
c1 = gmsh.model.getBoundary(s1,False,False,False)
c1Tags = list(set([i[1] for i in c1]))
c1 = [(1,i) for i in c1Tags]
p1 = gmsh.model.getBoundary(c1,False,False,False)
p1Tags = list(set([i[1] for i in p1]))
p1 = [(0,i) for i in p1Tags]
scp1 = []
scp1.extend([(2,i) for i in s1Tags])
scp1.extend([(1,i) for i in c1Tags])
scp1.extend([(0,i) for i in p1Tags])
nscp1 = len(scp1)
gmsh.model.mesh.renumberNodes(); gmsh.model.mesh.renumberElements()
N1TagsOld = []; N1CoordsOld = []
for i in range(nscp1):
    N1Coords = gmsh.model.mesh.getNodes(scp1[i][0],scp1[i][1],False,False)
    N1TagsOld.extend(list(N1Coords[0]))
    N1CoordsOld.extend(list(np.reshape(N1Coords[1],(-1,3),'C')))
for i in range(nscp1):
    N1Coords = gmsh.model.mesh.getNodes(scp1[i][0],scp1[i][1],False,False)
    N1Tags = N1Coords[0]; nN1 = len(N1Tags)
    N1Coords = np.reshape(N1Coords[1],(-1,3),'C')   
    for j in range(nN1):
        E1ToN1 = gmsh.model.mesh.getElementsByCoordinates(N1Coords[j,0],      \
                                                          N1Coords[j,1],      \
                                                          N1Coords[j,2],2)
        nE1ToN1 = len(E1ToN1)
        ALPHA = []; NORMAL = []
        for k in range(nE1ToN1):
            e1 = gmsh.model.mesh.getElement(E1ToN1[k])
            esTag = e1[3]; e1 = e1[1]
            if esTag in s1Tags:
                Ns1 = []
                JJ = [2,1,0]
                jj = np.where(e1 == N1Tags[j])[0][0]
                if jj == 1: JJ = list(np.roll(JJ,1))
                JJ.remove(jj)
                Ns1.append(N1CoordsOld[N1TagsOld.index(e1[0])])
                Ns1.append(N1CoordsOld[N1TagsOld.index(e1[1])])
                Ns1.append(N1CoordsOld[N1TagsOld.index(e1[2])])
                PP1 = Ns1[JJ[0]] - Ns1[jj]
                PP2 = Ns1[JJ[1]] - Ns1[jj]
                u1 = PP1/np.linalg.norm(PP1)
                u2 = PP2/np.linalg.norm(PP2)
                ALPHA.append(np.arccos(np.dot(u1,u2)))
                u1xu2 = np.cross(u1,u2)
                NORMAL.append(u1xu2/np.linalg.norm(u1xu2))
        ALPHA = np.asarray(ALPHA,np.float64,'C')
        NORMAL = np.asarray(NORMAL,np.float64,'C')
        normal = np.dot(ALPHA,NORMAL)/np.sum(ALPHA)
        normal = normal/np.linalg.norm(normal)
        n2Coords = np.asarray([N1Coords[j,0],N1Coords[j,1],N1Coords[j,2]],\
                              np.float64,'C')
        n2Coords = list(n2Coords + normal*hbl)
        try:
            gmsh.model.mesh.setNode(N1Tags[j],n2Coords,[0.,0.,0.])
        except:
            pass
gmsh.model.geo.synchronize()

Log = gmsh.logger.get(); gmsh.logger.stop()
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',2)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write('Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt.stl')
gmsh.finalize()

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt")
gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                    # Removes duplicate triangles when importing STL files.
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.mesh.classifySurfaces(pi,True,True,pi)                              # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.geo.synchronize()
gmsh.option.setNumber('Mesh.SurfaceEdges',1)
gmsh.option.setNumber('Mesh.SurfaceFaces',1)
gmsh.option.setNumber('Mesh.VolumeEdges',0)
gmsh.option.setNumber('Mesh.VolumeFaces',0)
gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
s1 = [(2,18),(2,19),(2,20),(2,21),(2,22)]

shellNames = []
shellGroups = {}
ss = gmsh.model.getEntities(2)
ns = len(ss)
for i in range(ns):
    sName = gmsh.model.getEntityName(ss[i][0],ss[i][1])                        # Attempt to get entity labels read from the content of the IGES files.
    sName = sName.split('/')                                                   # Extracts the names of BC groups from entity labels.
    if (len(sName) in [0,1]) or ((len(sName) > 1) and (len(sName[1]) == 0)):
        shellNames.append("A" + str(i + 1))                                    # Extracts names of surface BC groups from model part names.
    else:
        nName = shellNames.count(sName[1])
        if nName == 0:
            shellNames.append(sName[1])
        else:
            shellNames.append(sName[1] + "_" + str(i + 1))
for p in ss:
    sName = shellNames[ss.index(p)]                                             # Get entity labels read from the name of the STL files.
    gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")        # Prints names of all surface entities with successfuly identified labels/names of BC groups.
    if sName not in shellGroups:
        shellGroups[sName] = []
        shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
for sName,sTag in shellGroups.items():
    g = gmsh.model.addPhysicalGroup(2,sTag,                                   \
                                    shellNames.index(sName) + 1)               # Creates boundary surface groups.
    gmsh.model.setPhysicalName(2,g,sName)                                      # Assigns names to the boundary surface groups.

# create a boundary layer for all the surfaces through extrusion using the
# built-in CAD kernel: this creates topological entities that will be filled
# with a discrete geometry (a mesh extruded along the boundary normals) during
# mesh generation; in 2D more general boundary layer meshing constraints are
# also available through the BoundaryLayer Field - see
# 'naca_boundary_layer_2d.py'.
tbl = hbl / sum(np.logspace(0,nbl - 1,nbl,True,gbl))
n = np.linspace(1,1,nbl)
d = np.array([hbl - (tbl * sum(np.logspace(0,i - 1,i,True,gbl))) for i in range(nbl - 1,-1,-1)])
extbl = gmsh.model.geo.extrudeBoundaryLayer(s1,n,d,True,True)
gmsh.model.geo.synchronize()

bs = [i[1] for i in ss]
sl = gmsh.model.geo.addSurfaceLoop(bs)
v = gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

# generate the mesh
try:
    gmsh.model.mesh.generate(3)
except:
    pass

Log = gmsh.logger.get(); gmsh.logger.stop()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
remove('Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving_Alt.stl')