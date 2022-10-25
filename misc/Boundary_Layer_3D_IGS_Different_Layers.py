from sys import platform
if platform == "win32":
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
else:
    from Lib.gmsh.api import gmsh                                              # Locates the Gmsh source library directory.
import sys
from os import remove
import numpy as np
from math import pi

ms = 50.
nV = 2
sb = [(2,5),(2,6),(2,7),(2,12)]; nsb = len(sb)
nbl1 = 10; hbl1 = 20.; gbl1 = 1.2
nbl2 = 5

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_STL_Node_Moving")
gmsh.option.setNumber('Geometry.OCCImportLabels',1)                            # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
gmsh.option.setNumber('Geometry.OCCSewFaces',1)                                # Sews surfaces into shells in STEP, IGES and BRep geometries.
gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                              # Fixes shells and make solids in STEP, IGES and BRep geometries.
# gmsh.option.setNumber('Geometry.Tolerance',1e-3)
gmsh.model.occ.importShapes('Boundary_Layer_3D_IGS_Subdomain.V1.igs')
gmsh.model.occ.importShapes('Boundary_Layer_3D_IGS_Subdomain.V2.igs')
gmsh.model.occ.synchronize()
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)

k = 1
for i in range(1,3):
    for j in range(1+k,4):
        try:
            VV, _ = gmsh.model.occ.fragment([(3,i)],[(3,j)],-1,True,True)
            gmsh.model.occ.synchronize()
        except:
            pass

ss = gmsh.model.getEntities(2); ns = len(ss)
V = gmsh.model.getEntities(3); nV = len(V)
shellNames = []
shellGroups = {}
for i in range(ns):
    shellNames.append("A" + str(i + 1))                                        # Extracts names of surface BC groups from model part names.
for p in ss:
    sName = shellNames[ss.index(p)]                                            # Get entity labels read from the name of the STL files.
    gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")       # Prints names of all surface entities with successfuly identified labels/names of BC groups.
    if sName not in shellGroups:
        shellGroups[sName] = []
        shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
for sName,sTag in shellGroups.items():
    g = gmsh.model.addPhysicalGroup(2,sTag,shellNames.index(sName) + 1)        # Creates boundary surface groups.
    gmsh.model.setPhysicalName(2,g,sName)                                      # Assigns names to the boundary surface groups.

try:
    gmsh.model.mesh.generate(2)
except:
    pass

gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',2)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write('Tmp.Boundary_Layer_3D_STL_Node_Moving.stl')
gmsh.finalize()

# gmsh.initialize(sys.argv)
# gmsh.logger.start()
# gmsh.model.add("Boundary_Layer_3D_STL_Node_Moving")
# gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                    # Removes duplicate triangles when importing STL files.
# # gmsh.option.setNumber('Geometry.Tolerance',1e-3)
# gmsh.merge("Boundary_Layer_3D_STL_Node_Moving.stl")
# gmsh.model.geo.synchronize()
# gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
# gmsh.option.setNumber('Mesh.Algorithm',6)
# gmsh.option.setNumber('Mesh.Algorithm3D',10)
# gmsh.model.mesh.removeDuplicateNodes()
# gmsh.model.mesh.createTopology()
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.classifySurfaces(pi/4,True,True,pi)                            # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
# gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
# gmsh.model.geo.synchronize()

# ss = gmsh.model.getEntities(2); ns = len(ss)
# shellNames = []; shellGroups = {}
# for i in range(ns):
#     sName = gmsh.model.getEntityName(ss[i][0],ss[i][1])                        # Attempt to get entity labels read from the content of the IGES files.
#     nName = shellNames.count(sName[1])
#     if nName == 0:
#         shellNames.append(sName[1])
#     else:
#         shellNames.append(sName[1] + "_" + str(i + 1))
# for p in ss:
#     sName = shellNames[ss.index(p)]                                            # Get entity labels read from the name of the STL files.
#     gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")       # Prints names of all surface entities with successfuly identified labels/names of BC groups.
#     if sName not in shellGroups:
#         shellGroups[sName] = []
#         shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
# for sName,sTag in shellGroups.items():
#     g = gmsh.model.addPhysicalGroup(2,sTag,shellNames.index(sName) + 1)        # Creates boundary surface groups.
#     gmsh.model.setPhysicalName(2,g,sName)                                      # Assigns names to the boundary surface groups.

# try:
#     gmsh.model.mesh.generate(2)
# except:
#     pass

# gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
# gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',2)                          # Sets the Gmsh to save the shell groups to the output format.
# gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
# gmsh.write("Tmp.Boundary_Layer_3D_STL_Node_Moving.stl")
# gmsh.logger.stop()
# gmsh.finalize()

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_STL_Node_Moving.Tmp1")
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.removeEntities(sb,True)
sa = gmsh.model.getEntities(2); nsa = len(sa)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp1.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.remove()

gmsh.model.add("Boundary_Layer_3D_STL_Node_Moving.Tmp2")
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.removeEntities(sa,True)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp2.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.remove()

gmsh.model.add("Boundary_Layer_3D_STL_Node_Moving.Tmp")
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
gmsh.merge("Tmp2.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.mesh.classifySurfaces(pi,True,True,pi)                              # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.geo.synchronize()
sb = gmsh.model.getEntities(2); nsb = len(sb)
gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.
gmsh.option.setNumber("Mesh.SaveAll",1)                                        # Force Gmsh to write only the elements belonging to a Physical Group.
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('Mesh.SurfaceEdges',1)
gmsh.option.setNumber('Mesh.SurfaceFaces',1)
gmsh.option.setNumber('Mesh.VolumeEdges',0)
gmsh.option.setNumber('Mesh.VolumeFaces',0)

# create a boundary layer for all the surfaces through extrusion using the
# built-in CAD kernel: this creates topological entities that will be filled
# with a discrete geometry (a mesh extruded along the boundary normals) during
# mesh generation; in 2D more general boundary layer meshing constraints are
# also available through the BoundaryLayer Field - see
# 'naca_boundary_layer_2d.py'.
hbl1 = hbl1 / sum(np.logspace(0,nbl1 - 1,nbl1,True,gbl1))
n1 = np.linspace(1,1,nbl1)
d1 = np.array([-hbl1 * sum(np.logspace(0,i - 1,i,True,gbl1)) for i in range(1,nbl1 + 1)])
SbI = gmsh.model.geo.extrudeBoundaryLayer(sb,n1,d1,True); nSbI = len(SbI)
gmsh.model.geo.synchronize()
Sb1I = [SbI[i - 1] for i in range(nSbI) if SbI[i][0] == 3]
SbI = [i for i in SbI if i[0] != 3]

# generate the mesh
try:
    gmsh.model.mesh.generate(3)
except:
    pass

sbII = [(2,121)]
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.synchronize()
hbl2 = hbl1*gbl1**(nbl1-1)
n2 = np.linspace(1,1,nbl2 - 1)
d2 = np.array([-hbl2 * sum(np.logspace(0,i - 1,i,True,gbl1)) for i in range(1,nbl2 + 1)]); d2 = d2[1:]
SbII = gmsh.model.geo.extrudeBoundaryLayer(sbII,n2,d2,True); nSbII = len(SbII)
gmsh.model.geo.synchronize()
Sb1II = [SbII[i - 1] for i in range(nSbII) if SbII[i][0] == 3]
SbII = [i for i in SbII if i[0] != 3]

# generate the mesh
try:
    gmsh.model.mesh.generate(3)
except:
    pass

gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.synchronize()
gmsh.merge("Tmp1.Boundary_Layer_3D_STL_Node_Moving.stl")
gmsh.model.geo.synchronize()
s1 = [(2,i) for i in range(5,8)]; s1.extend([(2,i) for i in range(122,130)])
s2 = [(2,8),(2,122),(2,123),(2,130),(2,131)]
cb = []
sb1 = [i for i in s1 if i in sb]
sb2 = [i for i in s2 if i in sb]
cb.extend(gmsh.model.getBoundary(sb1,True,False,False))
cb.extend(gmsh.model.getBoundary(sb2,True,False,False))
cb = list(set(cb))
pb = gmsh.model.getBoundary(cb,False,False,False)
pb = list(set(pb))
cpb = []; cpb.extend(cb); cpb.extend(pb); ncpb = len(cpb)

sbTags = [i[1] for i in sb]
SbTags = list(set([i[1] for i in Sb])); SbTags.extend(sbTags)
Cb1Tags = gmsh.model.getBoundary(Sb1,False,False,False)
Cb1Tags = list(set([i[1] for i in Cb1Tags]))
Cb0Tags = gmsh.model.getBoundary(sb,False,False,False)
Cb0Tags = list(set([i[1] for i in Cb0Tags]))
for i in range(ncpb):
    Nb0xyz = gmsh.model.mesh.getNodes(cpb[i][0],cpb[i][1],False,False)
    Nb0Tags = Nb0xyz[0]; Nb0xyz = Nb0xyz[1]; nNb0 = len(Nb0Tags)
    Nb0xyz = np.reshape(Nb0xyz,(-1,3),'C')
    if cpb[i][0] == 1:
        sb01 = gmsh.model.getAdjacencies(1,cpb[i][1])[0]
        sb01 = [(2,[j for j in sb01 if j not in sbTags][0])]
        cb01Tags = gmsh.model.getBoundary(sb01,True,False,False)
        cb01Tags = [j[1] for j in cb01Tags]; cb01Tags.remove(cpb[i][1])
        cpb1Tag = [j for j in cb01Tags if j in Cb1Tags][0]
    else:
        cb01 = gmsh.model.getAdjacencies(0,cpb[i][1])[0]
        cb01 = [(1,[j for j in cb01 if j not in Cb0Tags][0])]
        pb01Tags = gmsh.model.getBoundary(cb01,True,False,False)
        cb01Tags = [j[1] for j in pb01Tags]; cb01Tags.remove(cpb[i][1])
        cpb1Tag = cb01Tags[0]
    Nb1xyz = gmsh.model.mesh.getNodes(cpb[i][0],cpb1Tag,False,False)
    Nb1Tags = Nb1xyz[0]; Nb1xyz = Nb1xyz[1]; nNb1 = len(Nb1Tags)
    Nb1xyz = np.reshape(Nb1xyz,(-1,3),'C')
    for j in range(nNb0):
        Eab = gmsh.model.mesh.getElementsByCoordinates(Nb0xyz[j,0],           \
                                                       Nb0xyz[j,1],           \
                                                       Nb0xyz[j,2],2)
        nEab = len(Eab)
        NaTags = []
        for k in range(nEab):
            ea = gmsh.model.mesh.getElement(Eab[k])
            easTag = ea[3]; ea = ea[1]
            if easTag not in SbTags:
                NaTags.extend(ea)
        NaTags = list(set(NaTags)); nNa = len(NaTags)
        for k in range(nNa):
            naxyz = gmsh.model.mesh.getNode(NaTags[k])[0]
            if (abs(Nb0xyz[j,:] - naxyz) <= 1e-5).all():
                Nb1Tona = np.linalg.norm(naxyz - Nb1xyz,None,1)
                jj = np.argmin(Nb1Tona)
                gmsh.model.mesh.setNode(NaTags[k],list(Nb1xyz[jj,:]),[0.,0.,0.])
gmsh.model.geo.synchronize()
gmsh.model.mesh.removeDuplicateNodes()

# gmsh.model.mesh.optimize("Netgen")
# gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.

# get "top" surfaces of the boundary layer
Sb1Tags = [i[1] for i in Sb1]
s1b1Tags = [55,77,99]
s2b1Tags = [121]
s1aTags = [i[1] for i in s1 if i not in sb]
s2aTags = [i[1] for i in s2 if i not in sb]
s1aTags.extend(s1b1Tags); s2aTags.extend(s2b1Tags)
sl1 = gmsh.model.geo.addSurfaceLoop(s1aTags)
sl2 = gmsh.model.geo.addSurfaceLoop(s2aTags)
V1 = gmsh.model.geo.addVolume([sl1])
V2 = gmsh.model.geo.addVolume([sl2])
gmsh.model.geo.synchronize()

shellNames = []
shellGroups = {}
ss = gmsh.model.getEntities(2); ns = len(ss)
for i in range(ns):
    shellNames.append("A" + str(i + 1))                                        # Generate names of surface BC groups from model part names.
for p in ss:
    sName = shellNames[ss.index(p)]                                            # Get entity labels read from the name of the STL files.
    gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")       # Prints names of all surface entities with successfuly identified labels/names of BC groups.
    if sName not in shellGroups:
        shellGroups[sName] = []
        shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
for sName,sTag in shellGroups.items():
    g = gmsh.model.addPhysicalGroup(2,sTag,shellNames.index(sName) + 1)        # Creates boundary surface groups.
    gmsh.model.setPhysicalName(2,g,sName)                                      # Assigns names to the boundary surface groups.
solidNames = []
solidGroups = {}
V = gmsh.model.getEntities(3); nV = len(V)
for i in range(nV):
    solidNames.append("V" + str(i + 1))                                        # Generate names of surface BC groups from model part names.
for p in V:
    VName = solidNames[V.index(p)]
    if VName:
        gmsh.logger.write("Entity " + str(p) + " has label " + VName,"info")   # Prints names of all volume entities with successfuly identified labels/names of BC groups.
        if VName not in solidGroups:
            solidGroups[VName] = []
        solidGroups[VName].append(p[1])
for VName,VTag in solidGroups.items():
    G = gmsh.model.addPhysicalGroup(3,VTag,solidNames.index(VName) + 1)        # Creates volume group.
    gmsh.model.setPhysicalName(3,G,VName)                                      # Assigns name to the volume group.\

# generate the mesh
try:
    gmsh.model.mesh.generate(3)
except:
    pass

# gmsh.model.mesh.optimize('HighOrderFastCurving')
# gmsh.model.mesh.optimize('HighOrder')

Log = gmsh.logger.get(); gmsh.logger.stop()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
remove('Tmp.Boundary_Layer_3D_STL_Node_Moving.stl')
remove('Tmp1.Boundary_Layer_3D_STL_Node_Moving.stl')
remove('Tmp2.Boundary_Layer_3D_STL_Node_Moving.stl')