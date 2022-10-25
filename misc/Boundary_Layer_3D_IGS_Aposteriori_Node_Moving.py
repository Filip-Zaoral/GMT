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
sba = [(2,5),(2,6),(2,7),(2,12)]; nsba = len(sba)
nbl = 10; Hbl = 30.; gbl = 1.2
hbl = Hbl / sum(np.logspace(0,nbl - 1,nbl,True,gbl))
# Hbl = hbl * sum(np.logspace(0,nbl - 1,nbl,True,gbl))
s1 = [(2,i) for i in range(1,12)]
s2 = [(2,i) for i in range(1,3)]; s2.extend([(2,i) for i in range(12,15)])

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_IGS_Aposteriori_Node_Moving")
gmsh.option.setNumber('Geometry.OCCImportLabels',1)                            # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
gmsh.option.setNumber('Geometry.OCCSewFaces',1)                                # Sews surfaces into shells in STEP, IGES and BRep geometries.
gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                              # Fixes shells and make solids in STEP, IGES and BRep geometries.
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
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
gmsh.write('Tmp.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl')
gmsh.finalize()

# gmsh.initialize(sys.argv)
# gmsh.logger.start()
# gmsh.model.add("Boundary_Layer_3D_IGS_Post_Node_Moving")
# gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                    # Removes duplicate triangles when importing STL files.
# # gmsh.option.setNumber('Geometry.Tolerance',1e-3)
# gmsh.merge("Boundary_Layer_3D_IGS_Post_Node_Moving.stl")
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
# gmsh.write("Tmp.Boundary_Layer_3D_IGS_Post_Node_Moving.stl")
# gmsh.logger.stop()
# gmsh.finalize()

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.Tmp2")
gmsh.merge("Tmp.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
ss = gmsh.model.getEntities(2); ns = len(ss)
sa = [i for i in ss if i not in sba]
gmsh.model.removeEntities(sa,True)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp2.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.remove()

gmsh.model.add("Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.Tmp2")
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
gmsh.merge("Tmp2.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.mesh.classifySurfaces(pi,True,True,pi)                              # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.geo.synchronize()
gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('Mesh.SurfaceEdges',1)
gmsh.option.setNumber('Mesh.SurfaceFaces',1)
gmsh.option.setNumber('Mesh.VolumeEdges',0)
gmsh.option.setNumber('Mesh.VolumeFaces',0)
sbb = gmsh.model.getEntities(2); nsbb = len(sbb)

sb1 = [(2,5),(2,6),(2,7)]
sb2 = [(2,8)]
n = np.linspace(1,1,nbl)
d = np.array([hbl * sum(np.logspace(0,i - 1,i,True,gbl)) for i in range(1,nbl + 1)]) / Hbl
NsbTags = []; EsbTags = []; NpcbTags = []
for i in sbb:                                                                  # Loop over all boundary layer surfaces.
    NsbTags += list(gmsh.model.mesh.getNodes(2,i[1],True,False)[0])            # Get tags of all the nodes on the current boundary layer surface.
    EsbTags += list(gmsh.model.mesh.getElements(2,i[1])[1][0])                 # Get tags of all the nodes on the current boundary layer surface.
NsbTags = list(set(NsbTags)); nNsb = len(NsbTags)                              # Remove duplicates.
EsbTags = list(set(EsbTags)); nEsb = len(EsbTags)                              # Remove duplicates.
cb1 = list(gmsh.model.getBoundary(sb1,True,False,False))
cb2 = list(gmsh.model.getBoundary(sb2,True,False,False))
cb = list(set(cb1 + cb2))
pb = list(set(gmsh.model.getBoundary(cb,False,False,False)))
pcb = pb + cb; npcb = len(pcb)
for i in pcb:
    NpcbTags += list(gmsh.model.mesh.getNodes(i[0],i[1],False,False)[0])       # Get tags of all the nodes on the current boundary layer point/curve.
NpcbTags = list(set(NpcbTags)); nNpcb = len(NpcbTags)                          # Remove duplicates.
XYZ = np.zeros((nEsb,3,3),np.float64,'C')
NORMALSN = np.zeros((nNsb,3),np.float64,'C')
NORMALSE = np.zeros((nEsb,3,3),np.float64,'C')
NsbToEsb = np.zeros((nEsb,3),np.int32,'C')
pcbXYZ = np.zeros((nNpcb,6),np.float64,'C')
for i in range(nsbb):                                                          # Loop over the boundary layer surfaces.
    Eb = gmsh.model.mesh.getElements(sbb[i][0],sbb[i][1])
    NbToEb = np.reshape(Eb[2],(-1,3),'C'); Eb = Eb[1][0]                       # Get tags of elements and their nodes of the current surface.
    nEb = len(Eb)
    for j in range(nEb):                                                       # Loop over the triangular elements forming the current boundary layer surface.
        idxEN = []; NEsb = []
        idxE = EsbTags.index(Eb[j])
        for k in range(3):
            idxEN.extend([NsbTags.index(NbToEb[j,k])])
            NsbToEsb[idxE,k] = idxEN[k]
            NEsb.extend([gmsh.model.mesh.getNode(NbToEb[j,k])[0]])
            XYZ[idxE,:,k] = NEsb[k]
        u1 = NEsb[2] - NEsb[0]; u2 = NEsb[1] - NEsb[0]                         # Calculate two oriented vectors, each spanning across one of the sides of the triangular element.
        u1xu2 = np.cross(u1,u2)                                                # Calculate normal of the element.
        u1xu2 = u1xu2/np.linalg.norm(u1xu2)                                    # Normalize the normal of the element.
        NORMALSN[idxEN,:] += u1xu2                                             # Calculate the Gouraud smoothed node normals from the elements normals.
        for k in range(3):
            try:
                idxN = NpcbTags.index(NbToEb[j,k])
                pcbXYZ[idxN,:3] = NEsb[k]
                pcbXYZ[idxN,3:] += u1xu2
            except:
                pass
XYZ = np.reshape(np.reshape(XYZ,(nEsb,1,9),'C'),(nEsb,9),'C')
for i in range(nNsb):
    NORMALSN[i,:] = NORMALSN[i,:] / np.linalg.norm(NORMALSN[i,:]) * Hbl
for i in range(nNpcb):
    pcbXYZ[i,3:] = pcbXYZ[i,3:] / np.linalg.norm(pcbXYZ[i,3:]) * Hbl
for i in range(nEsb):
    for j in range(3):
        NORMALSE[i,j,:] = NORMALSN[NsbToEsb[i,j],:]
NORMALSE = np.reshape(np.reshape(NORMALSE,(nEsb,1,9),'C'),(nEsb,9),'C')
XYZ = np.reshape(np.concatenate((XYZ,NORMALSE),1),(1,-1),'C')[0]
pcbXYZ = pcbXYZ[:,:3] + pcbXYZ[:,3:]
blvTag = gmsh.view.add('InflationLayersNormals')
gmsh.view.addListData(blvTag,'VT',nEsb,XYZ)
blvi = gmsh.view.getIndex(blvTag)
Sb = gmsh.model.geo.extrudeBoundaryLayer(sbb,n,d,True,False,blvi); nSb = len(Sb)
gmsh.model.geo.synchronize()
Vb = [i for i in Sb if i[0] == 3]
Sb1 = [Sb[i - 1] for i in range(nSb) if Sb[i][0] == 3]
Sb = list(set([i for i in Sb if i[0] != 3])); Sb.extend(sbb); nSb = len(Sb)
del (NsbTags,EsbTags,XYZ,NORMALSN,NORMALSE,NsbToEsb)

# Generate mesh for the boundary layers:
try:
    gmsh.model.mesh.generate(3)
except:
    pass
gmsh.model.mesh.removeDuplicateNodes(Vb)

s1 = [(2,i) for i in range(1,12)]
s2 = [(2,1),(2,2),(2,12),(2,13),(2,14)]
gmsh.model.add("Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.Tmp1")
# gmsh.option.setNumber('Geometry.Tolerance',1e-5)
gmsh.merge("Tmp.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
sbTags = [i[1] for i in sba]
sb1 = [i for i in s1 if i in sba]
sb2 = [i for i in s2 if i in sba]
cb1 = list(gmsh.model.getBoundary(sb1,True,False,False))
cb2 = list(gmsh.model.getBoundary(sb2,True,False,False))
cb = list(set(cb1 + cb2))
pb = list(set(gmsh.model.getBoundary(cb,False,False,False)))
pcb = cb + pb; npcb = len(pcb)
for i in range(npcb):
    Nbxyz = gmsh.model.mesh.getNodes(pcb[i][0],pcb[i][1],False,False)
    NbTags = Nbxyz[0]; Nbxyz = Nbxyz[1]
    Nbxyz = np.reshape(Nbxyz,(-1,3),'C'); nNb = len(NbTags)
    for j in range(nNb):
        NbToNa = np.linalg.norm(Nbxyz[j,:] - pcbXYZ,None,1)
        jj = np.argmin(NbToNa)
        gmsh.model.mesh.setNode(NbTags[j],list(pcbXYZ[jj,:]),[0.,0.,0.])
gmsh.model.geo.synchronize()
gmsh.model.removeEntities(sba,True)
gmsh.model.mesh.classifySurfaces(pi,False,True,pi)                             # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.mesh.optimize('Laplace2D',True,5)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp1.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.remove()
gmsh.model.setCurrent("Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.Tmp2")
gmsh.merge("Tmp1.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
ss = gmsh.model.getEntities(2); ns = len(ss)
gmsh.model.mesh.removeDuplicateNodes()

# Get the "top" surfaces of the boundary layers:
s1 = [(2,i) for i in range(5,8)]; s1.extend([(2,i) for i in range(122,130)])
s2 = [(2,8),(2,122),(2,123),(2,130),(2,131)]
s1b1Tags = [55,77,99]
s2b1Tags = [121]
Sb1Tags = [i[1] for i in Sb1]
s1aTags = [i[1] for i in s1 if i not in sbb]
s2aTags = [i[1] for i in s2 if i not in sbb]
s1aTags.extend(s1b1Tags); s2aTags.extend(s2b1Tags)
sl1 = gmsh.model.geo.addSurfaceLoop(s1aTags)
sl2 = gmsh.model.geo.addSurfaceLoop(s2aTags)
V1 = gmsh.model.geo.addVolume([sl1])
V2 = gmsh.model.geo.addVolume([sl2])
gmsh.model.geo.synchronize()
V = gmsh.model.getEntities(3); nV = len(V)

# Generate mesh for the volumes:
try:
    gmsh.model.mesh.generate(3)
except:
    pass

# Generate surface and volume elements groups:
shellNames = []
shellGroups = {}
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
    gmsh.model.setPhysicalName(3,G,VName)                                      # Assigns name to the volume group.

Log = gmsh.logger.get(); gmsh.logger.stop()
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
remove('Tmp.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl')
remove('Tmp1.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl')
remove('Tmp2.Boundary_Layer_3D_IGS_Aposteriori_Node_Moving.stl')