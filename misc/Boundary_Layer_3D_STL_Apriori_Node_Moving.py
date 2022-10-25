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
nbl = 10; Hbl = 35.; gbl = 1.3
hbl = Hbl / sum(np.logspace(0,nbl - 1,nbl,True,gbl))
# Hbl = hbl * sum(np.logspace(0,nbl - 1,nbl,True,gbl))

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving")
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('Geometry.OCCImportLabels',1)                            # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
gmsh.option.setNumber('Geometry.OCCSewFaces',1)                                # Sews surfaces into shells in STEP, IGES and BRep geometries.
gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                              # Fixes shells and make solids in STEP, IGES and BRep geometries.
# gmsh.option.setNumber('Geometry.Tolerance',1e-12)
gmsh.model.occ.importShapes('Boundary_Layer_3D_IGS_Subdomain.V1.igs')
gmsh.model.occ.importShapes('Boundary_Layer_3D_IGS_Subdomain.V2.igs')
gmsh.model.occ.synchronize()

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
gmsh.write('Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')
gmsh.finalize()

# gmsh.initialize(sys.argv)
# gmsh.logger.start()
# gmsh.model.add("Boundary_Layer_3D_STL_Prep_Node_Moving")
# # gmsh.option.setNumber('Geometry.Tolerance',1e-3)
# gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
# gmsh.option.setNumber('Mesh.Algorithm',6)
# gmsh.option.setNumber('Mesh.Algorithm3D',10)
# gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                    # Removes duplicate triangles when importing STL files.
# gmsh.merge("Boundary_Layer_3D_STL_Prep_Node_Moving.stl")
# gmsh.model.geo.synchronize()
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
# gmsh.write("Tmp.Boundary_Layer_3D_STL_Prep_Node_Moving.stl")
# gmsh.logger.stop()
# gmsh.finalize()

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
# gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving.Tmp1")
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
s1 = [(2,i) for i in range(1,12)]
s2 = [(2,i) for i in range(1,3)]; s2.extend([(2,i) for i in range(12,15)])
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()

sbTags = [i[1] for i in sb]
cb = list(set(gmsh.model.getBoundary(sb,False,False,False)))
pb = list(set(gmsh.model.getBoundary(cb,False,False,False)))
cpb = cb + pb; ncpb = len(cpb)
NbTagsOld = []; NbxyzOld = []
for i in range(nsb):
    Nbxyz = gmsh.model.mesh.getNodes(sb[i][0],sb[i][1],True,False)
    NbTagsOld.extend(list(Nbxyz[0]))
    NbxyzOld.extend(list(np.reshape(Nbxyz[1],(-1,3),'C')))
for i in range(ncpb):
    Nbxyz = gmsh.model.mesh.getNodes(cpb[i][0],cpb[i][1],False,False)
    NbTags = Nbxyz[0]; Nbxyz = Nbxyz[1]; nNb = len(NbTags)
    Nbxyz = np.reshape(Nbxyz,(-1,3),'C')
    for j in range(nNb):
        EbToNb = gmsh.model.mesh.getElementsByCoordinates(Nbxyz[j,0],         \
                                                          Nbxyz[j,1],         \
                                                          Nbxyz[j,2],2)
        nEbToNb = len(EbToNb)
        NORMAL = np.zeros((nEbToNb,3),np.float64,'C')
        for k in range(nEbToNb):
            eb = gmsh.model.mesh.getElement(EbToNb[k])
            esTag = eb[3]; eb = eb[1]
            if (esTag in sbTags) and (NbTags[j] in eb):
                Nsb = []
                Nsb.append(NbxyzOld[NbTagsOld.index(eb[0])])
                Nsb.append(NbxyzOld[NbTagsOld.index(eb[1])])
                Nsb.append(NbxyzOld[NbTagsOld.index(eb[2])])
                u1 = Nsb[2] - Nsb[0]; u2 = Nsb[1] - Nsb[0]
                u1xu2 = np.cross(u1,u2)
                NORMAL[k,:] = u1xu2/np.linalg.norm(u1xu2)                      # Normalize the normal of the element.
        NORMAL = np.sum(NORMAL,0)
        NORMAL = NORMAL/np.linalg.norm(NORMAL)
        nbxyz = np.asarray([Nbxyz[j,0],Nbxyz[j,1],Nbxyz[j,2]],np.float64,'C')
        nbxyz = list(nbxyz + NORMAL*Hbl)
        try:
            gmsh.model.mesh.setNode(NbTags[j],nbxyz,[0.,0.,0.])
        except:
            pass
gmsh.model.geo.synchronize()

# sbTags = [i[1] for i in sb]
# cb = list(set(gmsh.model.getBoundary(sb,False,False,False)))
# pb = list(set(gmsh.model.getBoundary(cb,False,False,False)))
# cpb = cb + pb; ncpb = len(cpb)
# NbTagsOld = []; NbxyzOld = []
# for i in range(nsb):
#     Nbxyz = gmsh.model.mesh.getNodes(sb[i][0],sb[i][1],True,False)
#     NbTagsOld.extend(list(Nbxyz[0]))
#     NbxyzOld.extend(list(np.reshape(Nbxyz[1],(-1,3),'C')))
# for i in range(ncpb):
#     Nbxyz = gmsh.model.mesh.getNodes(cpb[i][0],cpb[i][1],False,False)
#     NbTags = Nbxyz[0]; Nbxyz = Nbxyz[1]; nNb = len(NbTags)
#     Nbxyz = np.reshape(Nbxyz,(-1,3),'C')
#     for j in range(nNb):
#         EbToNb = gmsh.model.mesh.getElementsByCoordinates(Nbxyz[j,0],         \
#                                                           Nbxyz[j,1],         \
#                                                           Nbxyz[j,2],2)
#         nEbToNb = len(EbToNb)
#         NORMAL = []; ALPHA = []
#         for k in range(nEbToNb):
#             eb = gmsh.model.mesh.getElement(EbToNb[k])
#             esTag = eb[3]; eb = eb[1]
#             if (esTag in sbTags) and (NbTags[j] in eb):
#                 Nsb = []
#                 Nsb.append(NbxyzOld[NbTagsOld.index(eb[0])])
#                 Nsb.append(NbxyzOld[NbTagsOld.index(eb[1])])
#                 Nsb.append(NbxyzOld[NbTagsOld.index(eb[2])])
#                 JJ = [2,1,0]
#                 jj = np.where(eb == NbTags[j])[0][0]
#                 if jj == 1: JJ = list(np.roll(JJ,1))
#                 JJ.remove(jj)
#                 u1 = Nsb[JJ[0]] - Nsb[jj]; u2 = Nsb[JJ[1]] - Nsb[jj]
#                 u1 = u1/np.linalg.norm(u1); u2 = u2/np.linalg.norm(u2)
#                 ALPHA.append(np.arccos(np.dot(u1,u2)))
#                 u1xu2 = np.cross(u1,u2)
#                 NORMAL.append(u1xu2/np.linalg.norm(u1xu2))
#         NORMAL = np.sum(NORMAL,0)
#         NORMAL = NORMAL/np.linalg.norm(NORMAL)
#         ALPHA = np.asarray(ALPHA,np.float64,'C')
#         NORMAL = np.asarray(NORMAL,np.float64,'C')
#         normal = np.dot(ALPHA,NORMAL)/np.sum(ALPHA)
#         normal = normal/np.linalg.norm(normal)
#         nbxyz = np.asarray([Nbxyz[j,0],Nbxyz[j,1],Nbxyz[j,2]],np.float64,'C')
#         nbxyz = list(nbxyz + NORMAL*Hbl)
#         try:
#             gmsh.model.mesh.setNode(NbTags[j],nbxyz,[0.,0.,0.])
#         except:
#             pass
# gmsh.model.geo.synchronize()

gmsh.model.removeEntities(sb,True)
sa = gmsh.model.getEntities(2); nsa = len(sa)
# ss = gmsh.model.getEntities(2); ns1 = len(ss)
# sa = [i for i in ss if i not in sb]
gmsh.model.mesh.classifySurfaces(pi,False,True,pi)                             # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.mesh.optimize('Laplace2D',True,5)
# ss = gmsh.model.getEntities(2); ns2 = len(ss)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp1.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
# gmsh.merge("Tmp1.Boundary_Layer_3D_STL_Prep_Node_Moving.stl")
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.createTopology()
# gmsh.model.geo.synchronize()
# # saTmp = gmsh.model.getEntities(2); saTmp = saTmp[nsa:]
# # caTags = [57,60,61,62,63,64]
# # caTmpTags = [78,79,80,82,83,87]
# sbTags = [i + ns1 for i in sbTags]; sb = [(2,i) for i in sbTags]
# ssTmp = gmsh.model.getEntities(2); ssTmp = ssTmp[ns2:]
# sbTmpTags = [i + ns2 for i in sbTags]
# AT = [1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
# gmsh.model.mesh.setPeriodic(2,sbTags,sbTmpTags,AT)
# gmsh.model.mesh.generate(2)
# gmsh.model.geo.synchronize()
# gmsh.model.removeEntities(ssTmp,True)
# gmsh.write("Tmp1.Boundary_Layer_3D_STL_Prep_Node_Moving.stl")
gmsh.model.remove()

gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving.Tmp2")
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.removeEntities(sa,True)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp2.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.remove()

gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving.Tmp")
# gmsh.option.setNumber('Geometry.Tolerance',1e-12)
gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.
gmsh.option.setNumber("Mesh.SaveAll",1)                                        # Force Gmsh to write only the elements belonging to a Physical Group.
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('Mesh.SurfaceEdges',1)
gmsh.option.setNumber('Mesh.SurfaceFaces',1)
gmsh.option.setNumber('Mesh.VolumeEdges',0)
gmsh.option.setNumber('Mesh.VolumeFaces',0)
gmsh.merge("Tmp2.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.mesh.classifySurfaces(pi,True,True,pi)                              # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
gmsh.model.mesh.createGeometry()                                               # Creates a geometry for all the discrete curves and surfaces in the mesh.
gmsh.model.geo.synchronize()
sb = gmsh.model.getEntities(2); nsb = len(sb)

n = np.linspace(1,1,nbl)
d = np.array([hbl * sum(np.logspace(0,i - 1,i,True,gbl)) for i in range(1,nbl + 1)]) / Hbl
NsbTags = list(set(gmsh.model.mesh.getNodes(2,-1,True,False)[0]))              # Get tags of all the nodes of the boundary layer surfaces.
nNsb = len(NsbTags)
EsbTags = list(gmsh.model.mesh.getElements(2,-1)[1][0])
nEsb = len(EsbTags)
XYZ = np.zeros((nEsb,3,3),np.float64,'C')
NORMALSN = np.zeros((nNsb,3),np.float64,'C')
NORMALSE = np.zeros((nEsb,3,3),np.float64,'C')
NsbToEsb = np.zeros((nEsb,3),np.int32,'C')
for i in range(nsb):                                                           # Loop over the boundary layer surfaces.
    Eb = gmsh.model.mesh.getElements(sb[i][0],sb[i][1])
    NbToEb = np.reshape(Eb[2],(-1,3),'C'); Eb = Eb[1][0]                       # Get tags of elements and their nodes of the current surface.
    nEb = len(Eb)
    for j in range(nEb):                                                       # Loop over the triangular elements forming the current boundary layer surface.
        idx = EsbTags.index(Eb[j])
        idx1 = NsbTags.index(NbToEb[j,0]); NsbToEsb[idx,0] = idx1
        idx2 = NsbTags.index(NbToEb[j,1]); NsbToEsb[idx,1] = idx2
        idx3 = NsbTags.index(NbToEb[j,2]); NsbToEsb[idx,2] = idx3
        NEsb1 = gmsh.model.mesh.getNode(NbToEb[j,0])[0]
        NEsb2 = gmsh.model.mesh.getNode(NbToEb[j,1])[0]
        NEsb3 = gmsh.model.mesh.getNode(NbToEb[j,2])[0]
        XYZ[idx,:,0] = NEsb1
        XYZ[idx,:,1] = NEsb2
        XYZ[idx,:,2] = NEsb3
        u1 = NEsb3 - NEsb1; u2 = NEsb2 - NEsb1                                 # Calculate two oriented vectors, each spanning across one of the sides of the triangular element.
        u1xu2 = np.cross(u1,u2)                                                # Calculate normal of the element.
        u1xu2 = u1xu2/np.linalg.norm(u1xu2)                                    # Normalize the normal of the element.
        NORMALSN[idx1,:] += u1xu2                                              # Calculate the Gouraud smoothed node normals from the elements normals.
        NORMALSN[idx2,:] += u1xu2                                              # Calculate the Gouraud smoothed node normals from the elements normals.
        NORMALSN[idx3,:] += u1xu2                                              # Calculate the Gouraud smoothed node normals from the elements normals.
XYZ = np.reshape(np.reshape(XYZ,(nEsb,1,9),'C'),(nEsb,9),'C')
NORMALSN = np.array([NORMALSN[i,:] / np.linalg.norm(NORMALSN[i,:]) * Hbl for i in range(nNsb)])
for i in range(nEsb):
    NORMALSE[i,0,:] = NORMALSN[NsbToEsb[i,0],:]
    NORMALSE[i,1,:] = NORMALSN[NsbToEsb[i,1],:]
    NORMALSE[i,2,:] = NORMALSN[NsbToEsb[i,2],:]
NORMALSE = np.reshape(np.reshape(NORMALSE,(nEsb,1,9),'C'),(nEsb,9),'C')
XYZ = np.reshape(np.concatenate((XYZ,NORMALSE),1),(1,-1),'C')[0]
blvTag = gmsh.view.add('InflationLayersNormals')
gmsh.view.addListData(blvTag,'VT',nEsb,XYZ)
blvi = gmsh.view.getIndex(blvTag)
Sb = gmsh.model.geo.extrudeBoundaryLayer(sb,n,d,True,False,blvi); nSb = len(Sb)
gmsh.model.geo.synchronize()
Sb1 = [Sb[i - 1] for i in range(nSb) if Sb[i][0] == 3]
Sb = [i for i in Sb if i[0] != 3]; Sb.extend(sb)

# generate the mesh
try:
    gmsh.model.mesh.generate(3)
except:
    pass
    
gmsh.merge("Tmp1.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
s1 = [(2,i) for i in range(5,8)]; s1.extend([(2,i) for i in range(122,130)])
s2 = [(2,8),(2,122),(2,123),(2,130),(2,131)]
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.synchronize()

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
remove('Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')
remove('Tmp1.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')
remove('Tmp2.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')