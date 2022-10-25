from sys import platform
if platform == "win32":
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
else:
    from Lib.gmsh.api import gmsh                                              # Locates the Gmsh source library directory.
import sys
from os import remove
import numpy as np
from math import pi

ms = 40.
nV = 2
sb = [(2, 1), (2, 2), (2, 3), (2, 6), (2, 8), (2, 15), (2, 16), (2, 17), (2, 20), (2, 21), (2, 22), (2, 23), (2, 25), (2, 26), (2, 27), (2, 29), (2, 30), (2, 31), (2, 32), (2, 34), (2, 35), (2, 36), (2, 37), (2, 39), (2, 40), (2, 41), (2, 42), (2, 43), (2, 44), (2, 46), (2, 49), (2, 50), (2, 51), (2, 52), (2, 53), (2, 55), (2, 56), (2, 57), (2, 59), (2, 60)]
nsb = len(sb)
nbl = 3; Hbl = 2. ; gbl = 1.
hbl = Hbl / sum(np.logspace(0,nbl - 1,nbl,True,gbl))
# Hbl = hbl * sum(np.logspace(0,nbl - 1,nbl,True,gbl))

gmsh.initialize(sys.argv)
gmsh.logger.start()
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving")
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',True)
gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',30)                         # using the value as the target number of elements per 2*Pi radians.
gmsh.option.setNumber('Mesh.MinimumCircleNodes',6)
gmsh.option.setNumber('Mesh.MinimumCurveNodes',6)
gmsh.option.setNumber('Mesh.Optimize',False)
gmsh.option.setNumber('Mesh.OptimizeNetgen',False)
gmsh.option.setNumber('Mesh.Algorithm',6)
gmsh.option.setNumber('Mesh.Algorithm3D',10)
gmsh.option.setNumber('Geometry.OCCImportLabels',1)                            # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
gmsh.option.setNumber('Geometry.OCCSewFaces',1)                                # Sews surfaces into shells in STEP, IGES and BRep geometries.
gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                              # Fixes shells and make solids in STEP, IGES and BRep geometries.
gmsh.option.setNumber('Geometry.Tolerance',1e-6)
gmsh.model.occ.importShapes('Aircraft_Wing_Single_Flap_SpaceClaim_2.Domain.igs')
gmsh.model.occ.importShapes('Aircraft_Wing_Single_Flap_SpaceClaim_2.Trailing_Edge_Subdomain.igs')
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
gmsh.model.add("Boundary_Layer_3D_STL_Apriori_Node_Moving.Tmp1")
gmsh.merge("Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
gmsh.model.geo.synchronize()
gmsh.model.mesh.createTopology()
gmsh.model.geo.synchronize()
gmsh.model.removeEntities(sb,True)
sa = gmsh.model.getEntities(2); nsa = len(sa)
gmsh.option.setNumber('Mesh.Format',27)                                        # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                          # Sets the Gmsh to save the shell groups to the output format.
gmsh.option.setNumber('Mesh.Binary',1)                                         # Write mesh files in binary format (if possible)?
gmsh.write("Tmp1.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl")
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
gmsh.option.setNumber('Geometry.Tolerance',1e-9)
gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                                  # Mesh only entities that have no existing mesh.
gmsh.option.setNumber("Mesh.SaveAll",1)                                        # Force Gmsh to write only the elements belonging to a Physical Group.
gmsh.option.setNumber('Mesh.MeshSizeMax',ms)
gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',True)
gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',20)                         # using the value as the target number of elements per 2*Pi radians.
gmsh.option.setNumber('Mesh.MinimumCircleNodes',0)
gmsh.option.setNumber('Mesh.MinimumCurveNodes',0)
gmsh.option.setNumber('Mesh.Optimize',False)
gmsh.option.setNumber('Mesh.OptimizeNetgen',False)
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
# ALPHA = np.zeros((nNsb,1),np.float64,'C')
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
        u1 = NEsb3 - NEsb1; u2 = NEsb1 - NEsb2; u3 = NEsb2 - NEsb3             # Calculate two oriented vectors, each spanning across one of the sides of the triangular element.
        u1 = u1/np.linalg.norm(u1)
        u2 = u2/np.linalg.norm(u2)
        u3 = u3/np.linalg.norm(u3)
        u1xu2 = np.cross(u1,-u2)                                                # Calculate normal of the element.
        u1xu2 = u1xu2/np.linalg.norm(u1xu2)                                    # Normalize the normal of the element.
        a1 = np.arccos(np.dot(u1,-u2))
        a2 = np.arccos(np.dot(u2,-u3))
        a3 = np.arccos(np.dot(u3,-u1))
        NORMALSN[idx1,:] += a1 * u1xu2                                         # Calculate the Gouraud smoothed node normals from the elements normals.
        NORMALSN[idx2,:] += a2 * u1xu2                                         # Calculate the Gouraud smoothed node normals from the elements normals.
        NORMALSN[idx3,:] += a3 * u1xu2                                         # Calculate the Gouraud smoothed node normals from the elements normals.
NORMALSN = NORMALSN / np.linalg.norm(NORMALSN,None,1,True) * Hbl
XYZ = np.reshape(XYZ,(nEsb,1,9),'C')[:,0,:]
for i in range(nEsb):
    NORMALSE[i,0,:] = NORMALSN[NsbToEsb[i,0],:]
    NORMALSE[i,1,:] = NORMALSN[NsbToEsb[i,1],:]
    NORMALSE[i,2,:] = NORMALSN[NsbToEsb[i,2],:]
NORMALSE = np.reshape(NORMALSE,(nEsb,1,9),'C')[:,0,:]
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

Log = gmsh.logger.get(); gmsh.logger.stop()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
remove('Tmp.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')
remove('Tmp1.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')
remove('Tmp2.Boundary_Layer_3D_STL_Apriori_Node_Moving.stl')