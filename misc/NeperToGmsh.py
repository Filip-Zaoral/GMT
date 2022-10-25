
# Intitialization:
from sys import path,argv
import os
from glob import glob
from math import pi,sqrt
index = os.path.dirname(argv[0]).find(r'\Benchmarks')
if index == -1:
    index = os.path.dirname(argv[0]).find(r'\benchmarks')
path.append(os.path.dirname(argv[0])[:index])                                  # Locates the Gmsh Meshing Tool library directory.
from Lib import writeToLogFile
from Lib.gmsh.lib import gmsh                                                  # Locates the Gmsh source library directory.
from time import perf_counter
import numpy as np
from collections import deque
os.path.dirname(argv[0])                                                   # Locates the working (current) directory.
gmsh.initialize(argv)
gmsh.logger.start()
TPreP = perf_counter()

# Options:
gmsh.option.setNumber('Mesh.Binary',1)                                     # Write mesh files in binary format (if possible)?

# Geometry import:
gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                # Removes duplicate triangles when importing STL files.
gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
parts = glob('*.stl')                                                      # Obtains a list of names of all STL files in the working directory.
nParts = len(parts)                                                        # Number of model parts.
name = parts[0]
gmsh.model.add(name)                                                       # Starts a new model.
for i in range(nParts):
    gmsh.merge(parts[i])                                                   # Loads the STL model part file located in the working directory.
    s = gmsh.model.getEntities(2)
    ns1 = len(s)

# Preprocessing:
gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
xmin, ymin, zmin, xmax, ymax, zmax = \
gmsh.model.getBoundingBox(-1,-1)                                           # Gets the bounding box of the whole model.
d = sqrt((xmax - xmin) ** 2 + (ymax-ymin) ** 2 + (zmax-zmin) ** 2)         # Calculates the diagonal of the bounding box.
gmsh.logger.write("Model bounding box dimensions are:",'info')
gmsh.logger.write("lx = " + str("{:.4e}".format(xmax - xmin)),'info')
gmsh.logger.write("ly = " + str("{:.4e}".format(ymax - ymin)),'info')
gmsh.logger.write("lz = " + str("{:.4e}".format(zmax - zmin)),'info')
gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
gmsh.logger.write("d  = " + str("{:.4e}".format(d)),'info')
gmsh.logger.write("The value of the key 'GeometryTolerance' was calculate"\
                  "d automaticaly as " + str("{:.4e}".format(d * 1.e-5)), \
                  "warning")
gmsh.option.setNumber('Geometry.Tolerance',d * 1.e-8)                      # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
gmsh.option.setNumber('Geometry.ToleranceBoolean',d * 1.e-8)
# delta = d * 1.e-5
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
s = gmsh.model.getEntities(2)
ns1 = len(s)
gmsh.model.mesh.renumberNodes()
nN = len(gmsh.model.mesh.getNodes(-1,-1,False)[0])
nE = len(gmsh.model.mesh.getElementsByType(2,-1,0)[0])                     # Calculates the total number of triangles in whole STL geometry.
E = gmsh.model.mesh.getElementsByType(2,-1,0)[1]
E = np.reshape(E,(-1,3),'C')
NToE = [[] for i in range(nN)]
EToE = [set() for i in range(nE)]
Ie = -np.ones(nE,np.integer,'C')
Q = deque()
for i in range(nE):                                                        # Loops over every triangle in the STL geometry in order to identify all of its neighbouring nodes.
    for j in range(3):
        NToE[int(E[i,j]) - 1].append(i)
for i in range(nE):                                                        # Loops over every triangle in the STL geometry in order to identify all of its neighbouring triangles.
    for j in range(3):
        EToE[i].update(NToE[int(E[i,j]) - 1])
del NToE
k = 0
for j in range(nE):                                                        # Loops over all triangles in the STL geometry in order to find the total number of surface loops ("islands") in the model.
    if Ie[j] == -1:       
        Q.append(j)
        while Q:
            e = Q.pop()
            if Ie[e] == -1:
                Ie[e] = k
                for i in EToE[e]:
                    if Ie[i] == -1:
                        Q.insert(0,i)
        k += 1
nI = np.max(Ie) + 1
IeOld = Ie
Ie = [[] for i in range(nI)]
for i in range(nE):
    Ie[IeOld[i]].append(i + 1)
del (IeOld,EToE)

solidTags = []


Is = [[] for i in range(nI)]
IV = [[] for i in range(nI)]
nsa = [0,0]
for i in range(ns1):                                                       # Loops over all surfaces in the model geometry in order to link them to the surface loops ("islands") in the model.
    e2 = gmsh.model.mesh.getElementsByType(2,i + 1 + ns1,0)[0][0]
    for j in range(nI):                                                    # Loops over all surface loops ("islands") in the model.
        eInI = any(np.isin(Ie[j],e2))
        if eInI:
            nII = len(IV[j])
            if nII > 0:
                for k in range(nII):                                       # Loops over all volumes linked to those surface loops ("islands").
                    if type(solidTags[i]) is tuple:
                        if solidTags[i][0] == IV[j][k]:
                            Is[j][k].extend([i + 1 + ns1])
                            nsa[0] +=  1
                        if solidTags[i][1] == IV[j][k]:
                            Is[j][k].extend([i + 1 + ns1])
                            nsa[1] +=  1
                    else:
                        if solidTags[i] == IV[j][k]:
                            Is[j][k].extend([i + 1 + ns1])
                            nsa[0] +=  1
                if type(solidTags[i]) is tuple:
                    if nsa[0] == 0:
                        Is[j].append([i + 1 + ns1])
                        IV[j].extend([solidTags[i][0]])
                    if nsa[1] == 0:
                        Is[j].append([i + 1 + ns1])
                        IV[j].extend([solidTags[i][1]])
                else:
                    if nsa[0] == 0:
                        Is[j].append([i + 1 + ns1])
                        IV[j].extend([solidTags[i]])
                nsa = [0,0]
            else:
                if type(solidTags[i]) is tuple:
                    Is[j].append([i + 1 + ns1])
                    Is[j].append([i + 1 + ns1])
                    IV[j].extend([solidTags[i][0]])
                    IV[j].extend([solidTags[i][1]])
                else:
                    Is[j].append([i + 1 + ns1])
                    IV[j].extend([solidTags[i]])

# del Ie
# jj = 1
# for i in range(nSolids):                                                   # Loops over all volumes in the model.
#     ii = [[i1,i2] for i1,v1 in enumerate(IV)                          \
#                   for i2,v2 in enumerate(v1) if v2 == i + 1]               # Gets indicies of all surface loops ("islands") linked to this volume.
#     nii = len(ii)
#     if nii > 1:
#         BBox = np.zeros((nii,3),np.integer,'C')
#         for j in range(nii):                                               # Loops over all linked surface loops ("islands").
#             l = gmsh.model.geo.addSurfaceLoop(Is[ii[j][0]][ii[j][1]])
#             gmsh.model.geo.addVolume([l],jj)
#             gmsh.model.geo.synchronize()
#             xmin, ymin, zmin, xmax, ymax, zmax = \
#             gmsh.model.getBoundingBox(3,jj)
#             BBox[j,:] = [xmax - xmin,ymax - ymin,zmax - zmin]
#             jj += 1
#         iii = np.argmax(BBox[:,0])
# V = gmsh.model.getEntities(3)
# nV = len(V)                                                                # Total number of auxiliary volumes.
# nIV = sum([len(sublist) for sublist in IV])                                # Total number of surface loops ("islands") times auxiliary volumes.
# gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
# if nIV > nSolids:
#     gmsh.model.removeEntities(V,False)

# Writing of individual shells into their respective STL files:
for i in range(ns1):
    gmsh.model.addPhysicalGroup(2, [ns1 + i + 1], i + 1)
    gmsh.write(name + '.V' + str(solidTags[i]) + '.A' + str(i + 1) +      \
               '.stl')
    gmsh.model.removePhysicalGroups()
TPreP = perf_counter() - TPreP

# Write the message console outputs to the Log File:
gmsh.logger.write("Time elapsed for geometry preprocessing: "             \
                  + str("{:.4e}".format(TPreP)) + " s","info")             # Prints the value of time in seconds that it took to preprocess the whole model.
Log = gmsh.logger.get(); gmsh.logger.stop()
writeToLogFile.write(Log,'NeperToGmsh')                                    # The saved file can be located in the working directory.

# Launch the GUI to edit the boundary surface and volume groups:
if '-nopopup' not in argv:
    gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.

# Finalization:
gmsh.finalize()
