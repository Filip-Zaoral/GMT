
from sys import argv, platform
from os import chdir,rename,remove,path
from glob import glob
from math import pi,sqrt,floor
from Lib import writeToLogFile
if platform == "win32":
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
else:
    import gmsh                                                                # Locates the Gmsh source library directory.
from time import perf_counter
import numpy as np
from collections import deque
from multiprocessing import Pool
from functools import partial

def IGS(GCF,Log,TModel,TPreP):
    
    # Intitialization:
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    
    # Options:
    gmsh.logger.write("All length dimensions in case of IGES format are assum"\
                      "ed to be in milimeters",'info')
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    if GCF['GeometryTolerance'][0] is not None:
        delta = GCF['GeometryTolerance'][0]
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)               # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',                        # Automatically compute mesh element sizes from curvature,
                          GCF['MeshSizeFromCurvature'][0])                     # using the value as the target number of elements per 2*Pi radians.
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCircleNodes',\
                          GCF['MinimumCircleNodes'][0])                        # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurveNodes',\
                          GCF['MinimumCurveNodes'][0])                         # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    if any(GCF['InflationLayers'][0]) == False:
        gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])      # Sets the element order.
    else:
        gmsh.option.setNumber('Mesh.ElementOrder',1)                           # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',GCF['SecondOrderLinear'][0])# Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',GCF['Optimize'][0])                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen', GCF['OptimizeNetgen'][0])     # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",GCF['HighOrderOptimize'][0])# Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('General.NumThreads',GCF['MaxNumThreads'][0])        # Sets the maximum number of CPU threads to use for 2D/3D meshing.
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',1)                            # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.MaxNumThreads3D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of volumes.
    gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?
    gmsh.option.setNumber('Geometry.OCCSewFaces',1)                            # Sews surfaces into shells in STEP, IGES and BRep geometries.
    gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                          # Fixes shells and make solids in STEP, IGES and BRep geometries.
    gmsh.option.setNumber('Geometry.OCCImportLabels',1)                        # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
    gmsh.option.setNumber('Geometry.OCCParallel',1)                            # Use multi-threaded OpenCASCADE boolean operators.
    gmsh.option.setNumber('Geometry.SurfaceType',2)                            # Surface display type (0: cross, 1: wireframe, 2: solid). Wireframe and solid are not available with the built-in geometry kernel.
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.option.setNumber('Mesh.SurfaceEdges',1)                               # Display edges of surface mesh?
    gmsh.option.setNumber('Mesh.SurfaceFaces',1)                               # Display faces of surface mesh?
    gmsh.option.setNumber('Mesh.VolumeEdges',0)                                # Display edges of volume mesh?
    gmsh.option.setNumber('Mesh.VolumeFaces',0)                                # Display faces of volume mesh?
    gmsh.option.setNumber('Mesh.MeshOnlyEmpty',0)                              # Mesh only entities that have no existing mesh.
    
    # gmsh.option.setNumber('Mesh.RecombinationAlgorithm',2)                     # Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad).
    # gmsh.option.setNumber('Mesh.RecombineAll',True)                            # Apply recombination algorithm to all surfaces, ignoring per-surface spec.
    # gmsh.option.setNumber('Mesh.RecombineOptimizeTopology',30)                 # Number of topological optimization passes (removal of diamonds, ...) of recombined surface meshes.
    # gmsh.option.setNumber('Mesh.Recombine3DAll',False)                         # Apply recombination3D algorithm to all volumes, ignoring per-volume spec (experimental).
    # gmsh.option.setNumber('Mesh.Recombine3DLevel',0)                           # 3D recombination level (0: hex, 1: hex+prisms, 2: hex+prism+pyramids) (experimental).
    # gmsh.option.setNumber('Mesh.Recombine3DConformity',0)                      # 3D recombination level (0: hex, 1: hex+prisms, 2: hex+prism+pyramids) (experimental).
    
    # Geometry import:
    gmsh.model.add(name)                                                       # Starts a new model.
    parts = glob(name + '*.igs') + glob(name + '*.iges'); parts.sort()         # Attempts to obtain a list of names of all parts of the model with .igs extension.
    if parts == []:
        Log += gmsh.logger.get(); gmsh.logger.stop()
        Log.append("Error: No model geometry files in IGES format found in th"\
                   "e working directory")                                      # Raises error if no sufficient model geometry files are found in the working directory.
        writeToLogFile.write(Log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    nParts = len(parts)                                                        # Number of model parts.
    solidNames = []; solidTags = []; shellNames = []
    for i in range(nParts):
        part = parts[i].split('.')
        if part[-1].lower() == "iges":
            rename(path.join(GCF['WorkingDirectoryPath'][0],parts[i]),        \
                   path.join(GCF['WorkingDirectoryPath'][0],                  \
                   '.'.join(part[:-1]) + '.igs'))                              # Changes the extension of model geometry files from .iges to .igs (Gmsh doesnt recognize .iges extension).
        s0 = gmsh.model.getEntities(2); ns0 = len(s0)
        gmsh.model.occ.importShapes(parts[i])                                  # Loads IGES geometry file located in the working directory.
        gmsh.model.occ.synchronize()                                           # Synchronizes model data in the case of STEP, IGES and BRep geometries.
        s1 = gmsh.model.getEntities(2); ns1 = len(s1)
        if len(part) == 2:
            solidNames.append(part[0])                                         # Extracts names of surface BC groups from model part names.
            solidTags.append((solidNames.index(part[0]) + 1))                  # Assigns tags to those surfaces.
            gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                  # Fixes shells and make solids in STEP, IGES and BRep geometries.
        elif len(part) == 3:
            solidNames.append(part[1])                                         # Extracts names of surface BC groups from model part names.
            solidTags.append((solidNames.index(part[1]) + 1))                  # Assigns tags to those surfaces.
            gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                  # Fixes shells and make solids in STEP, IGES and BRep geometries.
        elif len(part) == 4:
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[2])                                     # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.
            gmsh.option.setNumber('Geometry.OCCMakeSolids',0)                  # Fixes shells and make solids in STEP, IGES and BRep geometries.
        elif len(part) == 5:
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            shellNames.append(part[3])                                         # Extracts names of BC groups from model part names.
            solidTags.append((solidNames.index(part[1]) + 1,\
                              solidNames.index(part[2]) + 1))                  # Assigns tags to those volumes.
            gmsh.option.setNumber('Geometry.OCCMakeSolids',0)                  # Fixes shells and make solids in STEP, IGES and BRep geometries.
        else:
            Log += gmsh.logger.get(); gmsh.logger.stop()
            Log.append("Error: " + parts[i] + " is not a correct filename for"\
                       " a model part")                                        # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
            writeToLogFile.write(Log,name)                                     # The saved file can be located in the working directory.
            raise Exception("Fatal error occured, see " + GCF['Name'][0] + "."\
                            "log file for details")
    
    # Preprocessing:
    nV = len(solidNames)
    if (len(part) == 4) or (len(part) == 5):                                   # Form shells into solids in case when geometry input files contain separate surfaces of respective boundary conditions.
        # k = 1
        # for i in range(1,ns1):
        #     for j in range(1+k,ns1+1):
        #         try:
        #             ss, sss = gmsh.model.occ.fragment([(2,i)],[(2,j)],-1,1,1)
        #             gmsh.model.occ.synchronize()
        #             nss = len(ss)
        #             if nss > 2:
        #                 for k in range(nss - 1):
        #                     shellNames.append(shellNames[i - 1])               # Extracts names of surface BC groups from model part names.
        #                     solidTags.append((solidTags[i - 1],               \
        #                                       solidTags[j - 1]))               # Assigns tags to those volumes.
        #             if (nss != 2) or (((2,i) not in ss) or ((2,j) not in ss)):
        #                 shellNames.pop(j - 1)
        #                 solidTags[i - 1] = (solidTags[i - 1],solidTags[j - 1])
        #                 solidTags.pop(j - 1)
        #         except:
        #             pass
        #     k += 1
        s = gmsh.model.getEntities(2); ns = len(s)
        Is = [[] for i in range(nV)]
        for i in range(ns):                                                    # Loop over all surfaces in the model.
            if type(solidTags[i]) is tuple:
                Is[solidTags[i][0] - 1].append(s[i][1])
                Is[solidTags[i][1] - 1].append(s[i][1])
            else:
                Is[solidTags[i] - 1].append(s[i][1])
        for i in range(nV):                                                    # Loop over all volumes in the model.
            gmsh.model.occ.addSurfaceLoop(Is[i],i + 1,True)
            gmsh.model.occ.addVolume([i + 1],i + 1)
            gmsh.model.occ.synchronize()                                       # Synchronizes model data in the case of STEP, IGES and BRep geometries.
    if any(GCF['InflationLayers'][0]) == False:
        P = 1 #P = min(nV,GCF['MaxNumThreads'][0])
    else:
        P = 1
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1,-1)      # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)     # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str("{:.4e}".format(xmax - xmin)) + " mm",    \
                      'info')
    gmsh.logger.write("ly = " + str("{:.4e}".format(ymax - ymin)) + " mm",    \
                      'info')
    gmsh.logger.write("lz = " + str("{:.4e}".format(zmax - zmin)) + " mm",    \
                      'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str("{:.4e}".format(d)) + " mm",'info')
    if GCF['GeometryTolerance'][0] is None:
        delta = d * 1.e-8
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str("{:.4e}".format(delta)) + " mm","warning")
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)) + " mm","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)) + " mm","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if one specified by the user is too large.
    if GCF['MeshSizeMin'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMin'. T"\
                          "he value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if none was specified by the user.
    elif GCF['MeshSizeMin'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMin' is too large. N"\
                          "ew value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if one specified by the user is too large.
    k = 1
    for i in range(1,nV):
        for j in range(1+k,nV+1):
            try:
                VV, _ = gmsh.model.occ.fragment([(3,i)],[(3,j)],-1,1,1)
                gmsh.model.occ.synchronize()
                if len(VV) > 2:
                    Log += gmsh.logger.get(); gmsh.logger.stop()
                    Log.append("Error: Volumes (3, " + str(i) + ") and (3, " +\
                                str(j) + ") are intersecting. Named volumes s"\
                               "hould be adjacent")
                    writeToLogFile.write(Log,name)                             # The saved file can be located in the working directory.
                    raise Exception("Fatal error occured, see " +             \
                                     GCF['Name'][0] + ".log file for details") # Raises error if number of the resultant volumes is larger than two (most probably due to intersection of the volumes).
            except:
                gmsh.logger.write("Boolean fragment on volumes (3, " + str(i) \
                                  + ") and (3, " + str(j) + ") failed. Named "\
                                  "volumes are not adjacent","warning")        # Raises error if the fragment funcion fails (most probably due to volumes not touching each other).
        k += 1
    if (len(part) == 4) or (len(part) == 5):
        gmsh.model.removeEntities(s1,True)
    s = gmsh.model.getEntities(2); ns = len(s)
    V = gmsh.model.getEntities(3); nV = len(V)
    if (len(part) == 2) or (len(part) == 3):
        for i in range(ns):
            sName = gmsh.model.getEntityName(s[i][0],s[i][1])                  # Attempt to get entity labels read from the content of the IGES files.
            sName = sName.split('/')                                           # Extracts the names of BC groups from entity labels.
            if (len(sName) in [0,1]) or ((len(sName) > 1) and                 \
                                         (len(sName[1]) == 0)):
                shellNames.append("A" + str(i + 1))                            # Extracts names of surface BC groups from model part names.
            else:
                nName = shellNames.count(sName[1])
                if nName == 0:
                    shellNames.append(sName[1])
                else:
                    shellNames.append(sName[1] + "_" + str(i + 1))
    if (P > 1) or (any(GCF['InflationLayers'][0]) == True):
        shellTags = [[] for i in range(ns)]
        if nV > 0:
            for i in range(nV):                                                # Forms the list of tags of individual surfaces, linking them to their respective volumes.
                idx = gmsh.model.getBoundary([(3,i + 1)],False,False,False)
                idx = [s.index(j) for j in idx]
                [shellTags[j].extend([i + 1]) for j in idx]
        else:
            [shellTags[j].extend([1]) for j in range(ns)]
        shellTags = [j[0] if len(j) == 1 else tuple(j) for j in shellTags]
    
    # Declaration of Mesh fields:
    nms = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshSurfaces'][0]])
    nmV = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshVolumes'][0]])
    if (len(nms) > 0) and (len(nmV) > 0):
        nmf = sum(nms | nmV)
    else:
        nmf = 0
    if nmf > 0:
        slm = [[] for i in range(nmf)]; slmTags = [[] for i in range(nmf)]
        for i in range(ns):
            for j in range(nmf):
                if shellNames[i] in GCF['LocalMeshSurfaces'][0][j]:
                    slm[j].append((s[i][0],s[i][1]))
                    slmTags[j].extend([s[i][1]])
        Vlm = [[] for i in range(nmf)]; VlmTags = [[] for i in range(nmf)]
        for i in range(nV):
            for j in range(nmf):
                if solidNames[i] in GCF['LocalMeshVolumes'][0][j]:
                    Vlm[j].append((V[i][0],V[i][1]))
                    VlmTags[j].extend([V[i][1]])
        for i in range(0,nmf):
            ii = 4 * (i + 1)
            lms = GCF['LocalMeshSize'][0][i]
            lmg = GCF['LocalMeshGrowthRate'][0][i]
            if Vlm[i]:
                slm[i].extend(gmsh.model.getBoundary(Vlm[i],False,False,False))
                slm[i] = list(set(slm[i]))
                slmTags[i] = [j[1] for j in slm[i]]
            clm = gmsh.model.getBoundary(slm[i],False,False,False)
            clm = list(set(clm)); clmTags = [j[1] for j in clm]
            plm = gmsh.model.getBoundary(clm,False,False,False)
            plm = list(set(plm)); plmTags = [j[1] for j in plm]
            gmsh.model.mesh.field.add("MathEval",ii - 3)
            gmsh.model.mesh.field.setString(ii - 3,"F",str(lms))
            gmsh.model.mesh.field.add("Restrict",ii - 2)
            gmsh.model.mesh.field.setNumbers(ii - 2,"InField",[ii - 3])
            gmsh.model.mesh.field.setNumbers(ii - 2,"VolumesList",VlmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 2,"PointsList",plmTags)
            gmsh.model.mesh.field.add("Distance",ii - 1)
            gmsh.model.mesh.field.setNumbers(ii - 1,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 1,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"PointsList",plmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"NumPointsPerCurve",[1e3])
            gmsh.model.mesh.field.add("MathEval",ii)
            gmsh.model.mesh.field.setString(ii,"F","(1-1/" + str(lmg) + ")*F" \
                                            + str(ii - 1) + "+" + str(lms) +  \
                                            "/" + str(lmg))
        gmsh.model.mesh.field.add("Min",4 * nmf + 1)
        gmsh.model.mesh.field.setNumbers(4 * nmf + 1,"FieldsList",            \
                                         list(range(2,4 * nmf + 1,2)))
        gmsh.model.mesh.field.setAsBackgroundMesh(4 * nmf + 1)
    
    # Meshing:
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    TMesh = perf_counter()
    try:
        if (P > 1) or (any(GCF['InflationLayers'][0]) == True):
            gmsh.model.mesh.generate(GCF['MeshDim'][0] - 1)                    # Generates the finite element mesh.
        else:
            gmsh.model.mesh.generate(GCF['MeshDim'][0])                        # Generates the finite element mesh.
    except:
        pass
    TMesh = perf_counter() - TMesh
    tModel = perf_counter() - tModel; TModel += tModel
    
    # Creation of Boundary surface and volume groups:
    tModel = perf_counter(); tPreP = perf_counter()
    if P > 1:
        nNs = np.zeros(ns,np.integer,'C')
        nNV = np.zeros(nV,np.integer,'C')
        nNP = [[] for i in range(P)]; PP = [[] for i in range(P)]
        for i in range(ns):                                                    # Calculates the number of mesh nodes for each and every surface.
            nNs[i] = gmsh.model.mesh.getNodes(2,s[i][1],False,False)[0].size
        for  i in range(nV):
            ii = [j for j,k in enumerate(shellTags) if                        \
                  ((type(k) is int) and (k == i + 1) or                       \
                  ((type(k) is tuple) and ((k[0] == i + 1) or                 \
                                           (k[1] == i + 1))))]
            nNV[i] = np.sum(nNs[ii])
        nNV[::-1].sort()
        for i in range(nV):                                                    # Calculates the number of volumes per each process with intention to balance the work load per process.
            nNP = [[sum(j)] for j in nNP]
            ii = nNP.index(min(nNP))
            nNP[ii].extend([nNV[i]])
            PP[ii].extend([i + 1])
        PP = list(map(sorted,PP))
    for i in range(P):                                                         # Generates physical groups and saves the mesh into a number P of MSH files.
        if P > 1:
            ii = [j for j,k in enumerate(shellTags) if                        \
                  ((type(k) is int) and (bool(np.isin(k,PP[i]))) or           \
                  ((type(k) is tuple) and ((bool(np.isin(k[0],PP[i],True))) or\
                                           (bool(np.isin(k[1],PP[i],True))))))]
            ss = np.asarray(s)[ii]; ss = [tuple(sublist) for sublist in ss]
        else:
            ss = s
        shellGroups = {}
        for p in ss:
            sName = shellNames[s.index(p)]                                     # Get entity labels read from the name of the STL files.
            gmsh.logger.write("Entity " + str(p) + " has label " + sName,     \
                              "info")                                          # Prints names of all surface entities with successfuly identified labels/names of BC groups.
            if sName not in shellGroups:
                shellGroups[sName] = []
            shellGroups[sName].append(p[1])                                    # Add names of BC groups to dictionary.
        for sName,sTag in shellGroups.items():
            g = gmsh.model.addPhysicalGroup(2,sTag,                           \
                                            shellNames.index(sName) + 1)       # Creates inflation surface groups.
            gmsh.model.setPhysicalName(2,g,sName)                              # Assigns names to the inflation surface groups.
        if (P == 1) and (GCF['MeshDim'][0] == 3):
            solidGroups = {}
            for p in V:
                VName = solidNames[V.index(p)]
                if VName:
                    gmsh.logger.write("Entity " + str(p) + " has label "  \
                                      + VName,"info")                          # Prints names of all volume entities with successfuly identified labels/names of BC groups.
                    if VName not in solidGroups:
                        solidGroups[VName] = []
                    solidGroups[VName].append(p[1])
            for VName,VTag in solidGroups.items():
                G = gmsh.model.addPhysicalGroup(3,VTag,solidNames.index(VName)\
                                                + 1)                           # Creates volume group.
                gmsh.model.setPhysicalName(3,G,VName)                          # Assigns name to the volume group.
        if (P == 1) and (any(GCF['InflationLayers'][0]) == False):
            gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])              # Write mesh files in binary format (if possible)?
            gmsh.option.setNumber("Mesh.SaveAll", 1)                           # Force Gmsh to write also the elements not belonging to any Physical Group.
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
            of = ['.msh','.unv','.vtk']
            gmsh.write(name + of[GCF['OutputFormat'][0] - 1])                  # The saved file can be located in the working directory.
            if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
                tView = perf_counter()
                gmsh.fltk.run()                                                # Launches the Gmsh GUI to see the result.
                tView = perf_counter() - tView
            else:
                tView = 0
        else:
            gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',2)              # Sets the Gmsh to save the shell groups to the output format.
            gmsh.option.setNumber('Mesh.Binary',1)                             # Write mesh files in binary format (if possible)?
            gmsh.option.setNumber("Mesh.SaveAll", 0)                           # Force Gmsh to write also the elements not belonging to any Physical Group.
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
            gmsh.write('Tmp' + str(i + 1) + '.' + name + '.stl')#".msh")       # The saved file can be located in the working directory.
            tView = 0
            # gmsh.fltk.run()                                                  # Launches the Gmsh GUI for debug purposes.
        sNames = list(shellGroups.keys())
        nNames = len(sNames)
        gmsh.model.removePhysicalGroups()
        for j in range(nNames):
            gmsh.model.removePhysicalName(sNames[j])
    tPreP = perf_counter() - tPreP; TPreP += tPreP - tView
    tModel = perf_counter() - tModel; TModel += tModel - tView
    
    # Write the message console outputs to the Log File:
    gmsh.logger.write("Time elapsed for preprocessing: "                      \
                      + str("{:.4e}".format(TPreP)) + " s","info")             # Prints the value of time in seconds that it took to preprocess the whole model.
    gmsh.logger.write("Time elapsed for meshing:       "                      \
                      + str("{:.4e}".format(TMesh)) + " s","info")             # Prints the value of time in seconds that it took to mesh the whole model.
    gmsh.logger.write("Total elapsed time:             "                      \
                      + str("{:.4e}".format(TModel)) + " s","info")            # Prints the value of time in seconds that it took to solve the whole task.
    Log += gmsh.logger.get(); gmsh.logger.stop()
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    # Finalization:
    gmsh.finalize()
    
    # Parallel meshing of volumes:
    if (P > 1) and (GCF['MeshDim'][0] == 3):
        Mesh3DParallel(GCF,Log,solidNames,shellNames,shellTags,TModel,P,PP)
    
    # Generate inflation layers:
    if any(GCF['InflationLayers'][0]) == True:
        Inflation(GCF,Log,TModel,TPreP,TMesh,solidNames,shellNames,shellTags, \
                  d,nmf)
    return

def STL(GCF,Log,TModel,TPreP):
    
    # Intitialization:
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    
    # Options:
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    facetAngle = GCF['STLFacetAngle'][0] * pi / 180.                           # Angle between two facets above which an edge is considered as sharp.
    curveAngle = GCF['STLCurveAngle'][0] * pi / 180.                           # Angle between two curve segments above which an edge is considered as sharp.
    if GCF['GeometryTolerance'][0] is not None:
        delta = GCF['GeometryTolerance'][0]
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)               # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',                        # Automatically compute mesh element sizes from curvature,
                          GCF['MeshSizeFromCurvature'][0])                     # using the value as the target number of elements per 2*Pi radians.
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCircleNodes',\
                          GCF['MinimumCircleNodes'][0])                        # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurveNodes',\
                          GCF['MinimumCurveNodes'][0])                         # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    if any(GCF['InflationLayers'][0]) == False:
        gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])      # Sets the element order.
    else:
        gmsh.option.setNumber('Mesh.ElementOrder',1)                           # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',GCF['SecondOrderLinear'][0])# Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',GCF['Optimize'][0])                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',GCF['OptimizeNetgen'][0])      # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",GCF['HighOrderOptimize'][0])# Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('General.NumThreads',GCF['MaxNumThreads'][0])        # Sets the maximum number of CPU threads to use for 2D/3D meshing.
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',1)                            # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.MaxNumThreads3D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of volumes.
    gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                # Removes duplicate triangles when importing STL files.
    gmsh.option.setNumber('Mesh.MeshOnlyEmpty',0)                              # Mesh only entities that have no existing mesh.
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.option.setNumber('Mesh.SurfaceEdges',1)                               # Display edges of surface mesh?
    gmsh.option.setNumber('Mesh.SurfaceFaces',1)                               # Display faces of surface mesh?
    gmsh.option.setNumber('Mesh.VolumeEdges',0)                                # Display edges of volume mesh?
    gmsh.option.setNumber('Mesh.VolumeFaces',0)                                # Display faces of volume mesh?
    
    # Geometry import:
    gmsh.model.add(name)                                                       # Starts a new model.
    parts = glob(name + '*.stl'); parts.sort()                                 # Obtains a list of names of all remeshed parts of the model.
    nParts = len(parts)                                                        # Number of model parts.
    solidNames = []; shellNames = []; shellTags = []
    for i in range(nParts):
        s0 = gmsh.model.getEntities(2)
        ns0 = len(s0)
        gmsh.merge(parts[i])                                                   # Loads the STL model part file located in the working directory.
        s1 = gmsh.model.getEntities(2)
        ns1 = len(s1)
        part = parts[i].split('.')
        if len(part) == 2:
            solidName = part[0] + "_" + str(i + 1)
            if part[0] not in solidNames:
                solidNames.append(solidName)                                   # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                sName = gmsh.model.getEntityName(s1[ns0 + j][0],s1[ns0 + j][1])# Attempt to get entity labels read from the content of the STL files.
                if len(sName) == 0:
                    shellNames.append("A" + str(ns0 + j + 1))                  # Extracts names of surface BC groups from model part names.
                else:
                    shellNames.append(sName)
                shellTags.append((solidNames.index(solidName) + 1))            # Assigns tags to those volumes.
        elif len(part) == 3:
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                sName = gmsh.model.getEntityName(s1[ns0 + j][0],s1[ns0 + j][1])# Attempt to get entity labels read from the content of the STL files.
                if len(sName) == 0:
                    shellNames.append("A" + str(ns0 + j + 1))                  # Extracts names of surface BC groups from model part names.
                else:
                    shellNames.append(sName)
                shellTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.   
        elif len(part) == 4:
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[2])                                     # Extracts names of surface BC groups from model part names.
                shellTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.
        elif len(part) == 5:
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            shellNames.append(part[3])                                         # Extracts names of BC groups from model part names.
            shellTags.append((solidNames.index(part[1]) + 1,\
                              solidNames.index(part[2]) + 1))                  # Assigns tags to those volumes.
        else:
            Log += gmsh.logger.get(); gmsh.logger.stop()
            Log.append("Error: " + parts[i] + " is not a correct filename for"\
                       " a model part")                                        # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
            writeToLogFile.write(Log,name)                                     # The saved file can be located in the working directory.
            raise Exception("Fatal error occured, see " + GCF['Name'][0] + "."\
                            "log file for details")
    
    # Preprocessing:
    nV = len(solidNames)
    if any(GCF['InflationLayers'][0]) == False:
        P = 1 #P = min(nV,GCF['MaxNumThreads'][0])
    else:
        P = 1
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1,-1)      # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)     # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str("{:.4e}".format(xmax - xmin)),'info')
    gmsh.logger.write("ly = " + str("{:.4e}".format(ymax - ymin)),'info')
    gmsh.logger.write("lz = " + str("{:.4e}".format(zmax - zmin)),'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str("{:.4e}".format(d)),'info')
    if GCF['GeometryTolerance'][0] is None:
        delta = d * 1.e-8
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str("{:.4e}".format(delta)),"warning")
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)),"warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)),"warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if one specified by the user is too large.
    if GCF['MeshSizeMin'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMin'. T"\
                          "he value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if none was specified by the user.
    elif GCF['MeshSizeMin'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMin' is too large. N"\
                          "ew value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if one specified by the user is too large.
    gmsh.model.mesh.removeDuplicateNodes()
    s = gmsh.model.getEntities(2)
    ns0 = len(s)
    gmsh.model.mesh.createTopology()
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    s = gmsh.model.getEntities(2)
    ns1 = len(s)
    if ns0 != ns1:
        Log += gmsh.logger.get(); gmsh.logger.stop()
        Log.append("Error: The Number of continuous shells doesnt match the n"\
                   "umber of model part files.Check if geometry in each model"\
                   " part file is manifold and continuous")                    # Raises error if the number of surfaces doesnt match the number of user defined surface BC groups.
        writeToLogFile.write(Log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if (GCF['STLRemesh'][0]):                                                  # Classifies all surfaces in the model and links the boundary surface groups to the newly created surfaces.
        gmsh.model.mesh.classifySurfaces(facetAngle,True,True,curveAngle)      # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
        gmsh.model.mesh.createGeometry()                                       # Creates a geometry for all the discrete curves and surfaces in the mesh.
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        s = gmsh.model.getEntities(2); ns2 = len(s)
        gmsh.model.add(name + '.Tmp')
        for i in range(nParts):                                                # Loads the original model parts for following linking.
            gmsh.merge(parts[i])
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.createTopology()
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        s = gmsh.model.getEntities(2)
        s2Ins1 = np.zeros((ns2,1),np.bool_,'C')
        sLink = {}
        for i in range(ns1):                                                   # Links the tags of new surfaces to the tags of original surfaces.
            gmsh.model.setCurrent(name + '.Tmp')
            E1 = gmsh.model.mesh.getNodesByElementType(2,i + 1,False)[1]
            E1 = np.reshape(E1,(-1,3,3),'C')
            if len(E1) == 1:
                np.append(E1,np.zeros((1,3,3),np.float64,'C'),0)
            for j in range(ns2):
                if s2Ins1[j] == False:
                    gmsh.model.setCurrent(name)
                    E2 = gmsh.model.mesh.getNodesByElementType(2,j + 1 + ns1, \
                                                               False)[1]
                    E2 = np.reshape(E2,(-1,3,3),'C')
                    E2InE1 = any(np.equal(np.greater_equal(E1,E2[0] - delta), \
                                          np.less_equal(E1,E2[0] + delta))    \
                                 .all(2).all(1))
                    if E2InE1:
                        s2Ins1[j] = True
                        if i + 1 not in sLink:
                            sLink[i + 1] = []
                        sLink[i + 1].append(j + 1 + ns1)
        shellNamesOld = shellNames; shellNames = []
        shellTagsOld = shellTags; shellTags = []
        for i in range(ns1):                                                   # Updates the lists with various tags and names belonging to the surfaces.
            shellNames.extend([shellNamesOld[i]] * len(sLink[i + 1]))
            shellTags.extend([shellTagsOld[i]] * len(sLink[i + 1]))
        gmsh.model.setCurrent(name + '.Tmp')
        gmsh.model.remove()
        gmsh.model.setCurrent(name)
        del (E1,E2)
    else:
        shellNamesOld = shellNames; shellTagsOld = shellTags
        ns2 = ns1; ns1 = 0
    if (P == 1) and (GCF['MeshDim'][0] == 3):
        gmsh.model.mesh.renumberNodes()
        nN = len(gmsh.model.mesh.getNodes(-1,-1,False)[0])
        E = gmsh.model.mesh.getElementsByType(2,-1,0)
        nE = len(E[0]); E = E[1]                                               # Calculates the total number of triangles in whole STL geometry.
        E = np.reshape(E,(-1,3),'C')
        NToE = [[] for i in range(nN)]
        EToE = [set() for i in range(nE)]
        Ie = -np.ones(nE,np.integer,'C')
        Q = deque()
        for i in range(nE):                                                    # Loops over every triangle in the STL geometry in order to identify all of its neighbouring nodes.
            for j in range(3):
                NToE[int(E[i,j]) - 1].append(i)
        for i in range(nE):                                                    # Loops over every triangle in the STL geometry in order to identify all of its neighbouring triangles.
            for j in range(3):
                EToE[i].update(NToE[int(E[i,j]) - 1])
        del (E,NToE)
        k = 0
        for j in range(nE):                                                    # Loops over all triangles in the STL geometry in order to find the total number of surface loops ("islands") in the model.
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
        Is = [[] for i in range(nI)]
        IV = [[] for i in range(nI)]
        nsa = [0,0]
        for i in range(ns2):                                                   # Loops over all surfaces in the model geometry in order to link them to the surface loops ("islands") in the model.
            e2 = gmsh.model.mesh.getElementsByType(2,i + 1 + ns1,0)[0][0]
            for j in range(nI):                                                # Loops over all surface loops ("islands") in the model.
                eInI = any(np.isin(Ie[j],e2))
                if eInI:
                    nII = len(IV[j])
                    if nII > 0:
                        for k in range(nII):                                   # Loops over all volumes linked to those surface loops ("islands").
                            if type(shellTags[i]) is tuple:
                                if shellTags[i][0] == IV[j][k]:
                                    Is[j][k].extend([i + 1 + ns1])
                                    nsa[0] +=  1
                                if shellTags[i][1] == IV[j][k]:
                                    Is[j][k].extend([i + 1 + ns1])
                                    nsa[1] +=  1
                            else:
                                if shellTags[i] == IV[j][k]:
                                    Is[j][k].extend([i + 1 + ns1])
                                    nsa[0] +=  1
                        if type(shellTags[i]) is tuple:
                            if nsa[0] == 0:
                                Is[j].append([i + 1 + ns1])
                                IV[j].extend([shellTags[i][0]])
                            if nsa[1] == 0:
                                Is[j].append([i + 1 + ns1])
                                IV[j].extend([shellTags[i][1]])
                        else:
                            if nsa[0] == 0:
                                Is[j].append([i + 1 + ns1])
                                IV[j].extend([shellTags[i]])
                        nsa = [0,0]
                    else:
                        if type(shellTags[i]) is tuple:
                            Is[j].append([i + 1 + ns1])
                            Is[j].append([i + 1 + ns1])
                            IV[j].extend([shellTags[i][0]])
                            IV[j].extend([shellTags[i][1]])
                        else:
                            Is[j].append([i + 1 + ns1])
                            IV[j].extend([shellTags[i]])
        del Ie
        L = {}
        jj = int(1E9) + 1
        Vx = []
        for i in range(nV):                                                    # Loops over all volumes in the model.
            ii = [[i1,i2] for i1,v1 in enumerate(IV)                          \
                          for i2,v2 in enumerate(v1) if v2 == i + 1]           # Gets indicies of all surface loops ("islands") linked to this volume.
            nii = len(ii)
            if nii > 1:
                BBox = np.zeros((nii,3),np.integer,'C')
                L[i + 1] = []
                for j in range(nii):                                           # Loops over all linked surface loops ("islands").
                    l = gmsh.model.geo.addSurfaceLoop(Is[ii[j][0]][ii[j][1]])
                    L[i + 1].append(l)
                    gmsh.model.geo.addVolume([l],jj)
                    gmsh.model.geo.synchronize()
                    xmin, ymin, zmin, xmax, ymax, zmax = \
                    gmsh.model.getBoundingBox(3,jj)
                    BBox[j,:] = [xmax - xmin,ymax - ymin,zmax - zmin]
                    Vx.append((3,jj))
                    jj += 1
                iii = np.argmax(BBox[:,0])
                if iii != 0:                                                   # Make sure that outermost boundary is on the first place in the surface loop list of the volume.
                    l0 = L[i + 1][iii]
                    L[i + 1].remove(l0); L[i + 1].insert(0,l0)
            else:
                L[i + 1] = []
                l = gmsh.model.geo.addSurfaceLoop(Is[ii[0][0]][ii[0][1]])
                L[i + 1].append(l)
        V = gmsh.model.getEntities(3)
        nIV = sum([len(sublist) for sublist in IV])                            # Total number of surface loops ("islands") times auxiliary volumes.
        for i in range(nV):                                                    # Loops over all volumes in the model.
            gmsh.model.geo.addVolume(L[i + 1],i + 1)                           # Creates a volume from the shell.
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        if nIV > nV:
            gmsh.model.removeEntities(Vx,False)
        s = gmsh.model.getEntities(2); s = s[:ns2]; ns = len(s)
        if P == 1:
            V = gmsh.model.getEntities(3); V = V[-nV:]; nV = len(V)
    
    # Declaration of Mesh fields:
    nms = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshSurfaces'][0]])
    nmV = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshVolumes'][0]])
    if (len(nms) > 0) and (len(nmV) > 0):
        nmf = sum(nms | nmV)
    else:
        nmf = 0
    if nmf > 0:
        slm = [[] for i in range(nmf)]; slmTags = [[] for i in range(nmf)]
        for i in range(ns):
            for j in range(nmf):
                if shellNames[i] in GCF['LocalMeshSurfaces'][0][j]:
                    slm[j].append((s[i][0],s[i][1]))
                    slmTags[j].extend([s[i][1]])
        Vlm = [[] for i in range(nmf)]; VlmTags = [[] for i in range(nmf)]
        if P == 1:
            for i in range(nV):
                for j in range(nmf):
                    if solidNames[i] in GCF['LocalMeshVolumes'][0][j]:
                        Vlm[j].append((V[i][0],V[i][1]))
                        VlmTags[j].extend([V[i][1]])
        for i in range(0,nmf):
            ii = 4 * (i + 1)
            lms = GCF['LocalMeshSize'][0][i]
            lmg = GCF['LocalMeshGrowthRate'][0][i]
            if Vlm[i]:
                slm[i].extend(gmsh.model.getBoundary(Vlm[i],False,False,False))
                slm[i] = list(set(slm[i]))
                slmTags[i] = [j[1] for j in slm[i]]
            clm = gmsh.model.getBoundary(slm[i],False,False,False)
            clm = list(set(clm)); clmTags = [j[1] for j in clm]
            plm = gmsh.model.getBoundary(clm,False,False,False)
            plm = list(set(plm)); plmTags = [j[1] for j in plm]
            gmsh.model.mesh.field.add("MathEval",ii - 3)
            gmsh.model.mesh.field.setString(ii - 3,"F",str(lms))
            gmsh.model.mesh.field.add("Restrict",ii - 2)
            gmsh.model.mesh.field.setNumbers(ii - 2,"InField",[ii - 3])
            gmsh.model.mesh.field.setNumbers(ii - 2,"VolumesList",VlmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 2,"PointsList",plmTags)
            gmsh.model.mesh.field.add("Distance",ii - 1)
            gmsh.model.mesh.field.setNumbers(ii - 1,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 1,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"PointsList",plmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"NumPointsPerCurve",[1e3])
            gmsh.model.mesh.field.add("MathEval",ii)
            gmsh.model.mesh.field.setString(ii,"F","(1-1/" + str(lmg) + ")*F" \
                                            + str(ii - 1) + "+" + str(lms) +  \
                                            "/" + str(lmg))
        gmsh.model.mesh.field.add("Min",4 * nmf + 1)
        gmsh.model.mesh.field.setNumbers(4 * nmf + 1,"FieldsList",            \
                                          list(range(2,4 * nmf + 1,2)))
        gmsh.model.mesh.field.setAsBackgroundMesh(4 * nmf + 1)
    
    # Meshing:
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    TMesh = perf_counter()
    try:
        if (P > 1) or (any(GCF['InflationLayers'][0]) == True):
            gmsh.model.mesh.generate(GCF['MeshDim'][0] - 1)                    # Generates the finite element mesh.
        else:
            gmsh.model.mesh.generate(GCF['MeshDim'][0])                        # Generates the finite element mesh.
    except:
        pass
    TMesh = perf_counter() - TMesh
    tModel = perf_counter() - tModel; TModel += tModel
    
    # Creation of Boundary surface and volume groups and export:
    tPreP = perf_counter()
    if P > 1:
        nNs = np.zeros(ns2,np.integer,'C')
        nNV = np.zeros(nV,np.integer,'C')
        nNP = [[] for i in range(P)]; PP = [[] for i in range(P)]
        for i in range(ns2):                                                   # Calculates the number of mesh nodes for each and every surface.
            nNs[i] = gmsh.model.mesh.getNodes(2,s[i][1],False,False)[0].size
        for  i in range(nV):
            ii = [j for j,k in enumerate(shellTags) if                        \
                  ((type(k) is int) and (k == i + 1) or                       \
                  ((type(k) is tuple) and ((k[0] == i + 1) or                 \
                                           (k[1] == i + 1))))]
            nNV[i] = np.sum(nNs[ii])
        nNV[::-1].sort()
        for i in range(nV):                                                    # Calculates the number of volumes per each process with intention to balance the work load per process.
            nNP = [[sum(j)] for j in nNP]
            ii = nNP.index(min(nNP))
            nNP[ii].extend([nNV[i]])
            PP[ii].extend([i + 1])
        PP = list(map(sorted,PP))
    for i in range(P):                                                         # Generates physical groups and saves the mesh into a number P of MSH files.
        if P > 1:
            ii = [j for j,k in enumerate(shellTags) if                        \
                  ((type(k) is int) and (bool(np.isin(k,PP[i]))) or           \
                  ((type(k) is tuple) and ((bool(np.isin(k[0],PP[i],True))) or\
                                           (bool(np.isin(k[1],PP[i],True))))))]
            ss = np.asarray(s)[ii]; ss = [tuple(sublist) for sublist in ss]
        else:
            ss = s
        shellGroups = {}
        for p in ss:
            sName = shellNames[p[1] - 1 - ns1]                                 # Get entity labels read from the name of the STL files.
            gmsh.logger.write("Entity " + str(p) + " has label " + sName,     \
                              "info")                                          # Prints names of all surface entities with successfuly identified labels/names of BC groups.
            if sName not in shellGroups:
                shellGroups[sName] = []
            shellGroups[sName].append(p[1])                                    # Add names of BC groups to dictionary.
        for sName,sTag in shellGroups.items():
            g = gmsh.model.addPhysicalGroup(2,sTag,                           \
                                            shellNamesOld.index(sName) + 1)    # Creates boundary surface groups.
            gmsh.model.setPhysicalName(2,g,sName)                              # Assigns names to the boundary surface groups.
        if (P == 1) and (GCF['MeshDim'][0] == 3):
            solidGroups = {}
            for p in V:
                VName = solidNames[p[1] - 1]
                if VName:
                    gmsh.logger.write("Entity " + str(p) + " has label "  \
                                      + VName,"info")                          # Prints names of all volume entities with successfuly identified labels/names of BC groups.
                    if VName not in solidGroups:
                        solidGroups[VName] = []
                    solidGroups[VName].append(p[1])
            for VName,VTag in solidGroups.items():
                G = gmsh.model.addPhysicalGroup(3,VTag,solidNames.index(VName)\
                                                + 1)                           # Creates volume group.
                gmsh.model.setPhysicalName(3,G,VName)                          # Assigns name to the volume group.
        if (P == 1) and (any(GCF['InflationLayers'][0]) == False):
            gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])              # Write mesh files in binary format (if possible)?
            gmsh.option.setNumber("Mesh.SaveAll",1)                            # Force Gmsh to write also the elements not belonging to any Physical Group.
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
            of = ['.msh','.unv','.vtk']
            gmsh.write(name + of[GCF['OutputFormat'][0] - 1])                  # The saved file can be located in the working directory.
            if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
                tView = perf_counter()
                gmsh.fltk.run()                                                # Launches the Gmsh GUI to see the result.
                tView = perf_counter() - tView
            else:
                tView = 0
        else:
            gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',2)              # Sets the Gmsh to save the shell groups to the output format.
            gmsh.option.setNumber('Mesh.Binary',1)                             # Write mesh files in binary format (if possible)?
            gmsh.option.setNumber("Mesh.SaveAll",0)                            # Force Gmsh to write also the elements not belonging to any Physical Group.
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
            gmsh.write('Tmp' + str(i + 1) + '.' + name + '.stl')#".msh")       # The saved file can be located in the working directory.
            tView = 0
            # gmsh.fltk.run()                                                  # Launches the Gmsh GUI for debug purposes.
        sNames = list(shellGroups.keys())
        nNames = len(sNames)
        gmsh.model.removePhysicalGroups()
        for j in range(nNames):
            gmsh.model.removePhysicalName(sNames[j])
    tPreP = perf_counter() - tPreP; TPreP += tPreP - tView
    tModel = perf_counter() - tModel; TModel += tModel - tView
    
    # Write the message console outputs to the Log File:
    gmsh.logger.write("Time elapsed for preprocessing: "                      \
                      + str("{:.4e}".format(TPreP)) + " s","info")             # Prints the value of time in seconds that it took to preprocess the whole model.
    gmsh.logger.write("Time elapsed for meshing:       "                      \
                      + str("{:.4e}".format(TMesh)) + " s","info")             # Prints the value of time in seconds that it took to mesh the whole model.
    gmsh.logger.write("Total elapsed time:             "                      \
                      + str("{:.4e}".format(TModel)) + " s","info")            # Prints the value of time in seconds that it took to solve the whole task.
    Log += gmsh.logger.get(); gmsh.logger.stop()
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    # Finalization:
    gmsh.finalize()
    
    # Parallel meshing of volumes:
    if (P > 1) and (GCF['MeshDim'][0] == 3):
        Mesh3DParallel(GCF,Log,solidNames,shellNamesOld,shellTagsOld,TModel,P,\
                       PP)
    
    # Generate inflation layers:
    if any(GCF['InflationLayers'][0]) == True:
        Inflation(GCF,Log,TModel,TPreP,TMesh,solidNames,shellNamesOld,        \
                  shellTagsOld,d,nmf)
    return
    
def Mesh3DParallel(GCF,Log,solidNames,shellNames,shellTags,TModel,P,PP):
    
    # Initialization:
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    
    # Options:
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    
    # Geometry import, preprocessing, meshing and export:
    I = range(P)
    with Pool(P) as pool:                                                      # Multiprocessing over individual parts of the model.
        TPara = perf_counter()
        log,t1,t2,t3 = zip(*pool.map(partial(STL3D,GCF,solidNames,shellNames, \
                                              shellTags,PP),I))                # Geometry import, preprocessing, meshing and export of multiple parts (shells) in parallel.
        TPara = perf_counter() - TPara; TModel += TPara
        for i in log: Log.extend(i)
    # t1 = []; t2 = []; t3 = []
    # for i in range(P):                                                       # Sequential processing over individual parts of the model for debug purposes.
    #     TPara = perf_counter()
    #     log,tt1,tt2,tt3 = STL3D(GCF,solidNames,shellNames,shellTags,PP,i)
    #     t1.extend([tt1]); t2.extend([tt2]); t3.extend([tt3])
    #     for i in log: Log.append(i)
    #     TPara = perf_counter() - TPara; TModel += TPara
    
    # Merging of partial solid meshes to form a singular mesh:
    # TPosP = perf_counter()
    # gmsh.initialize(argv)
    # gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    # gmsh.option.setNumber('Mesh.SurfaceEdges',1)                               # Display edges of surface mesh?
    # gmsh.option.setNumber('Mesh.SurfaceFaces',1)                               # Display faces of surface mesh?
    # gmsh.option.setNumber('Mesh.VolumeEdges',0)                                # Display edges of volume mesh?
    # gmsh.option.setNumber('Mesh.VolumeFaces',0)                                # Display faces of volume mesh?
    # gmsh.option.setNumber("Mesh.SaveAll", 1)                                   # Force Gmsh to write also the elements not belonging to any Physical Group.
    # gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?\
    # gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
    # of = ['.msh','.unv','.vtk']
    # for i in range(P):
    #     gmsh.merge(str(i + 1) + '.' + name + of[GCF['OutputFormat'][0] - 1])
    #     remove(str(i + 1) + '.' + name + of[GCF['OutputFormat'][0] - 1])
    # gmsh.model.mesh.removeDuplicateNodes()
    # gmsh.write(name + of[GCF['OutputFormat'][0] - 1])                          # The saved file can be located in the working directory.
    # TPosP = perf_counter() - TPosP; TModel += TPosP
    
    # Write the message console outputs to the Log File:
    Log.append("Info: Time elapsed in parallel per volume:")
    Log.append("Info: Process No.: Preprocess.: Meshing:     Total:       No."\
               " Volumes:")
    for i in range(P):                                                         # Prints the value of time in seconds that it took to preprocess, mesh and solve the whole model.
        Log.append("Info: " + str("{:12d}".format(i + 1)) + " " +             \
                              str("{:.4e}".format(t1[i])) + " s " +           \
                              str("{:.4e}".format(t2[i])) + " s " +           \
                              str("{:.4e}".format(t3[i])) + " s " +           \
                              str("{:12d}".format(len(PP[i]))))
    Log.append("Info: Total time elapsed for parallel processing: "           \
                + str("{:.4e}".format(TPara)) + " s")                          # Prints the value of time in seconds that it took to preprocess the whole model.
    Log.append("Info: Total elapsed time: " + str("{:.4e}".format(TModel)) +  \
               " s")                                                           # Prints the value of time in seconds that it took to solve the whole task.
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    # Launch the GUI to see the results:
    # if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
    #     gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.
    
    # Finalization:
    # gmsh.finalize()
    return
    
def STL3D(GCF,solidNamesOld,shellNamesOld,shellTags,PP,Pi):
    
    # Intitialization:
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    tModel = perf_counter(); tPreP = perf_counter()
    log = []
    gmsh.logger.start()
    
    # Options:
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    if GCF['GeometryTolerance'][0] is not None:
        delta = GCF['GeometryTolerance'][0]
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)               # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',                        # Automatically compute mesh element sizes from curvature,
                          GCF['MeshSizeFromCurvature'][0])                     # using the value as the target number of elements per 2*Pi radians.
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCircleNodes',\
                          GCF['MinimumCircleNodes'][0])                        # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurveNodes',\
                          GCF['MinimumCurveNodes'][0])                         # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',GCF['SecondOrderLinear'][0])# Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',GCF['Optimize'][0])                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',GCF['OptimizeNetgen'][0])      # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",GCF['HighOrderOptimize'][0])# Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?
    gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                              # Mesh only entities that have no existing mesh.
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.option.setNumber('Mesh.SurfaceEdges',1)                               # Display edges of surface mesh?
    gmsh.option.setNumber('Mesh.SurfaceFaces',1)                               # Display faces of surface mesh?
    gmsh.option.setNumber('Mesh.VolumeEdges',0)                                # Display edges of volume mesh?
    gmsh.option.setNumber('Mesh.VolumeFaces',0)                                # Display faces of volume mesh?
    
    # Geometry import:
    gmsh.model.add(name)                                                       # Starts a new model.
    parts = glob('Tmp' + str(Pi + 1) + '.' + name + '.stl')                    # Obtains a list of names of all remeshed parts of the model.
    shellNames = []
    gmsh.merge(parts[0])                                                       # Loads the STL model part file located in the working directory.
    remove('Tmp' + str(Pi + 1) + '.' + name + '.stl')
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    s1 = gmsh.model.getEntities(2)
    ns1 = len(s1)
    part = parts[0].split('.')
    if ((len(part) == 3) and (part[1] == name)):
        for i in range(ns1):
            sName = gmsh.model.getEntityName(s1[i][0],s1[i][1])                # Attepmts to extract names of surface BC groups from model part names.
            if len(sName) == 0:
                shellNames.append("A" + str(i + 1))                            # Extracts names of surface BC groups from model part names.
            else:
                shellNames.append(sName)
        ii = []
        ii = [i for i, j in enumerate(shellTags) if                           \
              ((type(j) is int)   and (np.isin(j,PP[Pi])))    or              \
              ((type(j) is tuple) and ((np.isin(j[0],PP[Pi])) or              \
                                       (np.isin(j[1],PP[Pi]))))]
        shellTags = [shellTags[idx] for idx in ii]
    else:
        log += gmsh.logger.get(); gmsh.logger.stop()
        log.append("Error: " + parts + " is not a correct filename for a mode"\
                   "l part")                                                   # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
        writeToLogFile.write(log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    
    # Preprocessing:
    nV = len(PP[Pi])
    gmsh.option.setNumber('General.NumThreads',                               \
                          max(1,floor(GCF['MaxNumThreads'][0] / nV + 1)))      # Sets the remaining number of CPU threads to use for 3D meshing in the case of HXT mesher.
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1,-1)      # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)     # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str("{:.4e}".format(xmax - xmin)),'info')
    gmsh.logger.write("ly = " + str("{:.4e}".format(ymax - ymin)),'info')
    gmsh.logger.write("lz = " + str("{:.4e}".format(zmax - zmin)),'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str("{:.4e}".format(d)),'info')
    if GCF['GeometryTolerance'][0] is None:
        delta = d * 1.e-8
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str("{:.4e}".format(delta)),"warning")
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)),"warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str("{:.4e}".format(d / 50)),"warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if one specified by the user is too large.
    if GCF['MeshSizeMin'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMin'. T"\
                          "he value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if none was specified by the user.
    elif GCF['MeshSizeMin'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMin' is too large. N"\
                          "ew value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if one specified by the user is too large.
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.renumberNodes(); gmsh.model.mesh.renumberElements()
    nN = len(gmsh.model.mesh.getNodes(-1,-1,False)[0])
    E = gmsh.model.mesh.getElementsByType(2,-1,0)
    N = E[1]; E = E[0]; nE = len(E)                                            # Calculates the total number of triangles in whole STL geometry.
    N = np.reshape(N,(-1,3),'C')
    gmsh.model.mesh.renumberNodes()
    NToE = [[] for i in range(nN)]
    EToE = [set() for i in range(nE)]
    Ie = -np.ones(nE,np.integer,'C')
    Q = deque()
    for i in range(nE):                                                        # Loops over every triangle in the STL geometry in order to identify all of its neighbouring nodes.
        for j in range(3):
            NToE[int(N[i,j]) - 1].append(i)
    for i in range(nE):                                                        # Loops over every triangle in the STL geometry in order to identify all of its neighbouring triangles.
        for j in range(3):
            EToE[i].update(NToE[int(N[i,j]) - 1])
    del (N,NToE)
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
        Ie[IeOld[i]].append(E[i])
    del IeOld
    Is = [[] for i in range(nI)]
    IV = [[] for i in range(nI)]
    nsa = [0,0]
    for i in range(ns1):                                                       # Loops over all surfaces in the model geometry in order to link them to the surface loops ("islands") in the model.
        e2 = gmsh.model.mesh.getElementsByType(2,s1[i][1],0)[0][0]
        for j in range(nI):                                                    # Loops over all surface loops ("islands") in the model.
            eInI = any(np.isin(Ie[j],e2))
            if eInI:
                nII = len(IV[j])
                if nII > 0:
                    for k in range(nII):                                       # Loops over all volumes linked to those surface loops ("islands").
                        if type(shellTags[s1[i][1] - 1]) is tuple:
                            if shellTags[s1[i][1] - 1][0] == IV[j][k]:
                                Is[j][k].extend([s1[i][1]])
                                nsa[0] +=  1
                            if shellTags[s1[i][1] - 1][1] == IV[j][k]:
                                Is[j][k].extend([s1[i][1]])
                                nsa[1] +=  1
                        else:
                            if shellTags[s1[i][1] - 1] == IV[j][k]:
                                Is[j][k].extend([s1[i][1]])
                                nsa[0] +=  1
                    if type(shellTags[s1[i][1] - 1]) is tuple:
                        if nsa[0] == 0:
                            Is[j].append([s1[i][1]])
                            IV[j].extend([shellTags[s1[i][1] - 1][0]])
                        if nsa[1] == 0:
                            Is[j].append([s1[i][1]])
                            IV[j].extend([shellTags[s1[i][1] - 1][1]])
                    else:
                        if nsa[0] == 0:
                            Is[j].append([s1[i][1]])
                            IV[j].extend([shellTags[s1[i][1] - 1]])
                    nsa = [0,0]
                else:
                    if type(shellTags[s1[i][1] - 1]) is tuple:
                        Is[j].append([s1[i][1]])
                        Is[j].append([s1[i][1]])
                        IV[j].extend([shellTags[s1[i][1] - 1][0]])
                        IV[j].extend([shellTags[s1[i][1] - 1][1]])
                    else:
                        Is[j].append([s1[i][1]])
                        IV[j].extend([shellTags[s1[i][1] - 1]])
    del Ie
    L = {}
    jj = int(1E9) + 1
    Vx = []
    for i in range(nV):                                                        # Loops over all volumes in the model.
        ii = [[i1,i2] for i1,v1 in enumerate(IV)                              \
                      for i2,v2 in enumerate(v1) if v2 == PP[Pi][i]]           # Gets indicies of all surface loops ("islands") linked to this volume.
        nii = len(ii)
        if nii > 1:
            BBox = np.zeros((nii,3),np.integer,'C')
            L[i + 1] = []
            for j in range(nii):                                               # Loops over all linked surface loops ("islands").
                l = gmsh.model.geo.addSurfaceLoop(Is[ii[j][0]][ii[j][1]])
                L[i + 1].append(l)
                gmsh.model.geo.addVolume([l],jj)
                gmsh.model.geo.synchronize()
                xmin, ymin, zmin, xmax, ymax, zmax = \
                gmsh.model.getBoundingBox(3,jj)
                BBox[j,:] = [xmax - xmin,ymax - ymin,zmax - zmin]
                Vx.append((3,jj))
                jj += 1
            iii = np.argmax(BBox[:,0])
            if iii != 0:
                l0 = L[i + 1][iii]
                L[i + 1].remove(l0); L[i + 1].insert(0,l0)
        else:
            L[i + 1] = []
            l = gmsh.model.geo.addSurfaceLoop(Is[ii[0][0]][ii[0][1]])
            L[i + 1].append(l)
    V = gmsh.model.getEntities(3)
    nIV = sum([len(sublist) for sublist in IV])                                # Total number of surface loops ("islands") times auxiliary volumes.
    for i in range(nV):                                                        # Loops over all volumes in the model.
        gmsh.model.geo.addVolume(L[i + 1],PP[Pi][i])                           # Creates a volume from the shell.
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    if nIV > nV:
        gmsh.model.removeEntities(Vx,False)
    s = gmsh.model.getEntities(2)
    V = gmsh.model.getEntities(3)
    
    # Creation of Boundary surface and volume groups for STL model file:
    shellGroups = {}
    for p in s:
        sName = shellNames[p[1] - 1]                                           # Get entity labels read from the name of the STL files.
        gmsh.logger.write("Entity " + str(p) + " has label " + sName,"info")   # Prints names of all surface entities with successfuly identified labels/names of BC groups.
        if sName not in shellGroups:
            shellGroups[sName] = []
        shellGroups[sName].append(p[1]);                                       # Add names of BC groups to dictionary.
    for sName,sTag in shellGroups.items():
        g = gmsh.model.addPhysicalGroup(2,sTag,shellNamesOld.index(sName) + 1) # Creates boundary surface groups.
        gmsh.model.setPhysicalName(2,g,sName)                                  # Assigns names to the boundary surface groups.
    solidGroups = {}
    for p in V:
        VName = solidNamesOld[p[1] - 1]
        gmsh.logger.write("Entity " + str(V[0]) + " has label " + VName,"info")# Prints names of all volume entities with successfuly identified labels/names of BC groups.
        if VName not in solidGroups:
            solidGroups[VName] = []
        solidGroups[VName].append(p[1])
    for VName,VTag in solidGroups.items():
        G = gmsh.model.addPhysicalGroup(3,VTag,solidNamesOld.index(VName) + 1) # Creates volume group.
        gmsh.model.setPhysicalName(3,G,VName)                                  # Assigns name to the volume group.
    
    # Meshing and export:
    tPreP = perf_counter() - tPreP
    tMesh = perf_counter()
    try:
        gmsh.model.mesh.generate(max(GCF['MeshDim'][0],2))                     # Generates the finite element mesh.
    except:
        pass
    tMesh = perf_counter() - tMesh
    gmsh.logger.write("***3D meshing time: " + str(tMesh) + "***","info")
    gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
    of = ['.msh','.unv','.vtk']
    gmsh.write(str(Pi + 1) + "." + name + of[GCF['OutputFormat'][0] - 1])      # The saved file can be located in the working directory.
    tModel = perf_counter() - tModel
    
    # Launch the GUI to see the results:
    if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
        gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.
    
    # Write the message console outputs to the Log File:
    log += gmsh.logger.get(); gmsh.logger.stop()
    
    # Finalization:
    gmsh.finalize()
    return log,tPreP,tMesh,tModel
    
def Inflation(GCF,Log,TModel,TPreP,TMesh,solidNames,shellNames,shellTags,d,   \
              nmf):
    
    # Initialization:
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    name = GCF['Name'][0]                                                      # Model Name.
    nil = GCF['InflationLayers'][0][-1]
    mil = GCF['InflationLayersMethod'][0][-1]
    til = GCF['InflationLayersThickness'][0][-1]
    gil = GCF['InflationLayersGrowthRate'][0][-1]
    if mil != 2:
        hil = til / sum(np.logspace(0,nil - 1,nil,True,gil)); dil = til        # In case that input inflation layer thickness is total layer thickness, this line calculates first layer thickness.
    else:
        dil = til * sum(np.logspace(0,nil - 1,nil,True,gil)); hil = til        # In case that input inflation layer thickness is first layer thickness, this line calculates total layer thickness.
    
    # Options:
    gmsh.model.add(name)
    if GCF['GeometryTolerance'][0] is not None:
        delta = GCF['GeometryTolerance'][0]
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)               # Geometrical tolerance for boolean operations.
    else:
        delta = 1e-8
        gmsh.option.setNumber('Geometry.Tolerance',delta)                      # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',delta)               # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    else:
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    else:
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',GCF['SecondOrderLinear'][0])# Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',GCF['Optimize'][0])                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',GCF['OptimizeNetgen'][0])      # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",GCF['HighOrderOptimize'][0])# Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('General.NumThreads',GCF['MaxNumThreads'][0])        # Sets the maximum number of CPU threads to use for 2D/3D meshing.
    gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                              # Mesh only entities that have no existing mesh.
    gmsh.option.setNumber("Mesh.SaveAll",1)                                    # Force Gmsh to write only the elements belonging to a Physical Group.
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.option.setNumber('Mesh.SurfaceEdges',1)                               # Display edges of surface mesh?
    gmsh.option.setNumber('Mesh.SurfaceFaces',1)                               # Display faces of surface mesh?
    gmsh.option.setNumber('Mesh.VolumeEdges',0)                                # Display edges of volume mesh?
    gmsh.option.setNumber('Mesh.VolumeFaces',0)                                # Display faces of volume mesh?
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',1)                            # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.MaxNumThreads3D',GCF['MaxNumThreads'][0])      # Sets the maximum number of CPU threads to use for meshing of volumes.
    
    # Extract and export the inflation layer surfaces:
    gmsh.model.add(name + ".Tmpb")
    gmsh.merge("Tmp1." + name + ".stl")                                        # Loads the STL model file located in the working directory.
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    gmsh.model.mesh.createTopology()
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    sba = []
    s = gmsh.model.getEntities(2); ns = len(s)
    shellNames = np.asarray(shellNames,np.str)
    shellTags = np.asarray(shellTags,np.object)
    saNames = np.empty((ns,1),np.object,'C')
    nV = len(solidNames)
    for i in s:
        sName = gmsh.model.getEntityName(i[0],i[1])                            # Extracts names of surface BC groups from model part names.
        ii = list(np.where(shellNames == sName)[0])
        if sName in GCF['InflationLayersSurfaces'][0][-1]:
            sba.append((i[0],i[1]))
        else:
            saNames[ii] = sName
    sa = [i for i in s if i not in sba]; nsa = len(sa)
    saNames = saNames[saNames != np.array(None)]
    gmsh.model.removeEntities(sa,True)
    gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                      # Sets the Gmsh to save the shell groups to the output format.
    gmsh.option.setNumber('Mesh.Binary',1)                                     # Write mesh files in binary format (if possible)?
    gmsh.write("Tmpb." + name + ".stl")
    gmsh.model.remove()
    gmsh.model.setCurrent(name)
    
    # Import the inflation layer surfaces:
    gmsh.merge("Tmpb." + name + ".stl")                                        # Loads the STL model file located in the working directory.
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    gmsh.model.mesh.createTopology()
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    sbb = gmsh.model.getEntities(2); nsb = len(sbb)
    if GCF['InputFormat'][0] == 2:
        facetAngle = GCF['STLFacetAngle'][0] * pi / 180.                       # Angle between two facets above which an edge is considered as sharp.
        curveAngle = GCF['STLCurveAngle'][0] * pi / 180.                       # Angle between two curve segments above which an edge is considered as sharp.
        gmsh.model.mesh.classifySurfaces(facetAngle,False,False,curveAngle)    # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
        gmsh.model.mesh.createGeometry()                                       # Creates a geometry for all the discrete curves and surfaces in the mesh.
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        sbb2 = gmsh.model.getEntities(2); nsb2 = len(sbb2)
        gmsh.model.add(name + '.Tmp')
        gmsh.merge("Tmpb." + name + ".stl")
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        gmsh.model.mesh.createTopology()
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        s2Ins1 = np.zeros((nsb2,1),np.bool_,'C')
        sLink = {}
        for i in range(nsb):                                                   # Links the tags of new surfaces to the tags of original surfaces.
            gmsh.model.setCurrent(name + '.Tmp')
            E1 = gmsh.model.mesh.getNodesByElementType(2,i + 1,False)[1]
            E1 = np.reshape(E1,(-1,3,3),'C')
            if len(E1) == 1:
                np.append(E1,np.zeros((1,3,3),np.float64,'C'),0)
            for j in range(nsb2):
                if s2Ins1[j] == False:
                    gmsh.model.setCurrent(name)
                    E2 = gmsh.model.mesh.getNodesByElementType(2,j + 1 + nsb, \
                                                               False)[1]
                    E2 = np.reshape(E2,(-1,3,3),'C')
                    E2InE1 = any(np.equal(np.greater_equal(E1,E2[0] - delta), \
                                          np.less_equal(E1,E2[0] + delta))    \
                                 .all(2).all(1))
                    if E2InE1:
                        s2Ins1[j] = True
                        if i + 1 not in sLink:
                            sLink[i + 1] = []
                        sLink[i + 1].append(j + 1 + nsb)
        surfaceTags = np.zeros((ns),np.int32,'C')
        for i in sbb:
            sName = gmsh.model.getEntityName(i[0],i[1])                        # Extracts names of surface BC groups from model part names.
            ii = list(np.where(shellNames == sName)[0])
            surfaceTags[ii] = i[1]
        shellNamesOld = shellNames.copy()
        shellTagsOld = shellTags.copy()
        surfaceTagsOld = surfaceTags.copy()
        shellNames = [[i] for i in shellNames]
        shellTags = [[i] for i in shellTags]
        surfaceTags = [[i] for i in surfaceTags]
        shellNames = list(shellNames)
        shellTags = list(shellTags)
        surfaceTags = list(surfaceTags)
        sib = [[] for i in range(nV)]
        for i in sbb:
            ii = np.where(surfaceTagsOld == i[1])[0][0]
            shellNames.pop(ii); shellTags.pop(ii); surfaceTags.pop(ii)
            shellNames.insert(ii,[shellNamesOld[ii]] * len(sLink[i[1]]))
            shellTags.insert(ii,[shellTagsOld[ii]] * len(sLink[i[1]]))
            surfaceTags.insert(ii,sLink[i[1]])
            if type(shellTagsOld[ii]) is tuple:
                sib[shellTagsOld[ii][0] - 1].extend([(2,j)                    \
                                                     for j in sLink[i[1]]])
                sib[shellTagsOld[ii][1] - 1].extend([(2,j)                    \
                                                     for j in sLink[i[1]]])
            else:
                sib[shellTagsOld[ii] - 1].extend([(2,j) for j in sLink[i[1]]])
        shellNames = [j for i in shellNames for j in i]
        shellTags = [j for i in shellTags for j in i]
        surfaceTags = [[j] for i in surfaceTags for j in i]
        shellNames = np.asarray(shellNames,np.str)
        shellTags = np.asarray(shellTags,np.object)
        surfaceTags = np.asarray(surfaceTags,np.int)
        gmsh.model.setCurrent(name + '.Tmp')
        gmsh.model.remove()
        gmsh.model.setCurrent(name)
        del (E1,E2)
    else:
        gmsh.model.mesh.classifySurfaces(pi,False,False,pi)                    # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
        gmsh.model.mesh.createGeometry()                                       # Creates a geometry for all the discrete curves and surfaces in the mesh.
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        surfaceTags = np.zeros((ns),np.int32,'C')
        sib = [[] for i in range(nV)]
        for i in sbb:
            sName = gmsh.model.getEntityName(i[0],i[1])                        # Extracts names of surface BC groups from model part names.
            ii = np.where(shellNames == sName)[0][0]
            surfaceTags[ii] = i[1] + nsb
            if type(shellTags[ii]) is tuple:
                sib[shellTags[ii][0] - 1].append((2,i[1] + nsb))
                sib[shellTags[ii][1] - 1].append(i)
            else:
                sib[shellTags[ii] - 1].append((2,i[1] + nsb))
    remove("Tmpb." + name + ".stl")
    sbb = gmsh.model.getEntities(2); nsbb = len(sbb)
    
    # Generate inflation layer volumes:
    gmsh.model.mesh.renumberNodes(); gmsh.model.mesh.renumberElements()
    Nsb = gmsh.model.mesh.getNodes(2,-1,True,False)
    idx = np.unique(Nsb[0],True)[1]                                            # Get indicies of all the nodes on all the boundary layer surfaces.
    Nsb = np.reshape(Nsb[1],(-1,3),'C')[idx,:]                                 # Get sorted coords of all the nodes on all the boundary layer surfaces.
    ENsb = np.reshape(gmsh.model.mesh.getElements(2,-1)[2],(-1,3),'C') - 1     # Get indicies of all the nodes in all the elements on all the boundary layer surfaces.
    nNsb = np.shape(Nsb)[0]; nEsb = np.shape(ENsb)[0]
    M1 = np.transpose(Nsb[ENsb,:],(0,2,1))
    XYZ = np.reshape(M1,(nEsb,9,1),'C')[:,:,0]                                 # Obtain coordinates of all the nodes in all the elements on all the boundary layer surfaces in needed format.
    if mil == 3:
        NToEi = [[] for i in range(nNsb)]; NToEj = [[] for i in range(nNsb)]
        [NToEi[ENsb[i,j]].extend([i]) for j in range(3) for i in range(nEsb)]  # Get element indicies for all the nodes on all the boundary layer surfaces.
        [NToEj[ENsb[i,j]].extend([j]) for j in range(3) for i in range(nEsb)]  # Get local indicies for all the nodes in all the elements on all the boundary layer surfaces.
        vec = np.vectorize(len); NORMALn = vec(np.asarray(NToEi,np.object,'C'))# Calculate to how many elements given node belongs for all the nodes on all the boundary layer surfaces.
        M1 = -np.linalg.norm(M1 - np.roll(M1,1,2),None,1,False) * dil          # Obtain local inflation layer thickness as a fracton of element edge length.
        M1 += np.roll(M1,2,1)                                                  # Distribute the local inflation layer thickness between all the element nodes.
        NORMALNabs = np.asarray([np.sum(M1[NToEi[i],NToEj[i]])                \
                                for i in range(nNsb)],np.float64,'C')          # Summ all the contributions to the local nodal inflation layer thickness from all the respective elements.
        NORMALNabs = NORMALNabs / NORMALn / 2                                  # Calculate the final local nodal inflation layer thickness as a mean value from all the contributions.
        NORMALE = NORMALNabs[ENsb]                                             # Reorganize the local nodal inflation layer thickness values into needed format.
        del (Nsb,ENsb,NORMALNabs,NORMALn,NToEi,NToEj,M1)
        # NToEi = [[] for i in range(nNsb)]; NToEj = [[] for i in range(nNsb)]
        # [NToEi[ENsb[i,j]].extend([i]) for j in range(3) for i in range(nEsb)]  # Get element indicies for all the nodes on all the boundary layer surfaces.
        # [NToEj[ENsb[i,j]].extend([j]) for j in range(3) for i in range(nEsb)]  # Get local indicies for all the nodes in all the elements on all the boundary layer surfaces.
        # vec = np.vectorize(len); NORMALn = vec(np.asarray(NToEi,np.object,'C'))# Calculate to how many elements given node belongs for all the nodes on all the boundary layer surfaces.
        # M1 -= np.roll(M1,1,2)                                                  # Calculate lengths of all edges of all the elements on all the boundary layer surfaces.
        # M2 = np.cross(M1[:,:,0],M1[:,:,1],axis=1)                              # Calculate normal of each element on all the boundary layer surfaces.
        # M2 = M2 / np.linalg.norm(M2,None,1,True)                               # Normalize the normal of each element on all the boundary layer surfaces.
        # M3 = M1 / np.linalg.norm(M1,None,1,True)                               # Calculate unit tangents to each edge of all the elements on all the boundary layer surfaces.
        # M3 = np.roll(np.arccos(np.sum(M3 * np.roll(M3,1,2),1,keepdims=1)),1,2) # Calculate nodal angles for each node of each element on all the boundary layer surfaces.
        # M2 = np.repeat(M2[:,:,None],3,2) * np.repeat(M3,3,1)                   # Weight the unit normals of each element on all the boundary layer surfaces.
        # M1 = -np.linalg.norm(M1,None,1,False) * dil                            # Obtain local inflation layer thickness as a fracton of element edge length.
        # M1 += np.roll(M1,2,1)                                                  # Distribute the local inflation layer thickness between all the element nodes.
        # NORMALNabs = np.asarray([np.sum(M1[NToEi[i],NToEj[i]])                \
        #                         for i in range(nNsb)],np.float64,'C')          # Summ all the contributions to the local nodal inflation layer thickness from all the respective elements.
        # NORMALN = np.asarray([np.sum(M2[NToEi[i],:,NToEj[i]],0)               \
        #                       for i in range(nNsb)],np.float64,'C')             # Summ all the contributions to the unit nodal normal from all the respective elements.
        # NORMALN = NORMALN / np.linalg.norm(NORMALN,None,1,True)                # Normalize the unit nodal normals.
        # NORMALNabs = NORMALNabs / NORMALn / 2                                  # Calculate the final local nodal inflation layer thickness as a mean value from all the contributions.
        # NORMALN = NORMALN * NORMALNabs[:,None]                                 # Scale the unit nodal normals to the value of local inflation layer thickness.
        # NORMALE = np.reshape(NORMALN[ENsb,:],(nEsb,9,1),'C')[:,:,0]            # Reorganize the local nodal inflation layer thickness values into needed format.
        # del (Nsb,ENsb,NToEi,NToEj,NORMALn,M1,M2,NORMALN,NORMALNabs,)
    else:
        NORMALE = np.zeros((nEsb,3),np.float64,'C')
        NORMALE[:,:] = -dil
        del (Nsb,ENsb,M1)
    XYZ = np.reshape(np.concatenate((XYZ,NORMALE),1),(1,-1),'C')[0]            # Combine nodal coordinates and local nodal inflation layer thickness values into needed format.
    blvTag = gmsh.view.add('InflationLayersNormals')
    gmsh.view.addListData(blvTag,'ST',nEsb,XYZ)
    blvi = gmsh.view.getIndex(blvTag)
    Nil = np.linspace(1,1,nil)
    Hil = np.array([hil * sum(np.logspace(0,i - 1,i,True,gil))                \
                    for i in range(1,nil + 1)]) / dil
    Sb = gmsh.model.geo.extrudeBoundaryLayer(sbb,Nil,Hil,True,False,blvi)
    nSb = len(Sb)
    gmsh.model.geo.synchronize()
    VbTags = [Sb[i][1] for i in range(nSb) if Sb[i][0] == 3]; nVb = len(VbTags)
    Sb1 = [Sb[i - 1] for i in range(nSb) if Sb[i][0] == 3]
    Sb = list(set([i for i in Sb if i[0] != 3])); nSb = len(Sb)
    del (NORMALE,XYZ)
    
    # Declaration of Mesh fields:
    idx = np.nonzero(surfaceTags)[0]; s = [(2,i) for i in surfaceTags]
    if nmf > 0:
        slm = [[] for i in range(nmf)]; slmTags = [[] for i in range(nmf)]
        for i in idx:
            for j in range(nmf):
                if shellNames[i] in GCF['LocalMeshSurfaces'][0][j]:
                    slm[j].append((s[i][0],s[i][1]))
                    slmTags[j].extend([s[i][1]])
                if solidNames[shellTags[i] - 1] in GCF['LocalMeshVolumes'][0][j]:
                    slm[j].append((s[i][0],s[i][1]))
                    slmTags[j].extend([s[i][1]])
        for i in range(0,nmf):
            ii = 4 * (i + 1)
            lms = GCF['LocalMeshSize'][0][i]
            lmg = GCF['LocalMeshGrowthRate'][0][i]
            clm = gmsh.model.getBoundary(slm[i],False,False,False)
            clm = list(set(clm)); clmTags = [j[1] for j in clm]
            plm = gmsh.model.getBoundary(clm,False,False,False)
            plm = list(set(plm)); plmTags = [j[1] for j in plm]
            gmsh.model.mesh.field.add("MathEval",ii - 3)
            gmsh.model.mesh.field.setString(ii - 3,"F",str(lms))
            gmsh.model.mesh.field.add("Restrict",ii - 2)
            gmsh.model.mesh.field.setNumbers(ii - 2,"InField",[ii - 3])
            gmsh.model.mesh.field.setNumbers(ii - 2,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 2,"PointsList",plmTags)
            gmsh.model.mesh.field.add("Distance",ii - 1)
            gmsh.model.mesh.field.setNumbers(ii - 1,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 1,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"PointsList",plmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"NumPointsPerCurve",[1e3])
            gmsh.model.mesh.field.add("MathEval",ii)
            gmsh.model.mesh.field.setString(ii,"F","(1-1/" + str(lmg) + ")*F" \
                                            + str(ii - 1) + "+" + str(lms) +  \
                                            "/" + str(lmg))
        gmsh.model.mesh.field.add("Min",4 * nmf + 1)
        gmsh.model.mesh.field.setNumbers(4 * nmf + 1,"FieldsList",            \
                                          list(range(2,4 * nmf + 1,2)))
        gmsh.model.mesh.field.setAsBackgroundMesh(4 * nmf + 1)
    
    # Meshing of the inflation layer volumes:
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    tMesh = perf_counter()
    try:
        gmsh.model.mesh.generate(GCF['MeshDim'][0])                            # Generates the finite element mesh.
    except:
        pass
    gmsh.model.mesh.removeDuplicateNodes()
    tMesh = perf_counter() - tMesh; TMesh += tMesh
    tPreP = perf_counter()
    
    # Extract, modify and export the remaining surfaces:
    Nabxyz = []; Nabd = []; sb01 = []
    sb0 = [[j for j in i if j in sbb] for i in sib]
    if len(saNames) > 0:
        gmsh.model.add(name + ".Tmpa")
        gmsh.merge("Tmp1." + name + ".stl")                                    # Loads the STL model file located in the working directory.
        remove("Tmp1." + name + ".stl")
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.createTopology()
        gmsh.model.geo.synchronize()
        gmsh.model.setCurrent(name)
        cb = [gmsh.model.getBoundary(i,True,False,False) for i in sb0]
        cb = list(set([j for i in cb for j in i]))
        pb = list(set(gmsh.model.getBoundary(cb,False,False,False)))
        pcb = pb + cb; npcb = len(pcb)
        sbaTags = [i[1] for i in sba]; sbbTags = [i[1] for i in sbb]
        Cb1Tags = gmsh.model.getBoundary(Sb1,False,False,False)
        Cb1Tags = list(set([i[1] for i in Cb1Tags]))
        Cb0Tags = gmsh.model.getBoundary(sbb,False,False,False)
        Cb0Tags = list(set([i[1] for i in Cb0Tags]))
        for i in range(npcb):                                                  # Loops over all curves and points of boundary of the inflation layer surfaces with the rest of the model.
            Nb0xyz = gmsh.model.mesh.getNodes(pcb[i][0],pcb[i][1],False,False) # Obtain coordinates and tags of the inflation layer boundary nodes.
            Nb0Tags = Nb0xyz[0]; Nb0xyz = np.reshape(Nb0xyz[1],(-1,3),'C')
            nNb0 = len(Nb0Tags)
            if pcb[i][0] == 1:                                                 # In case of point, obtain tag of the extruded "top" node.
                sbi = gmsh.model.getAdjacencies(1,pcb[i][1])[0]
                sbi = [(2,j) for j in sbi if j not in sbbTags]
                sb01.extend(sbi)
                cb01Tags = gmsh.model.getBoundary(sbi,True,False,False)
                cb01Tags = [j[1] for j in cb01Tags]; cb01Tags.remove(pcb[i][1])
                pcb1Tag = [j for j in cb01Tags if j in Cb1Tags][0]
            else:                                                              # In case of curve, obtain tags of the extruded "top" nodes.
                cb01 = gmsh.model.getAdjacencies(0,pcb[i][1])[0]
                cb01 = [(1,j) for j in cb01 if j not in Cb0Tags]
                pb01Tags = gmsh.model.getBoundary(cb01,True,False,False)
                cb01Tags = [j[1] for j in pb01Tags]; cb01Tags.remove(pcb[i][1])
                pcb1Tag = cb01Tags[0]
            Nb1xyz = gmsh.model.mesh.getNodes(pcb[i][0],pcb1Tag,False,False)[1]# Obtain coordinates of the extruded "top" nodes.
            Nb1xyz = np.reshape(Nb1xyz,(-1,3),'C')
            gmsh.model.setCurrent(name + ".Tmpa")
            for j in range(nNb0):                                              # Loop over all previously selected nodes. 
                Eab = gmsh.model.mesh.getElementsByCoordinates(Nb0xyz[j,0],   \
                                                               Nb0xyz[j,1],   \
                                                               Nb0xyz[j,2],2)  # Obtain tags of all elements one of the nodes of which is located on the given coordinates.
                nEab = len(Eab)
                NaTags = []
                for k in range(nEab):                                          # Loop over all previously selected elements.
                    ea = gmsh.model.mesh.getElement(Eab[k])                    # Obtain surface tag and node tags of the given element.
                    easTag = ea[3]; ea = ea[1]
                    if easTag not in sbaTags:                                  # Select only nodes of those elements, which do not belong to the inflation layer surfaces.
                        NaTags.extend(ea)
                NaTags = list(set(NaTags)); nNa = len(NaTags)                  # Remove duplicities.
                for k in range(nNa):                                           # Loop over previously selected nodes.
                    naxyz = gmsh.model.mesh.getNode(NaTags[k])[0]              # Obtain coordinates of the given node.
                    if (abs(Nb0xyz[j,:] - naxyz) <= delta).all():              # Searches for the node that is located at the same location as the node belonging to the boundary of the inflation layer surfaces.
                        Nbd = Nb1xyz - naxyz
                        Nb1Tona = np.linalg.norm(Nbd,None,1)                   # Calculates distances from the given node to all extruded "top" nodes to find the correct (closest) one.
                        jj = np.argmin(Nb1Tona)
                        Nabxyz.extend(naxyz)
                        Nabd.extend(Nbd[jj,:])
                        gmsh.model.mesh.setNode(NaTags[k],list(Nb1xyz[jj,:]), \
                                                [0.,0.,0.])                    # Translates the node to the location of the extruded "top" node for the mesh to be stitched by removing duplicate nodes later.
            gmsh.model.setCurrent(name)
        gmsh.model.setCurrent(name + ".Tmpa")
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        siab = [[] for i in range(nV)]
        saa = []; caa = []; NaaTags = []; Naaxyz = []
        for i in sba:
            sName = gmsh.model.getEntityName(i[0],i[1])                        # Extracts names of surface BC groups from model part names.
            ii = np.where(shellNames == sName)[0][0]
            if type(shellTags[ii]) is tuple:
                siab[shellTags[ii][0] - 1].append(i)
                siab[shellTags[ii][1] - 1].append(i)
            else:
                siab[shellTags[ii] - 1].append(i)
        cab = [gmsh.model.getBoundary(i,True,False,False) for i in siab]
        cab = list(set([j for i in cab for j in i]))
        pb = list(set(gmsh.model.getBoundary(cab,False,False,False)))
        [saa.extend(list(gmsh.model.getAdjacencies(i[0],i[1])[0]))for i in cab]
        saa = list(set(saa))
        saa = [(2,i) for i in saa if i not in sbaTags]
        if saa:
            # tMorph = perf_counter()
            ca = list(set(gmsh.model.getBoundary(saa,False,False,False)))
            ca = [i for i in ca if i not in cab]; cab = []
            for i in ca:
                nsi = gmsh.model.getAdjacencies(i[0],i[1])[0]
                nsi = [j for j in nsi if (2,j) not in saa]
                cab.append(i) if nsi else caa.append(i)
            cab = list(set(cab)); caa = list(set(caa))
            pab = list(set(gmsh.model.getBoundary(cab,False,False,False)))
            paa = list(set(gmsh.model.getBoundary(caa,False,False,False)))
            pab = [i for i in pab if i not in pb]; pb += pab
            paa = [i for i in paa if i not in pb]
            pcab = pab + cab; pcaa = paa + caa
            for i in pcab:
                Nai = gmsh.model.mesh.getNodes(i[0],i[1],False,False)[1]
                Nabxyz.extend(Nai)
                Nabd.extend(np.zeros((Nai.shape[0]),np.float64,'C'))
            for i in pcaa:
                Nai = gmsh.model.mesh.getNodes(i[0],i[1],False,False)
                NaaTags.extend(Nai[0]); Naaxyz.extend(Nai[1])
            for i in saa:
                Nai = gmsh.model.mesh.getNodes(i[0],i[1],False,False)
                NaaTags.extend(Nai[0]); Naaxyz.extend(Nai[1])
            Nabxyz = np.array(Nabxyz); Nabd = np.array(Nabd)
            NaaTags = np.array(NaaTags); Naaxyz = np.array(Naaxyz)
            Nabxyz = np.reshape(Nabxyz,(-1,3),'C')
            Nabd = np.reshape(Nabd,(-1,3),'C')
            Naaxyz = np.reshape(Naaxyz,(-1,3),'C')
            nNaa = Naaxyz.shape[0]; nNab = Nabxyz.shape[0]
            M1 = np.transpose(np.repeat(Nabxyz[:,:,None],nNab,axis=2),(2,0,1)) # Repeat control point coordinates in 3rd dimension and transpose accordingly.
            M2 = np.transpose(Nabxyz[:,:,None],(0,2,1))                        # Transpose add 3rd dimension and transpose accordingly.
            M1 = np.linalg.norm(M1 - M2,None,2)                                # Substract the matrices from previous steps and calculate Frobenius norm for every element.
            # M1 = 1 / np.sqrt(1 + (0.02 * np.linalg.norm(M1 - M2,None,2)) ** 2)      # Substract the matrices from previous steps and calculate Frobenius norm for every element.
            C = np.concatenate((np.ones((nNab,1),np.float64,'C'),Nabxyz),1)    # Prepend ones to control points cooridnates.
            # C = np.concatenate((C,Nabxyz ** 2),1)
            # C = np.concatenate((C,Nabxyz[:,[0]] * Nabxyz[:,[1]]),1)
            # C = np.concatenate((C,Nabxyz[:,[1]] * Nabxyz[:,[2]]),1)
            # C = np.concatenate((C,Nabxyz[:,[0]] * Nabxyz[:,[2]]),1)
            A = np.block([[M1,C],[C.T,np.zeros((4,4),np.float64,'C')]])        # Assemble left-hand side.
            B = np.concatenate((Nabd,np.zeros((4,3),np.float64,'C')),0)        # Assemble right-hand side.
            # tMorphSolve = perf_counter()
            Q = np.linalg.solve(A,B)                                           # Solve for gammas and betas.
            # tMorphSolve = perf_counter() - tMorphSolve
            GAMMAS = Q[0:nNab,:]; BETAS = Q[nNab:,:]                           # Extract values of gammas and betas.
            M1 = np.repeat(Naaxyz[:,:,None],nNab,axis=2)                       # Repeat controled point coordinates in 3rd dimension.
            M2 = np.transpose(Nabxyz[:,:,None],(2,1,0))                        # Repeat control point coordinates in 3rd dimension and transpose accordingly.
            GAMMAS = np.transpose(GAMMAS[:,:,None],(2,1,0))                    # Repeat gamma values in 3rd dimension and transpose accordingly.
            M1 = np.linalg.norm(M1 - M2,None,1)                                # Calculate Frobenius norm of difference of controled and control point coordinates.
            # M1 = 1 / np.sqrt(1 + (0.02 * np.linalg.norm(M1 - M2,None,1)) ** 2)      # Calculate Frobenius norm of difference of controled and control point coordinates.
            M1 = np.transpose(np.repeat(M1[:,:,None],3,axis=2),(0,2,1))        # Further repeat resultant values from previous step in 3rd dimension and transpose accordingly.
            M1 = np.sum(M1 * GAMMAS,2)                                         # Further multiply resultant values from preivous step with gamma values and sum elements in 3rd dimension.
            M2 = np.concatenate((np.ones((nNaa,1),np.float64,'C'),Naaxyz),1)   # Prepend ones to controled points coordinates.
            # M2 = np.concatenate((M2,Naaxyz ** 2),1)
            # M2 = np.concatenate((M2,Naaxyz[:,[0]] * Naaxyz[:,[1]]),1)
            # M2 = np.concatenate((M2,Naaxyz[:,[1]] * Naaxyz[:,[2]]),1)
            # M2 = np.concatenate((M2,Naaxyz[:,[0]] * Naaxyz[:,[2]]),1)
            M2 = np.transpose(np.repeat(M2[:,:,None],3,axis=2),(0,2,1))        # Further repeat values from previous step in 3rd dimension and transpose accordingly.
            BETAS = np.transpose(BETAS[:,:,None],(2,1,0))                      # Add 3rd dimension to beta values and transpose accordingly.
            M1 += np.sum(M2 * BETAS,2)                                         # Multiply the values from pre-previous step with beta values and sum elements in 3rd dimension.
            Naaxyz += M1                                                       # Add the values from previous step to the controled point coordinates to obtain their morphed coordinates.
            for i in range(nNaa):                                              # Loop over all the controled points.
                gmsh.model.mesh.setNode(NaaTags[i],Naaxyz[i,:],[0.,0.,0.])     # Translate the controll points to their morphed coordinates.
            gmsh.model.geo.synchronize()                                       # Synchronizes model data in the case of STL geometries.
            # tMorph = perf_counter() - tMorph
            del (Nabxyz,Naaxyz,Nabd,M1,M2,A,B,C,Q,GAMMAS,BETAS)
        gmsh.model.removeEntities(sba,True)
        # gmsh.model.mesh.classifySurfaces(pi,False,False,pi)                    # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
        # gmsh.model.mesh.createGeometry()                                       # Creates a geometry for all the discrete curves and surfaces in the mesh.
        # gmsh.model.mesh.optimize('Laplace2D',True,10)
        gmsh.option.setNumber('Mesh.StlOneSolidPerSurface',1)                  # Sets the Gmsh to save the shell groups to the output format.
        gmsh.option.setNumber('Mesh.Binary',1)                                 # Write mesh files in binary format (if possible)?
        gmsh.write("Tmpa." + name + ".stl")
        gmsh.model.remove()
        gmsh.model.setCurrent(name)
        gmsh.merge("Tmpa." + name + ".stl")
        remove("Tmpa." + name + ".stl")
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.removeDuplicateNodes()
        
    # Preprocessing:
    sa = gmsh.model.getEntities(2); sa = sa[nVb + nSb:]; nsa = len(sa)
    for i in range(nsa):
        sName = saNames[i]                                                     # Extracts names of surface BC groups from model part names.
        ii = list(np.where(shellNames == sName)[0])
        surfaceTags[ii] = sa[i][1]
        for j in ii:
            if type(shellTags[j]) is tuple:
                sib[shellTags[j][0] - 1].append(sa[i])
                sib[shellTags[j][1] - 1].append(sa[i])
            else:
                sib[shellTags[j] - 1].append(sa[i])
    shellNames = shellNames.tolist()    
    sib1Tags = [[Sb1[VbTags.index(gmsh.model.getAdjacencies(j[0],j[1])[0][0])]\
                 [1] for j in i] for i in sb0]
    siaTags = [[j[1] for j in i if j not in sbb] for i in sib]
    [siaTags[i].extend(sib1Tags[i]) for i in range(nV)]
    for i in range(nV):
        l = gmsh.model.geo.addSurfaceLoop(siaTags[i])
        gmsh.model.geo.addVolume([l])
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
    jj = 0
    for i in range(nVb,0,-1):
        for j in range(nV):
            try:
                sib1Tags[j].index(Sb1[i - 1][1])
                solidNames = [solidNames[j + jj] + "_I" + str(VbTags[i - 1])] \
                             + solidNames
                jj += 1
            except:
                pass
    s = gmsh.model.getEntities(2); s = s[:nVb] + s[nVb + nSb:]; ns = len(s)
    V = gmsh.model.getEntities(3); nV = len(V)
    
    # Creation of Boundary surface and volume groups:
    shellGroups = {}
    for p in s:
        sName = shellNames[np.where(surfaceTags == p[1])[0][0]]                # Extracts names of surface BC groups from model part names.
        if sName not in shellGroups:
            shellGroups[sName] = []
        shellGroups[sName].append(p[1]);                                       # Add names of BC groups to dictionary.
    for sName,sTag in shellGroups.items():
        g = gmsh.model.addPhysicalGroup(2,sTag,shellNames.index(sName) + 1)    # Creates boundary surface groups.
        gmsh.model.setPhysicalName(2,g,sName)                                  # Assigns names to the boundary surface groups.
    solidGroups = {}
    for p in V:
        VName = solidNames[p[1] - 1]
        if VName not in solidGroups:
            solidGroups[VName] = []
        solidGroups[VName].append(p[1])
    for VName,VTag in solidGroups.items():
        G = gmsh.model.addPhysicalGroup(3,VTag,solidNames.index(VName) + 1)    # Creates volume group.
        gmsh.model.setPhysicalName(3,G,VName)                                  # Assigns name to the volume group.
    
    # Declaration of Mesh fields:
    s = [(2,i) for i in surfaceTags]
    if nmf > 0:
        slm = [[] for i in range(nmf)]; slmTags = [[] for i in range(nmf)]
        for i in range(ns):
            for j in range(nmf):
                if shellNames[i] in GCF['LocalMeshSurfaces'][0][j]:
                    slm[j].append((s[i][0],s[i][1]))
                    slmTags[j].extend([s[i][1]])
        Vlm = [[] for i in range(nmf)]; VlmTags = [[] for i in range(nmf)]
        for i in range(nV):
            for j in range(nmf):
                if solidNames[i] in GCF['LocalMeshVolumes'][0][j]:
                    Vlm[j].append((V[i][0],V[i][1]))
                    VlmTags[j].extend([V[i][1]])
        for i in range(0,nmf):
            ii = 4 * (i + nmf + 1) + 1
            lms = GCF['LocalMeshSize'][0][i]
            lmg = GCF['LocalMeshGrowthRate'][0][i]
            if Vlm[i]:
                slm[i].extend(gmsh.model.getBoundary(Vlm[i],False,False,False))
                slm[i] = list(set(slm[i]))
                slmTags[i] = [j[1] for j in slm[i]]
            clm = gmsh.model.getBoundary(slm[i],False,False,False)
            clm = list(set(clm)); clmTags = [j[1] for j in clm]
            plm = gmsh.model.getBoundary(clm,False,False,False)
            plm = list(set(plm)); plmTags = [j[1] for j in plm]
            gmsh.model.mesh.field.add("MathEval",ii - 3)
            gmsh.model.mesh.field.setString(ii - 3,"F",str(lms))
            gmsh.model.mesh.field.add("Restrict",ii - 2)
            gmsh.model.mesh.field.setNumbers(ii - 2,"InField",[ii - 3])
            gmsh.model.mesh.field.setNumbers(ii - 2,"VolumesList",VlmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 2,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 2,"PointsList",plmTags)
            gmsh.model.mesh.field.add("Distance",ii - 1)
            gmsh.model.mesh.field.setNumbers(ii - 1,"SurfacesList",slmTags[i])
            gmsh.model.mesh.field.setNumbers(ii - 1,"CurvesList",clmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"PointsList",plmTags)
            gmsh.model.mesh.field.setNumbers(ii - 1,"NumPointsPerCurve",[1e3])
            gmsh.model.mesh.field.add("MathEval",ii)
            gmsh.model.mesh.field.setString(ii,"F","(1-1/" + str(lmg) + ")*F" \
                                            + str(ii - 1) + "+" + str(lms) +  \
                                            "/" + str(lmg))
        gmsh.model.mesh.field.add("Min",8 * nmf + 2)
        gmsh.model.mesh.field.setNumbers(8 * nmf + 2,"FieldsList",            \
                                         list(range(4 * nmf + 3,8 * nmf + 2,2)))
        gmsh.model.mesh.field.setAsBackgroundMesh(8 * nmf + 2)
    
    # Meshing of the remaining volumes:
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    tMesh = perf_counter()
    try:
        gmsh.model.mesh.generate(GCF['MeshDim'][0])                            # Generates the finite element mesh.
    except:
        pass
    tMesh = perf_counter() - tMesh; TMesh += tMesh
    gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?
    gmsh.option.setNumber("Mesh.SaveAll", 1)                                   # Force Gmsh to write also the elements not belonging to any Physical Group.
    gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
    of = ['.msh','.unv','.vtk']
    gmsh.write(name + of[GCF['OutputFormat'][0] - 1])                          # The saved file can be located in the working directory.
    tModel = perf_counter() - tModel; TModel += tModel
    
    # Write the message console outputs to the Log File:
    gmsh.logger.write("Time elapsed for preprocessing: "                      \
                      + str("{:.4e}".format(TPreP)) + " s","info")             # Prints the value of time in seconds that it took to preprocess the whole model.
    gmsh.logger.write("Time elapsed for meshing:       "                      \
                        + str("{:.4e}".format(TMesh)) + " s","info")           # Prints the value of time in seconds that it took to mesh the whole model.
    gmsh.logger.write("Total elapsed time:             "                      \
                      + str("{:.4e}".format(TModel)) + " s","info")            # Prints the value of time in seconds that it took to solve the whole task.
    Log += gmsh.logger.get(); gmsh.logger.stop()
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
        gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.
    
    # Finalization:
    gmsh.finalize()
    return