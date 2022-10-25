
def IGS(GCF,Log,TModel,TPreP):
    
    # Intitialization:
    from sys import argv
    from os import chdir, rename
    from glob import glob
    from math import sqrt
    from Lib import writeToLogFile
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
    from time import perf_counter
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    
    # Start logging outputs of the message console:
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    
    # Options:
    gmsh.logger.write("All length dimensions in case of IGES format are assum"\
                      "ed to be in milimeters.",'info')
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    # gmsh.option.setNumber('Geometry.OCCTargetUnit','')
    if GCF['GeometryTolerance'][0] is not None:
        gmsh.option.setNumber('Geometry.Tolerance',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
                                                                               # 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',
                          GCF['MeshSizeFromPoints'][0])                        # Compute mesh element sizes from values given at geometry points?
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',
                          GCF['MeshSizeFromCurvature'][0])                     # Automatically compute mesh element sizes from curvature?
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCirclePoints',\
                          GCF['MinimumCirclePoints'][0])                       # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurvePoints',\
                          GCF['MinimumCurvePoints'][0])                        # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumElementsPerTwoPi',\
                          GCF['MinimumElementsPerTwoPi'][0])                   # Sets the minimum number of elements per 2*Pi radians when the mesh size is adapted to the curvature.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',\
                          GCF['SecondOrderLinear'][0])                         # Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',\
                          GCF['Optimize'][0])                                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',\
                          GCF['OptimizeNetgen'][0])                            # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",\
                          GCF['HighOrderOptimize'][0])                         # Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',GCF['MaxNumThreads1D'][0])    # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads2D'][0])    # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.MaxNumThreads3D',GCF['MaxNumThreads3D'][0])    # Sets the maximum number of CPU threads to use for meshing of volumes.
    gmsh.option.setNumber('Mesh.ScalingFactor',GCF['ScalingFactor'][0])        # Global scaling factor applied to the saved mesh.
    gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?
    
    # Geometry import:
    gmsh.option.setNumber('Geometry.OCCSewFaces',1)                            # Sews surfaces into shells in STEP, IGES and BRep geometries.
    gmsh.option.setNumber('Geometry.OCCMakeSolids',1)                          # Fixes shells and make solids in STEP, IGES and BRep geometries.
    gmsh.option.setNumber('Geometry.OCCImportLabels',1)                        # Imports entity labels and colors labels from imported geometry.gmsh.option.setString('Geometry.OCCTargetUnit', 'M').
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.model.add(name)                                                       # Starts a new model.
    parts = glob(name + '*.igs')                                               # Attempts to obtain a list of names of all parts of the model with .igs extension.
    IgesNotIgs = 0
    if parts == []:
        parts = glob(name + '*.iges')                                          # Attempts to obtain a list of names of all parts of the model with .iges extension.
        IgesNotIgs = 1
    if parts == []:
        Log += gmsh.logger.get(); gmsh.logger.stop()
        Log.append("Error: No model geometry files in IGES format found in th"\
                   "e working directory.")                                     # Raises error if no sufficient model geometry files are found in the working directory.
        writeToLogFile.write(Log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    nParts = len(parts)                                                        # Number of model parts.
    solidNames = []
    solidTags = []
    for i in range(nParts):
        part = parts[i].split('.')
        if IgesNotIgs == 1:
            rename(GCF['WorkingDirectoryPath'][0] + '\\' + parts[i],          \
                   GCF['WorkingDirectoryPath'][0] + '\\' + '.'.join(part[:-1])\
                   + '.igs')                                                   # Changes the extension of model geometry files from .iges to .igs (Gmsh doesnt recognize .iges extension).
        if len(part) == 2:
            solidNames.append(part[0])                                         # Extracts names of surface BC groups from model part names.
            solidTags.append((solidNames.index(part[0]) + 1))                  # Assigns tags to those surfaces.
        elif len(part) == 3:
            solidNames.append(part[1])                                         # Extracts names of surface BC groups from model part names.
            solidTags.append((solidNames.index(part[1]) + 1))                  # Assigns tags to those surfaces.
        else:
            Log += gmsh.logger.get(); gmsh.logger.stop()
            Log.append("Error: " + parts[i] + " is not a correct filename for"\
                       " a model part.")                                       # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
            writeToLogFile.write(Log,name)                                     # The saved file can be located in the working directory.
            raise Exception("Fatal error occured, see " + GCF['Name'][0] + "."\
                            "log file for details.")
        gmsh.model.occ.importShapes(parts[i])                                  # Loads IGES geometry file located in the working directory.
    
    # Preprocessing:
    nSolids = len(solidNames)
    gmsh.model.occ.synchronize()                                               # Synchronizes model data in the case of STEP, IGES and BRep geometries.
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1,-1)      # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax-ymin) ** 2 + (zmax-zmin) ** 2)         # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str(xmax - xmin) + " mm,",'info')
    gmsh.logger.write("ly = " + str(ymax - ymin) + " mm,",'info')
    gmsh.logger.write("lz = " + str(zmax - zmin) + " mm.",'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str(d) + " mm.",'info')
    if GCF['GeometryTolerance'][0] is None:
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str(d * 1.e-5) + " mm.","warning")
        gmsh.option.setNumber('Geometry.Tolerance',d * 1.e-5)                  # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str(d / 50) + " mm.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str(d / 50) + " mm.","warning")
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
    for i in range(1,nSolids):
        for j in range(1+k,nSolids+1):
            try:
                VV, _ = gmsh.model.occ.fragment([(3,i)],[(3,j)],-1,True,True)
                if len(VV) > 2:
                    Log += gmsh.logger.get(); gmsh.logger.stop()
                    Log.append("Error: Volumes (3, " + str(i) + ") and (3, " +\
                                str(j) + ") are intersecting. Named volumes s"\
                               "hould be adjacent.")
                    writeToLogFile.write(Log,name)                             # The saved file can be located in the working directory.
                    raise Exception("Fatal error occured, see " +             \
                                     GCF['Name'][0] + ".log file for details.")# Raises error if number of the resultant volumes is larger than two (most probably due to intersection of the volumes).
            except:
                gmsh.logger.write("Boolean fragment on volumes (3, " + str(i) \
                                  + ") and (3, " + str(j) + ") failed. Named "\
                                  "volumes are not adjacent.","warning")       # Raises error if the fragment funcion fails (most probably due to volumes not touching each other).
        k += 1
    gmsh.model.occ.synchronize()                                               # Synchronizes model data in the case of STEP, IGES and BRep geometries.
    s = gmsh.model.getEntities(2)
    V = gmsh.model.getEntities(3)
    
    # Creation of Boundary surface groups for IGES model file:
    shellGroups = {}
    for p in s:
        sName = gmsh.model.getEntityName(p[0],p[1])                            # Gets entity labels read from IGES file.
        if sName:
            gmsh.logger.write("Entity " + str(p) + " has label " + sName + "."\
                              ,"info")                                         # Prints names of all surface entities with successfuly identified labels/names of BC groups.
            sName = sName.split('/')                                           # Extracts the names of BC groups from entity labels.
            if sName[1] not in shellGroups:
                shellGroups[sName[1]] = []
            shellGroups[sName[1]].append(p[1])                                 # Adds names of BC groups to dictionary.
    solidGroups = {}
    for p in V:
        VName = solidNames[p[1] - 1]
        if VName:
            gmsh.logger.write("Entity " + str(p) + " has label " + VName + "."\
                              ,"info")                                         # Prints names of all volume entities with successfuly identified labels/names of BC groups.
            if VName not in solidGroups:
                solidGroups[VName] = []
            solidGroups[VName].append(p[1])
    for shellName, shellTag in shellGroups.items():
        g = gmsh.model.addPhysicalGroup(2,shellTag)                            # Creates boundary surface groups.
        gmsh.model.setPhysicalName(2,g,shellName)                              # Assigns names to the boundary surface groups.
    for solidName, solidTag in solidGroups.items():
        G = gmsh.model.addPhysicalGroup(3,solidTag)                            # Creates volume group.
        gmsh.model.setPhysicalName(3,G,solidName)                              # Assigns name to the volume group.
    
    # Meshing and export:
    gmsh.option.setNumber('Mesh.Format',1)                                     # Sets the mesh output format (1: .msh, 2: .unv, 10: auto).
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    try:
        TMesh = perf_counter()
        gmsh.model.mesh.generate(3)                                            # Generates the finite element mesh.
        gmsh.write(name + '.msh')                                              # The saved file can be located in the working directory.
        TMesh = perf_counter() - TMesh
    except:
        pass
    tModel = perf_counter() - tModel; TModel += tModel
    
    
    # Optionally launch the GUI to see the results:
    if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
        gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.
    
    # Write the message console outputs to the Log File:
    gmsh.logger.write("Time elapsed for preprocessing: " + str(TPreP)         \
                        + " seconds.","info")                                  # Prints the value of time in seconds that it took to preprocess the whole model.
    gmsh.logger.write("Time elapsed for meshing:       " + str(TMesh)         \
                        + " seconds.","info")                                  # Prints the value of time in seconds that it took to mesh the whole model.
    gmsh.logger.write("Total elapsed time:             " + str(TModel)        \
                        + " seconds.","info")                                  # Prints the value of time in seconds that it took to solve the whole task.
    Log += gmsh.logger.get(); gmsh.logger.stop()
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    # Finalization:
    gmsh.finalize()
    return

def RemeshSTLParallel(GCF,parts,i):
    
    # Intitialization:
    from sys import argv
    from os import chdir
    from math import pi, sqrt
    from Lib import writeToLogFile
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
    from time import perf_counter
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    
    # Start logging outputs of the message console:
    tModel = perf_counter(); tPreP = perf_counter()
    Log = []
    gmsh.logger.start()
    
    # Options:
    gmsh.logger.write("All length dimensions in case of STL format are assume"\
                      "d to be in arbitrary units.",'info')
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.<Surface Name>.stl).
    facetAngle = GCF['STLFacetAngle'][0] * pi / 180.                           # Angle between two facets above which an edge is considered as sharp.
    curveAngle = GCF['STLCurveAngle'][0] * pi / 180.                           # Angle between two curve segments above which an edge is considered as sharp.
    if GCF['GeometryTolerance'][0] is not None:
        gmsh.option.setNumber('Geometry.Tolerance',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
                                                                               # 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',
                          GCF['MeshSizeFromPoints'][0])                        # Compute mesh element sizes from values given at geometry points?
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',
                          GCF['MeshSizeFromCurvature'][0])                     # Automatically compute mesh element sizes from curvature?
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCirclePoints',\
                          GCF['MinimumCirclePoints'][0])                       # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurvePoints',\
                          GCF['MinimumCurvePoints'][0])                        # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumElementsPerTwoPi',\
                          GCF['MinimumElementsPerTwoPi'][0])                   # Sets the minimum number of elements per 2*Pi radians when the mesh size is adapted to the curvature.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements? (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',\
                          GCF['SecondOrderLinear'][0])                         # Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',\
                          GCF['Optimize'][0])                                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',\
                          GCF['OptimizeNetgen'][0])                            # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",\
                          GCF['HighOrderOptimize'][0])                         # Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',GCF['MaxNumThreads1D'][0])    # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads2D'][0])    # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.ScalingFactor',GCF['ScalingFactor'][0])        # Global scaling factor applied to the saved mesh.
    gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                # Removes duplicate triangles when importing STL files?
    gmsh.model.add(name)                                                       # Starts a new model.
    
    # Geometry import and preprocessing:
    gmsh.merge(parts[i])                                                       # Loads the STL model part file located in the working directory.
    part = parts[i].split('.')
    if ((len(part) < 2) or (len(part) > 5) or (part[0] != name)):
        Log += gmsh.logger.get(); gmsh.logger.stop()
        Log.append("Error: " + parts[i] + " is not a correct filename for a m"\
                   "odel part.")                                               # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
        writeToLogFile.write(Log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    xmin, ymin, zmin, xmax, ymax, zmax = \
    gmsh.model.getBoundingBox(-1,-1)                                           # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax-ymin) ** 2 + (zmax-zmin) ** 2)         # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model part bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str(xmax - xmin) + ",",'info')
    gmsh.logger.write("ly = " + str(ymax - ymin) + ",",'info')
    gmsh.logger.write("lz = " + str(zmax - zmin) + ".",'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str(d) + ".",'info')
    if GCF['GeometryTolerance'][0] is None:
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str(d * 1.e-5) + ".","warning")
        gmsh.option.setNumber('Geometry.Tolerance',d * 1.e-5)                  # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str(d / 50) + ".","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str(d / 50) + ".","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if one specified by the user is too large.
    if GCF['MeshSizeMin'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMin'. T"\
                          "he value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if none was specified by the user.
    elif GCF['MeshSizeMin'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMin' is too large. N"\
                          "ew value was thus automaticaly set to 0.","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                           # Calculates the default value of the 'MeshSizeMin' parameter if one specified by the user is too large.
    gmsh.model.mesh.classifySurfaces(facetAngle,True,True,curveAngle)          # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
    gmsh.model.mesh.createGeometry()                                           # Creates a geometry for all the discrete curves and surfaces in the mesh.
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.

    # Meshing and export:
    gmsh.option.setNumber('Mesh.Format',27)                                    # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
    tPreP = perf_counter() - tPreP
    try:
        tMesh = perf_counter()
        gmsh.model.mesh.generate(2)                                            # Generates the finite element mesh.
        gmsh.write('New.' + '.'.join(parts[i].split('.')[:-1]) + '.stl')       # The saved file can be located in the working directory.
        tMesh = perf_counter() - tMesh
    except:
        pass
    tModel = perf_counter() - tModel
    
    # Launch the GUI to see the results:
    # if '-nopopup' not in argv:
    #     gmsh.fltk.run()                                                      # Launches the Gmsh GUI to see the result.
    
    # Write the message console outputs to the Log File:
    Log += gmsh.logger.get(); gmsh.logger.stop()
    
    # Finalization:
    gmsh.finalize()
    return Log,tModel,tPreP,tMesh

def RemeshSTL(GCF,Log,TModel,TPreP):
    
    # Intitialization:
    from os import chdir
    from multiprocessing import Pool
    from functools import partial
    from glob import glob 
    from time import perf_counter
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    
    # Options:
    tModel = perf_counter(); tPreP = perf_counter()
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.<Surface Name>.stl).
    P = GCF['STLMaxNumCores'][0]                                               # Sets the maximum number of CPU cores to use for meshing of individual shells during remeshing (STLRemesh = True).
    
    # Geometry import, preprocessing, meshing and export:
    parts = glob(name + '*.stl')                                               # Obtains a list of names of all parts of the model.
    nParts = len(parts)                                                        # Number of model parts.
    I = range(nParts)
    tModel = perf_counter() - tModel; TModel += tModel;
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    with Pool(P) as pool:                                                      # Multiprocessing over individual parts of the model.
        L,t1,t2,t3 = zip(*pool.map(partial(RemeshSTLParallel,GCF,parts),I))    # Geometry import, preprocessing, meshing and export of multiple parts (shells) at the same time.
        for i in L: Log += i;
        TModel += max(t1); TPreP += max(t2); TMesh = max(t3)
    
    # Finalization:
    return Log,TModel,TPreP,TMesh

def RemeshSTLSeriall(GCF,Log,TModel,TPreP):
        
    # Intitialization:
    from sys import argv
    from os import chdir
    from math import pi, sqrt
    from Lib import writeToLogFile
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
    from glob import glob
    from time import perf_counter
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    
    # Start logging outputs of the message console:
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    
    # Options:
    gmsh.logger.write("All length dimensions in case of STL format are assume"\
                      "d to be in arbitrary units.",'info')
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.<Surface Name>.stl).
    facetAngle = GCF['STLFacetAngle'][0] * pi / 180.                           # Angle between two facets above which an edge is considered as sharp.
    curveAngle = GCF['STLCurveAngle'][0] * pi / 180.                           # Angle between two curve segments above which an edge is considered as sharp.
    if GCF['GeometryTolerance'][0] is not None:
        gmsh.option.setNumber('Geometry.Tolerance',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm',GCF['MeshAlgorithm'][0])            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
                                                                               # 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',
                          GCF['MeshSizeFromPoints'][0])                        # Compute mesh element sizes from values given at geometry points?
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',
                          GCF['MeshSizeFromCurvature'][0])                     # Automatically compute mesh element sizes from curvature?
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCirclePoints',\
                          GCF['MinimumCirclePoints'][0])                       # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurvePoints',\
                          GCF['MinimumCurvePoints'][0])                        # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumElementsPerTwoPi',\
                          GCF['MinimumElementsPerTwoPi'][0])                   # Sets the minimum number of elements per 2*Pi radians when the mesh size is adapted to the curvature.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements? (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',\
                          GCF['SecondOrderLinear'][0])                         # Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',\
                          GCF['Optimize'][0])                                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',\
                          GCF['OptimizeNetgen'][0])                            # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",\
                          GCF['HighOrderOptimize'][0])                         # Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('Mesh.MaxNumThreads1D',GCF['MaxNumThreads1D'][0])    # Sets the maximum number of CPU threads to use for meshing of edges.
    gmsh.option.setNumber('Mesh.MaxNumThreads2D',GCF['MaxNumThreads2D'][0])    # Sets the maximum number of CPU threads to use for meshing of surfaces.
    gmsh.option.setNumber('Mesh.ScalingFactor',GCF['ScalingFactor'][0])        # Global scaling factor applied to the saved mesh.
    gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                # Removes duplicate triangles when importing STL files?
    gmsh.model.add(name)                                                       # Starts a new model.
    
    # Geometry import, preprocessing, meshing and export:
    parts = glob(name + '*.stl')                                               # Obtains a list of names of all parts of the model.
    nParts = len(parts)                                                        # Number of model parts.
    tModel = perf_counter() - tModel
    tPreP = perf_counter() - tPreP
    TMesh = 0
    for i in range(nParts):                                                    # Loop over individual parts of the model.
    
        # Geometry import and preprocessing:
        tModel = perf_counter(); tPreP = perf_counter()
        gmsh.merge(parts[i])                                                   # Loads the STL model part file located in the working directory.
        part = parts[i].split('.')
        if ((len(part) < 2) or (len(part) > 5) or (part[0] != name)):
            Log += gmsh.logger.get(); gmsh.logger.stop()
            Log.append("Error: " + parts[i] + " is not a correct filename for"\
                       " a model part.")                                       # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
            writeToLogFile.write(Log,name)                                     # The saved file can be located in the working directory.
            raise Exception("Fatal error occured, see " + GCF['Name'][0] + "."\
                            "log file for details.")
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        xmin, ymin, zmin, xmax, ymax, zmax = \
        gmsh.model.getBoundingBox(-1,-1)                                       # Gets the bounding box of the whole model.
        d = sqrt((xmax - xmin) ** 2 + (ymax-ymin) ** 2 + (zmax-zmin) ** 2)     # Calculates the diagonal of the bounding box.
        gmsh.logger.write("Model part bounding box dimensions are:",'info')
        gmsh.logger.write("lx = " + str(xmax - xmin) + ",",'info')
        gmsh.logger.write("ly = " + str(ymax - ymin) + ",",'info')
        gmsh.logger.write("lz = " + str(zmax - zmin) + ".",'info')
        gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
        gmsh.logger.write("d  = " + str(d) + ".",'info')
        if GCF['GeometryTolerance'][0] is None:
            gmsh.logger.write("No value was provided for the key 'GeometryTol"\
                              "erance'. The value was thus calculated automat"\
                              "icaly as " + str(d * 1.e-5) + ".","warning")
            gmsh.option.setNumber('Geometry.Tolerance',d * 1.e-5)              # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
        if GCF['MeshSizeMax'][0] is None:
            gmsh.logger.write("No value was provided for the key 'MeshSizeMax"\
                              "'. The value was thus calculated automaticaly "\
                              "as " + str(d / 50) + ".","warning")
            gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                  # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
        elif GCF['MeshSizeMax'][0] > d:
            gmsh.logger.write("The value of the key 'MeshSizeMax' is too larg"\
                              "e. New value was thus calculated automaticaly "\
                              "as " + str(d / 50) + ".","warning")
            gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                  # Calculates the default value of the 'MeshSizeMax' parameter if one specified by the user is too large.
        if GCF['MeshSizeMin'][0] is None:
            gmsh.logger.write("No value was provided for the key 'MeshSizeMin"\
                              "'. The value was thus automaticaly set to 0.", \
                              "warning")
            gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                       # Calculates the default value of the 'MeshSizeMin' parameter if none was specified by the user.
        elif GCF['MeshSizeMin'][0] > d:
            gmsh.logger.write("The value of the key 'MeshSizeMin' is too larg"\
                              "e. New value was thus automaticaly set to 0.", \
                              "warning")
            gmsh.option.setNumber('Mesh.MeshSizeMin',0.)                       # Calculates the default value of the 'MeshSizeMin' parameter if one specified by the user is too large.
        gmsh.model.mesh.classifySurfaces(facetAngle,True,True,curveAngle)     # Splits surfaces and their boundaries based on angles between neighbouring facets and curve segments.
        gmsh.model.mesh.createGeometry()                                       # Creates a geometry for all the discrete curves and surfaces in the mesh.
        gmsh.model.geo.synchronize()                                           # Synchronizes model data in the case of STL geometries.
        
        # Meshing and export:
        gmsh.option.setNumber('Mesh.Format',27)                                # Sets the mesh output format (1: .msh, 2: .unv, 10: auto, 27: .stl).
        tPreP = perf_counter() - tPreP; TPreP += tPreP
        try:
            tMesh = perf_counter()
            gmsh.model.mesh.generate(2)                                        # Generates the finite element mesh.
            gmsh.write('New.' + '.'.join(parts[i].split('.')[:-1]) + '.stl')   # The saved file can be located in the working directory.
            tMesh = perf_counter() - tMesh; TMesh += tMesh
        except:
            pass
        tModel = perf_counter() - tModel; TModel += tModel
        
        # Launch the GUI to see the results:
        # if '-nopopup' not in argv:
        #     gmsh.fltk.run()                                                  # Launches the Gmsh GUI to see the result.
        
        # Clear all model data in Gmsh:
        gmsh.clear()                                                           # Clears all the model data.
    
    # Write the message console outputs to the Log File:
    Log += gmsh.logger.get(); gmsh.logger.stop()
    
    # Finalization:
    return Log,TModel,TPreP,TMesh

def STL(GCF,Log,TModel,TPreP,TMesh):
    
    # Intitialization:
    from sys import argv
    from os import chdir
    from glob import glob
    from math import sqrt
    from Lib import writeToLogFile
    from Lib.gmsh.lib import gmsh                                              # Locates the Gmsh source library directory.
    from time import perf_counter
    chdir(GCF['WorkingDirectoryPath'][0])                                      # Locates the working directory.
    gmsh.initialize(argv)
    
    # Start logging outputs of the message console:
    tModel = perf_counter(); tPreP = perf_counter()
    gmsh.logger.start()
    
    # Options:
    name = GCF['Name'][0]                                                      # Model Name (<Model Name>.<Volume Name>.igs).
    if GCF['GeometryTolerance'][0] is not None:
        gmsh.option.setNumber('Geometry.Tolerance',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance.
        gmsh.option.setNumber('Geometry.ToleranceBoolean',\
                              GCF['GeometryTolerance'][0])                     # Geometrical tolerance for boolean operations.
    gmsh.option.setNumber('Mesh.Algorithm3D',GCF['MeshAlgorithm3D'][0])        # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',
                          GCF['MeshSizeFromPoints'][0])                        # Compute mesh element sizes from values given at geometry points?
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',
                          GCF['MeshSizeFromCurvature'][0])                     # Automatically compute mesh element sizes from curvature?
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',
                          GCF['MeshSizeExtendFromBoundary'][0])                # Extend computation of mesh element sizes from the boundaries into the interior?
    if GCF['MeshSizeMax'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax',GCF['MeshSizeMax'][0])        # Sets the maximum global element size.
    if GCF['MeshSizeMin'][0] is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin',GCF['MeshSizeMin'][0])        # Sets the minimum global element size.
    gmsh.option.setNumber('Mesh.MinimumCirclePoints',\
                          GCF['MinimumCirclePoints'][0])                       # Sets the minimum number of nodes used to mesh circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumCurvePoints',\
                          GCF['MinimumCurvePoints'][0])                        # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.
    gmsh.option.setNumber('Mesh.MinimumElementsPerTwoPi',\
                          GCF['MinimumElementsPerTwoPi'][0])                   # Sets the minimum number of elements per 2*Pi radians when the mesh size is adapted to the curvature.
    gmsh.option.setNumber('Mesh.ElementOrder',GCF['ElementOrder'][0])          # Sets the element order.
    gmsh.option.setNumber('Mesh.SecondOrderIncomplete',\
                          GCF['SecondOrderIncomplete'][0])                     # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?
    gmsh.option.setNumber('Mesh.SecondOrderLinear',\
                          GCF['SecondOrderLinear'][0])                         # Should second order nodes, and nodes generated through subdivision algorithms, simply be created by linear interpolation?
    gmsh.option.setNumber('Mesh.Optimize',\
                          GCF['Optimize'][0])                                  # Optimize the mesh to improve the quality of tetrahedral elements?
    gmsh.option.setNumber('Mesh.OptimizeNetgen',\
                          GCF['OptimizeNetgen'][0])                            # Optimize the mesh using Netgen to improve the quality of tetrahedral elements?
    gmsh.option.setNumber("Mesh.HighOrderOptimize",\
                          GCF['HighOrderOptimize'][0])                         # Optimize high-order meshes (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)?
    gmsh.option.setNumber('Mesh.MaxNumThreads3D',GCF['MaxNumThreads3D'][0])    # Sets the maximum number of CPU threads to use for meshing of volumes.
    gmsh.option.setNumber('Mesh.ScalingFactor',GCF['ScalingFactor'][0])        # Global scaling factor applied to the saved mesh.
    gmsh.option.setNumber('Mesh.Binary',GCF['Binary'][0])                      # Write mesh files in binary format (if possible)?
    
    # Geometry import:
    gmsh.option.setNumber('Mesh.StlRemoveDuplicateTriangles',1)                # Removes duplicate triangles when importing STL files.
    gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)                              # Mesh only entities that have no existing mesh.
    gmsh.option.setNumber('Mesh.ColorCarousel',2)                              # Colors the mesh (0: by element type, 1: by elementary entity, 2: by physical group, 3: by mesh partition).
    gmsh.model.add(name)                                                       # Starts a new model.
    parts = glob('New.' + name + '*.stl')                                      # Obtains a list of names of all remeshed parts of the model.
    if len(parts) == 0:
        parts = glob(name + '*.stl')                                           # Obtains a list of names of all original parts of the model.
    nParts = len(parts)                                                        # Number of model parts.
    shellNames = []
    solidNames = []
    solidTags = []
    shellNamesFromFile = []
    for i in range(nParts):
        s0 = gmsh.model.getEntities(2)
        ns0 = len(s0)
        gmsh.merge(parts[i])                                                   # Loads the STL model part file located in the working directory.
        s1 = gmsh.model.getEntities(2)
        ns1 = len(s1)
        part = parts[i].split('.')
        if ((len(part) == 2) and (part[0] == name)):
            if part[0] not in solidNames:
                solidNames.append(part[0] + "_" + str(i + 1))                  # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[0] + "_" + str(ns0 + j + 1))            # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[0] + "_" + str(i + 1))\
                                   + 1))                                       # Assigns tags to those volumes.
                shellNamesFromFile.append(1)
        elif ((len(part) == 3) and (part[0] == 'New') and (part[1] == name)):
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[1] + "_" + str(ns0 + j + 1))            # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.
                shellNamesFromFile.append(1)
        elif ((len(part) == 3) and (part[0] == name)):
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[1] + "_" + str(ns0 + j + 1))            # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.
                shellNamesFromFile.append(1)
        elif ((len(part) == 4) and (part[0] == 'New') and (part[1] == name)):
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[2] + "_" + str(ns0 + j + 1))            # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[2]) + 1))              # Assigns tags to those volumes.
                shellNamesFromFile.append(1)
        elif ((len(part) == 4) and (part[0] == name)):
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[2])                                     # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[1]) + 1))              # Assigns tags to those volumes.
                shellNamesFromFile.append(0)
        elif ((len(part) == 5) and (part[0] == 'New') and (part[1] == name)):
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            for j in range(ns1 - ns0):
                shellNames.append(part[3])                                     # Extracts names of surface BC groups from model part names.
                solidTags.append((solidNames.index(part[2]) + 1))              # Assigns tags to those volumes.
                shellNamesFromFile.append(0)
        elif ((len(part) == 5) and (part[0] == name)):
            if part[1] not in solidNames:
                solidNames.append(part[1])                                     # Extracts names of volume BC groups from model part names.
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            shellNames.append(part[3])                                         # Extracts names of BC groups from model part names.
            solidTags.append((solidNames.index(part[1]) + 1,\
                              solidNames.index(part[2]) + 1))                  # Assigns tags to those volumes.
            shellNamesFromFile.append(0)
        elif ((len(part) == 6) and (part[0] == 'New') and (part[1] == name)):
            if part[2] not in solidNames:
                solidNames.append(part[2])                                     # Extracts names of volume BC groups from model part names.
            if part[3] not in solidNames:
                solidNames.append(part[3])                                     # Extracts names of volume BC groups from model part names.
            shellNames.append(part[4])                                         # Extracts names of BC groups from model part names.
            solidTags.append((solidNames.index(part[2]) + 1,\
                              solidNames.index(part[3]) + 1))                  # Assigns tags to those volumes.
            shellNamesFromFile.append(0)
        else:
            Log += gmsh.logger.get(); gmsh.logger.stop()
            Log.append("Error: " + parts[i] + " is not a correct filename for"\
                       " a model part.")                                       # Raises error if the model part name doesnt have the correct format (<Model Name>.<Volume Name>.<Surface Name>.stl).
            writeToLogFile.write(Log,name)                                     # The saved file can be located in the working directory.
            raise Exception("Fatal error occured, see " + GCF['Name'][0] + "."\
                            "log file for details.")
    
    # Preprocessing:
    nSolids = len(solidNames)
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    xmin, ymin, zmin, xmax, ymax, zmax = \
    gmsh.model.getBoundingBox(-1,-1)                                           # Gets the bounding box of the whole model.
    d = sqrt((xmax - xmin) ** 2 + (ymax-ymin) ** 2 + (zmax-zmin) ** 2)         # Calculates the diagonal of the bounding box.
    gmsh.logger.write("Model bounding box dimensions are:",'info')
    gmsh.logger.write("lx = " + str(xmax - xmin) + ",",'info')
    gmsh.logger.write("ly = " + str(ymax - ymin) + ",",'info')
    gmsh.logger.write("lz = " + str(zmax - zmin) + ".",'info')
    gmsh.logger.write("Diagonal dimension of the bounding box is:",'info')
    gmsh.logger.write("d  = " + str(d) + ".",'info')
    if GCF['GeometryTolerance'][0] is None:
        gmsh.logger.write("No value was provided for the key 'GeometryToleran"\
                          "ce'. The value was thus calculated automaticaly as"\
                          " " + str(d * 1.e-5) + ".","warning")
        gmsh.option.setNumber('Geometry.Tolerance',d * 1.e-5)                  # Calculates the default value of the 'GeometryTolerance' parameter if none was specified by the user.
    if GCF['MeshSizeMax'][0] is None:
        gmsh.logger.write("No value was provided for the key 'MeshSizeMax'. T"\
                          "he value was thus calculated automaticaly as "     \
                           + str(d / 50) + ".","warning")
        gmsh.option.setNumber('Mesh.MeshSizeMax',d / 50.)                      # Calculates the default value of the 'MeshSizeMax' parameter if none was specified by the user.
    elif GCF['MeshSizeMax'][0] > d:
        gmsh.logger.write("The value of the key 'MeshSizeMax' is too large. N"\
                          "ew value was thus calculated automaticaly as "     \
                           + str(d / 50) + ".","warning")
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
                   " part file is continuous.")                                # Raises error if the number of surfaces doesnt match the number of user defined surface BC groups.
        writeToLogFile.write(Log,name)                                         # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    for i in range(nSolids):                                                   # Loop over all volumes in the model.
        ss = []
        for j in range(ns1):                                                   # Loop over all surfaces belonging to these volumes.
            if type(solidTags[j]) is tuple:
                if solidTags[j][0] == i + 1:
                    ss.append(j + 1)
                if solidTags[j][1] == i + 1:
                    ss.append(j + 1)
            else:
                if solidTags[j] == i + 1:
                    ss.append(j + 1)
        l = gmsh.model.geo.addSurfaceLoop(ss)                                  # Creates shell from all surfaces of the volume. WARNING: Creates single watertight shell, currently doesnt support cavities.
        gmsh.model.geo.addVolume([l])                                          # Creates a volume from the shell.
    gmsh.model.geo.synchronize()                                               # Synchronizes model data in the case of STL geometries.
    s = gmsh.model.getEntities(2)
    V = gmsh.model.getEntities(3)
    
    # Creation of Boundary surface groups for STL model file:
    shellGroups = {}
    for p in s:
        if shellNamesFromFile[p[1] - 1] == 1:
            sName = gmsh.model.getEntityName(p[0],p[1])                        # Attempt to get entity labels read from the content of the STL files.
            if len(sName) == 0:
                sName = shellNames[p[1] - 1]                                   # Get entity labels read from the name of the STL files.
        else:
            sName = shellNames[p[1] - 1]                                       # Get entity labels read from the name of the STL files.
        gmsh.logger.write("Entity " + str(p) + " has label " + sName + ".","i"\
                          "nfo")                                               # Prints names of all surface entities with successfuly identified labels/names of BC groups.
        if sName not in shellGroups:
            shellGroups[sName] = []
        shellGroups[sName].append(p[1])                                        # Add names of BC groups to dictionary.
    solidGroups = {}
    for p in V:
        VName = solidNames[p[1] - 1]
        if VName:
            gmsh.logger.write("Entity " + str(p) + " has label " + VName + "."\
                              ,"info")                                         # Prints names of all volume entities with successfuly identified labels/names of BC groups.
            if VName not in solidGroups:
                solidGroups[VName] = []
            solidGroups[VName].append(p[1])
    for shellName, shellTag in shellGroups.items():
        g = gmsh.model.addPhysicalGroup(2,shellTag)                            # Creates boundary surface groups.
        gmsh.model.setPhysicalName(2,g,shellName)                              # Assigns names to the boundary surface groups.
    for solidName, solidTag in solidGroups.items():
        G = gmsh.model.addPhysicalGroup(3,solidTag)                            # Creates volume group.
        gmsh.model.setPhysicalName(3,G,solidName)                              # Assigns name to the volume group.
    
    # Meshing and export:
    gmsh.option.setNumber('Mesh.Format',1)                                     # Sets the mesh output format (1: .msh, 2: .unv, 10: auto).
    tPreP = perf_counter() - tPreP; TPreP += tPreP
    try:
        tMesh = perf_counter()
        gmsh.model.mesh.generate(3)                                            # Generates the finite element mesh.
        gmsh.write(name + '.msh')                                              # The saved file can be located in the working directory.
        tMesh = perf_counter() - tMesh; TMesh += tMesh
    except:
        pass
    tModel = perf_counter() - tModel; TModel += tModel
    
    # Launch the GUI to see the results:
    if ('-nopopup' not in argv) and (GCF['LaunchGmsh'][0] == 1):
        gmsh.fltk.run()                                                        # Launches the Gmsh GUI to see the result.
    
    # Write the message console outputs to the Log File:
    gmsh.logger.write("Time elapsed for preprocessing: " + str(TPreP)         \
                       + " seconds.","info")                                   # Prints the value of time in seconds that it took to preprocess the whole model.
    gmsh.logger.write("Time elapsed for meshing:       " + str(TMesh)         \
                       + " seconds.","info")                                   # Prints the value of time in seconds that it took to mesh the whole model.
    gmsh.logger.write("Total elapsed time:             " + str(TModel)        \
                       + " seconds.","info")                                   # Prints the value of time in seconds that it took to solve the whole task.
    Log += gmsh.logger.get(); gmsh.logger.stop()
    writeToLogFile.write(Log,name)                                             # The saved file can be located in the working directory.
    
    # Finalization:
    gmsh.finalize()
    return
