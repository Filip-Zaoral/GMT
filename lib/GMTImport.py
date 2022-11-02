
def Import():
    
    # Intitialization:
    from glob import glob
    from sys import platform
    from os import getenv, path
    from distutils.util import strtobool
    from lib import GMTMesh, writeToLogFile
    from time import perf_counter
    import re
    import numpy as np
    
    # Declaration of default GMT configuration options:
    Log = []
    GCF = {
           'WorkingDirectoryPath' : [None,3],                                  # Absolute path to the working directory.                                                               [string]
           'Name' : [None,3],                                                  # Model Name (IGS: <Model Name>.<Volume Name>.igs, STL: <Model Name>.<Volume Name>.<Surface Name>.stl). [string]
           'InputFormat' : [1,1],                                              # Format of the input model geometry files (1: IGES, 2: STL).                                           [integer]
           'STLRemesh' : [False,0],                                            # Parametrize and remesh geometry of the model when using STL file format (Format = 2)?                 [boolean]
           'STLFacetAngle' : [45.,2],                                          # Angle between two facets above which an edge is considered as sharp.                                  [float]
           'STLCurveAngle' : [180.,2],                                         # Angle between two curve segments above which an edge is considered as sharp.                          [float]
           'GeometryTolerance' : [None,2],                                     # Geometrical tolerance.                                                                                [float]
           'MeshDim' : [3,1],                                                  # Sets the number of dimensions of the finite element mesh (1: 1D, 2: 2D, 3: 3D).                       [integer]
           'MeshAlgorithm' : [6,1],                                            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
                                                                               # 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).           [integer]
           'MeshAlgorithm3D' : [10,1],                                         # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).      [integer]
           'MeshSizeFromCurvature' : [6,1],                                    # Automatically compute mesh element sizes from curvature,
                                    										   # using the value as the target number of elements per 2*Pi radians.                                    [integer]
           'MeshSizeExtendFromBoundary' : [True,0],                            # Extend computation of mesh element sizes from the boundaries into the interior?                       [boolean]
           'MeshSizeMax' : [None,2],                                           # Sets the maximum global element size.                                                                 [float]
           'MeshSizeMin' : [None,2],                                           # Sets the minimum global element size.                                                                 [float]
           'MinimumCircleNodes' : [0,1],                                       # Sets the minimum number of nodes used to mesh circles and ellipses.                                   [integer]
           'MinimumCurveNodes' : [0,1],                                        # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.         [integer]
           'ElementOrder' : [1,1],                                             # Sets the element order.                                                                               [integer]
           'SecondOrderIncomplete' : [True,0],                                 # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?                          [boolean]
           'SecondOrderLinear' : [False,0],                                    # Should second order nodes, and nodes generated through subdivision
                                    										   # algorithms, simply be created by linear interpolation?                                                [boolean]
           'Optimize' : [True,0],                                              # Optimize the mesh to improve the quality of tetrahedral elements?                                     [boolean]
           'OptimizeNetgen' : [False,0],                                       # Optimize the mesh using Netgen library to improve the quality of tetrahedral elements?                [boolean]
           'HighOrderOptimize' : [4,1],                                        # High-order mesh optimization algorithm (0: none, 1: optimization,
                                                                               # 2: elastic+optimization, 3: elastic, 4: fast curving).                                                [integer]
           'MaxNumThreads' : [1,1],                                            # Sets the maximum number of CPU cores to use for meshing of individual solids.                         [integer]
           'OutputFormat' : [1,1],                                             # Format of the output model mesh files (1: MSH, 2: UNV, 3 (16): VTK).                                  [integer]
           'Binary' : [False,0],                                               # Write mesh files in binary format (if possible)?                                                      [boolean]
           'LaunchGmsh' : [False,0],                                           # Launch Gmsh to see the result after finishing th task?                                                [boolean]
           'InflationLayersSurfaces' :  [[],4],                                # List of surfaces on which to generate inflation layers.                                               [list]
           'InflationLayers' : [[],1],            					           # Number of inflation layers to generate on the selected surfaces.                                      [integer]
           'InflationLayersMethod' : [[],1],                                   # Method of specification of inflation layers (1: Total thickess, 2: First layer thickness, 
                                                                               # 3: Total aspect ratio, 4: First layer aspect ratio, 5: Last layer transition ratio).                  [integer]
           'InflationLayersThickness' : [[],2],                                # Inflation layer thickness parameter according to specified method.                                    [float]
           'InflationLayersGrowthRate' : [[],2],                               # Rate of change of thickness of two neighbouring inflation layers, belonging to the selected surface.  [float]
           'LocalMeshSurfaces' : [[],4],                                       # List of surfaces on whitch to enforce local element size.                                             [list]
           'LocalMeshVolumes' : [[],4],                                        # List of volumes on whitch to enforce local element size.                                              [list]
           'LocalMeshSize' : [[],2],                                           # Target local size of the finite elements on the selected surfaces/volumes.                            [float]
           'LocalMeshGrowthRate' : [[],2]                                      # Rate of change of size of the neighbouring elements near the selected sufraces/volumes.               [float]
           }
    
    # Import of configuration options from the GMT Configuration File:
    TModel = perf_counter(); TPreP = perf_counter()
    SignedParameters = ['inflationlayersgrowthrate','localmeshgrowthrate']
    LocalParameters = ['inflationlayers','inflationlayersmethod',             \
                       'inflationlayersthickness','inflationlayersgrowthrate',\
                       'inflationlayerssurfaces','localmeshsize',             \
                       'localmeshgrowthrate','localmeshsurfaces',             \
                       'localmeshvolumes']
    LocalParametersSV = ['localmeshsurfaces','localmeshvolumes',              \
                         'inflationlayerssurfaces']
    GCFName = glob(path.join(getenv('GMTPATH'),                               \
                             'GMT_Configuration_File.gcf'))[0]                 # Obtains a list of names of all GMT Configuration Files in the installation directory.
    GCFRaw = open(GCFName,'r')                                                 # Opens the GMT Configuration File for reading.
    GCFLines = GCFRaw.readlines()                                              # Reads the GMT Configuration File.
    GCFRaw.close()                                                             # Closes the GMT Configuration File.
    nGCFLines = len(GCFLines)
    for i in range(nGCFLines):
        GCFLines[i] = GCFLines[i].replace(" ","").replace("\t","")
    for Key, Value in GCF.items():                                             # Loops over key - value pairs in the dictionary of default GMT configuration options.
        for i in range(nGCFLines):                                             # Loops over key - value pairs in the user specified GMT configuration file.
            GCFLine = GCFLines[i].splitlines()[0]
            idxA = GCFLine.find(Key)                                           # Looks for keys in the GMT configuration file.
            if idxA == 0:                                                      # Checks for lines with keys.
                idxA = GCFLine.find('=') + 1                                   # Looks for the begining index of the value entry in the line.
                idxB = GCFLine.find('#')                                       # Looks for the endg index of the value entry in the line.
                Label = GCFLine[:idxA - 1].lower()
                if idxB != -1:
                    Parameter = GCFLine[idxA:idxB]                             # Extracts the value entry when line contains comment and filters out space characters.
                else:
                    Parameter = GCFLine[idxA:]                                 # Extracts the value entry oterwise and filters out space characters.
                if Parameter and Label == Key.lower():                         # Checks for empty key value.
                    if Value[1] == 0:                                          # Checks for data type (0:boolean, 1: integer, 2: float, 3: string, 4: list).
                        try: 
                            if Label in LocalParameters:
                                Value[0].extend([bool(strtobool(Parameter))])  # Adds another user defined integer key value in case of local parameter.
                            else:
                                Value[0] = bool(strtobool(Parameter))          # Overwrites the default key value with user defined boolean key value.
                        except:
                            Log.append("Error: Invalid value for configuratio"\
                                       "n key '" + Key + "'. The value should"\
                                       " be either 'True', 'Yes, 'On', '1', '"\
                                       "False', 'No', 'Off', or '0'")          # Raises error if the key value is of invalid format.
                            writeToLogFile.write(Log,GCF['Name'][0])           # The saved file can be located in the working directory.
                            raise Exception("Fatal error occured, see " +     \
                                            GCF['Name'][0] + ".log file for d"\
                                            "etails")
                    elif Value[1] == 1:                                        # Checks for data type (0:boolean, 1: integer, 2: float, 3: string, 4: list).
                        try: 
                            if Label in LocalParameters:
                                if Label not in SignedParameters:
                                    Value[0].extend([abs(int(Parameter))])     # Adds another user defined integer key value in case of local parameter.
                                else:    
                                    Value[0].extend([int(Parameter)])          # Adds another user defined integer key value in case of local parameter.
                            else:
                                Value[0] = abs(int(Parameter))                 # Overwrites the default key value with user defined integer key value.
                        except:
                            Log.append("Error: Invalid value for configuratio"\
                                       "n key '" + Key + "'. The value should"\
                                       " be a single integer number")          # Raises error if the key value is of invalid format.
                            writeToLogFile.write(Log,GCF['Name'][0])           # The saved file can be located in the working directory.
                            raise Exception("Fatal error occured, see " +     \
                                            GCF['Name'][0] + ".log file for d"\
                                            "etails")
                    elif Value[1] == 2:                                        # Checks for data type (0:boolean, 1: integer, 2: float, 3: string, 4: list).
                        try: 
                            if Label in LocalParameters:
                                if Label not in SignedParameters:
                                    Value[0].extend([abs(float(Parameter))])   # Adds another user defined integer key value in case of local parameter.
                                else:
                                    Value[0].extend([float(Parameter)])        # Adds another user defined integer key value in case of local parameter.
                            else:
                                Value[0] = abs(float(Parameter))               # Overwrites the default key value with user defined floating point key value.
                        except:
                            Log.append("Error: Invalid value for configuratio"\
                                       "n key '" + Key + "'. The value should"\
                                       " be a single floating point number")   # Raises error if the key value is of invalid format.
                            writeToLogFile.write(Log,GCF['Name'][0])           # The saved file can be located in the working directory.
                            raise Exception("Fatal error occured, see " +     \
                                            GCF['Name'][0] + ".log file for d"\
                                            "etails")
                    elif Value[1] == 3:                                        # Checks for data type (0:boolean, 1: integer, 2: float, 3: string, 4: list).
                        if Label in LocalParameters:
                            Value[0].extend([Parameter])                       # Adds another user defined integer key value in case of local parameter.
                        else:
                            Value[0] = Parameter                               # Overwrites the default key value with user defined string key value.
                    else:
                        separator = [',',';',';,',',;']
                        r = re.compile(r'\b(' + '|'.join(separator) + r')\b')
                        idxCD = [j.start() for j in r.finditer(Parameter)]     # Looks for the indicies of separators.
                        if idxCD == []:
                            Value[0].append([Parameter])
                        else:
                            nP = len(idxCD) + 1
                            idxCD.insert(0,-1); idxCD.append(len(Parameter))
                            Parameter = [Parameter[idxCD[j] + 1:idxCD[j + 1]] \
                                         for j in range(nP)]
                            Value[0].append(Parameter)
                    GCF[Key] = Value                                           # Writes the key value to the dictionary of GMT configuration options.
                elif (Label in LocalParametersSV) and (Label == Key.lower()):
                    Value[0].append([])
                    GCF[Key] = Value                                           # Writes the key value to the dictionary of GMT configuration options.
        i += 1
    
    # Validation of input key values for configuration options:
    if GCF['WorkingDirectoryPath'][0] is None:
        Log.append("Error: The GMT Configuration File " + GCFName + " is miss"\
                   "ing an absolute path to the working directory")            # Raises error if the GMT Configuration File lacks an absolute path to the working directory.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if GCF['Name'][0] is None:
        Log.append("Error: The GMT Configuration File " + GCFName + " is miss"\
                   "ing a name of the model")                                  # Raises error if the GMT Configuration File lacks a name of the model.
        writeToLogFile.write(Log,"Log")                                        # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if GCF['MeshDim'][0] not in [0,1,2,3]:
        Log.append("Error: Value of the key 'MeshDim' should be either "\
                   "'0', '1', '2', or '3'")                                    # Raises error if the parameter 'MeshDim' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if GCF['MeshAlgorithm'][0] not in [1,2,3,5,6,7,8,9]:
        Log.append("Error: Value of the key 'MeshAlgorithm' should be either "\
                   "'1', '2', '3', '5', '6', '7', '8', or '9'")                # Raises error if the parameter 'MeshAlgorithm' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if GCF['MeshAlgorithm3D'][0] not in [1,3,4,7,9,10]:
        Log.append("Error: Value of the key 'MeshAlgorithm3D' should be eithe"\
                   "r '1', '3', '4', '7', '9', or '10'")                       # Raises error if the parameter 'MeshAlgorithm3D' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if GCF['HighOrderOptimize'][0] not in [0,1,2,3,4]:
        Log.append("Error: Value of the key 'HighOrderOptimize' should be eit"\
                   "her '0', '1', '2', '3', or '4'")                           # Raises error if the parameter 'HighOrderOptimize' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    ni = sum([1 for i in GCF['InflationLayersSurfaces'][0] if len(i) > 0])
    if len(GCF['InflationLayers'][0]) != ni:
        Log.append("Error: The number of local parameters 'InflationLayers' i"\
                   "s not equal to the number of local parameters 'InflationL"\
                   "ayersSurfaces'")                                           # Raises error if the number of parameters 'InflationLayers' is not equal to the number of parameters 'InflationLayersSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if len(GCF['InflationLayersMethod'][0]) != ni:
        Log.append("Error: The number of local parameters 'InflationLayersMet"\
                   "hod' is not equal to the number of local parameters 'Infl"\
                   "ationLayersSurfaces'")                                     # Raises error if the number of parameters 'InflationLayersMethod' is not equal to the number of parameters 'InflationLayersSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if len(GCF['InflationLayersThickness'][0]) != ni:
        Log.append("Error: The number of local parameters 'InflationLayersThi"\
                   "ckness' is not equal to the number of local parameters 'I"\
                   "nflationLayersSurfaces'")                                  # Raises error if the number of parameters 'InflationLayersThickness' is not equal to the number of parameters 'InflationLayersSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if len(GCF['InflationLayersGrowthRate'][0]) != ni:
        Log.append("Error: The number of local parameters 'InflationLayersGro"\
                   "wthRate' is not equal to the number of local parameters '"\
                   "InflationLayersSurfaces'")                                 # Raises error if the number of parameters 'InflationLayersGrowthRate' is not equal to the number of parameters 'InflationLayersSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    nms = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshSurfaces'][0]])
    nmV = np.asarray([1 if len(i) > 0 else 0 for i in                         \
                      GCF['LocalMeshVolumes'][0]])
    if (len(nms) > 0) and (len(nmV) > 0):
        nm = sum(nms | nmV)
    else:
        nm = 0
    if len(GCF['LocalMeshSize'][0]) != nm:
        Log.append("Error: The number of local parameters 'LocalMeshSize' is "\
                   "not equal to the number of local parameters 'LocalMeshSur"\
                   "faces'")                                                   # Raises error if the number of parameters 'LocalMeshSize' is not equal to the number of parameters 'LocalMeshSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    if len(GCF['LocalMeshGrowthRate'][0]) != nm:
        Log.append("Error: The number of local parameters 'LocalMeshGrowthRat"\
                   "e' is not equal to the number of local parameters 'LocalM"\
                   "eshSurfaces'")                                             # Raises error if the number of parameters 'LocalMeshGrowthRate' is not equal to the number of parameters 'LocalMeshSurfaces'.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    
    # Additional treatment of local parameters values:
    if ni > 0:
        S = GCF['InflationLayersSurfaces'][0]
        n = len(S)
        S[0] = list(set(S[0]))
        # S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
        #       for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['InflationLayersSurfaces'][0] = S
    if nm > 0:
        S = GCF['LocalMeshSurfaces'][0]
        n = len(S)
        S[0] = list(set(S[0]))
        # S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
        #       for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['LocalMeshSurfaces'][0] = S
        S = GCF['LocalMeshVolumes'][0]
        n = len(S)
        S[0] = list(set(S[0]))
        # S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
        #       for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['LocalMeshVolumes'][0] = S
    
    # Launch of the respective meshing module:
    TModel = perf_counter() - TModel; TPreP = perf_counter() - TPreP
    if GCF['InputFormat'][0] == 1:
        GMTMesh.IGS(GCF,Log,TModel,TPreP)
    elif GCF['InputFormat'][0] == 2:
        GMTMesh.STL(GCF,Log,TModel,TPreP)
    else:
        Log.append("Error: Value of the key 'Format' should be either '1' or "\
                   "'2'")                                                      # Raises error if the parameter 'Format' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details")
    
    #Finalization:
    return
