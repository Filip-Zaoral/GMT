
def Import():
    
    # Intitialization:
    from glob import glob
    from distutils import util
    from Lib import GmshMesh, writeToLogFile
    from time import perf_counter
    
    # Declaration of default Gmsh configuration options:
    Log = []
    GCF = {
           'WorkingDirectoryPath' : (None,3),                                  # Absolute path to the working directory.                                                               [string]
           'Name' : (None,3),                                                  # Model Name (IGS: <Model Name>.<Volume Name>.igs, STL: <Model Name>.<Volume Name>.<Surface Name>.stl). [string]
           'Format' : (1,1),                                                   # Format of the model files (1: IGES, 2: STL).                                                          [integer]
           'STLRemesh' : (False,0),                                            # Parametrize and remesh geometry of the model when using STL file format (Format = 2)?                 [boolean]
           'STLMaxNumCores' : (1,1),                                           # Sets the maximum number of CPU cores to use for meshing of individual shells during remeshing
                                                                               # (STLRemesh = True).                                                                                   [integer]
           'STLFacetAngle' : (40.,2),                                          # Angle between two facets above which an edge is considered as sharp.                                  [float]
           'STLCurveAngle' : (40.,2),                                          # Angle between two curve segments above which an edge is considered as sharp.                          [float]
           'GeometryTolerance' : (None,2),                                     # Geometrical tolerance.                                                                                [float]
           'MeshAlgorithm' : (2,1),                                            # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
                                                                               # 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms).           [integer]
           'MeshAlgorithm3D' : (1,1),                                          # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT).      [integer]
           'MeshSizeFromPoints' : (False,0),                                   # Compute mesh element sizes from values given at geometry points?                                      [boolean]
           'MeshSizeFromCurvature' : (False,0),                                # Automatically compute mesh element sizes from curvature?                                              [boolean]
           'MeshSizeExtendFromBoundary' : (True,0),                            # Extend computation of mesh element sizes from the boundaries into the interior?                       [boolean]
           'MeshSizeMax' : (None,2),                                           # Sets the maximum global element size.                                                                 [float]
           'MeshSizeMin' : (None,2),                                           # Sets the minimum global element size.                                                                 [float]
           'MinimumCirclePoints' : (8,1),                                      # Sets the minimum number of nodes used to mesh circles and ellipses.                                   [integer]
           'MinimumCurvePoints' : (3,1),                                       # Sets the minimum number of points used to mesh curves other than lines, circles and ellipses.         [integer]
           'MinimumElementsPerTwoPi' : (6,1),                                  # Sets the minimum number of elements per 2*Pi radians when the mesh size is adapted to the curvature.  [integer]
           'ElementOrder' : (1,1),                                             # Sets the element order.                                                                               [integer]
           'SecondOrderIncomplete' : (True,0),                                 # Create incomplete second order elements (8-node quads, 20-node hexas, etc.)?                          [boolean]
           'SecondOrderLinear' : (False,0),                                    # Should second order nodes, and nodes generated through subdivision
                                    										   # algorithms, simply be created by linear interpolation?                                                [boolean]
           'Optimize' : (True,0),                                              # Optimize the mesh to improve the quality of tetrahedral elements?                                     [boolean]
           'OptimizeNetgen' : (False,0),                                       # Optimize the mesh using Netgen library to improve the quality of tetrahedral elements?                [boolean]
           'HighOrderOptimize' : (1,1),                                        # High-order mesh optimization algorithm (0: none, 1: optimization,
                                                                               # 2: elastic+optimization, 3: elastic, 4: fast curving).                                                [integer]
           'MaxNumThreads1D' : (1,1),                                          # Sets the maximum number of CPU threads to use for meshing of edges.                                   [integer]
           'MaxNumThreads2D' : (1,1),                                          # Sets the maximum number of CPU threads to use for meshing of surfaces.                                [integer]
           'MaxNumThreads3D' : (1,1),                                          # Sets the maximum number of CPU threads to use for meshing of volumes.                                 [integer]
           'ScalingFactor' : (1.,2),                                           # Global scaling factor applied to the saved mesh.                                                      [float]
           'Binary' : (False,0),                                               # Write mesh files in binary format (if possible)?                                                      [boolean]
           'LaunchGmsh' : (False,0)                                            # Launch Gmsh after finishing the task?                                                                 [boolean]
           }
    
    # Import of configuration options from the Gmsh Configuration File:
    TModel = perf_counter(); TPreP = perf_counter()
    GCFList = glob('*.gcf')                                                    # Obtains a list of names of all Gmsh Configuration Files in the installation directory.
    if len(GCFList) == 0:                                                      # Checks if there is any Gmsh Configuration File in the installation directory at all.
        Log.append("Error: No Gmsh Configuration File (.gcf) found in the roo"\
                   "t folder of the Gmsh Meshing Tool library.")               # Raise error if there is no Gmsh Configuration File in the installation directory.
        writeToLogFile.write(Log,"Log")                                        # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    elif len(GCFList) > 1:                                                     # Checks if there are more than one Gmsh Configuration File in the installation directory.
        Log.append("Warning: Multiple following Gmsh Configuration Files foun"\
                   "d in the library directory:")
        Log += ["Warning:  - " + i for i in GCFList]
        Log.append(str(GCFList[0]) + " was selected to continue.")
        GCFName = GCFList[0]                                                   # Obtain name of the chosen Gmsh Configuration File.
    else:
        GCFName = GCFList[0]                                                   # Obtain name of Gmsh Configuration File otherwise.
    GCFRaw = open(GCFName,'r')                                                 # Opens the Gmsh Configuration File for reading.
    GCFLines = GCFRaw.readlines()                                              # Reads the Gmsh Configuration File.
    GCFRaw.close()                                                             # Closes the Gmsh Configuration File.
    nGCFLines = len(GCFLines)
    for Key, Value in GCF.items():                                             # Loops over key - value pairs in the dictionary of default Gmsh configuration options.
        for i in range(nGCFLines):                                             # Loops over key - value pairs in the user specified Gmsh configuration file.
            idxA = GCFLines[i].find(Key)                                       # Looks for keys in the Gmsh configuration file.
            if idxA != -1:                                                     # Checks for lines with keys.
                break
        if idxA != -1:
            idxA = GCFLines[i].find('=') + 1                                   # Looks for the begining index of the value entry in the line.
            idxB = GCFLines[i].find('#')                                       # Looks for the endg index of the value entry in the line.
            if idxB != -1:
                Parameter = GCFLines[i][idxA:idxB].replace(" ","")             # Extracts the value entry when line contains comment and filters out space characters.
            else:
                Parameter = GCFLines[i][idxA:].replace(" ","")                 # Extracts the value entry oterwise and filters out space characters.
            if Parameter:                                                      # Checks for empty key value.
                if Value[1] == 0:                                              # Checks for data type (0:boolean, 1: integer, 2: float, 3: string).
                    try: 
                        Value = list(Value)
                        Value[0] = bool(util.strtobool(Parameter))             # Overwrites the default key value with user defined boolean key value.
                        Value = tuple(Value)
                    except:
                        Log.append("Error: Invalid value for configuration ke"\
                                   "y '" + Key + "'. The value should be eith"\
                                   "er 'True', 'Yes, 'On', '1', 'False', 'No'"\
                                   ", 'Off', or '0'.")                         # Raises error if the key value is of invalid format.
                        writeToLogFile.write(Log,GCF['Name'][0])               # The saved file can be located in the working directory.
                        raise Exception("Fatal error occured, see " + GCF['Na'\
                                        'me'][0] + ".log file for details.")
                elif Value[1] == 1:                                            # Checks for data type (0:boolean, 1: integer, 2: float, 3: string).
                    try: 
                        Value = list(Value)
                        Value[0] = abs(int(Parameter))                         # Overwrites the default key value with user defined integer key value.
                        Value = tuple(Value)
                    except:
                        Log.append("Error: Invalid value for configuration ke"\
                                   "y '" + Key + "'. The value should be sing"\
                                   "le integer number.")                       # Raises error if the key value is of invalid format.
                        writeToLogFile.write(Log,GCF['Name'][0])               # The saved file can be located in the working directory.
                        raise Exception("Fatal error occured, see " + GCF['Na'\
                                        'me'][0] + ".log file for details.")
                elif Value[1] == 2:                                            # Checks for data type (0:boolean, 1: integer, 2: float, 3: string).
                    try: 
                        Value = list(Value)
                        Value[0] = abs(float(Parameter))                       # Overwrites the default key value with user defined floating point key value.
                        Value = tuple(Value)
                    except:
                        Log.append("Error: Invalid value for configuration ke"\
                                   "y '" + Key + "'. The value should be sing"\
                                   "le floating point number.")                # Raises error if the key value is of invalid format.
                        writeToLogFile.write(Log,GCF['Name'][0])               # The saved file can be located in the working directory.
                        raise Exception("Fatal error occured, see " + GCF['Na'\
                                        'me'][0] + ".log file for details.")
                else:                                                          # Checks for data type (0:boolean, 1: integer, 2: float, 3: string).
                    Value = list(Value)
                    Value[0] = Parameter                                       # Overwrites the default key value with user defined string key value.
                    Value = tuple(Value)
                GCF[Key] = Value                                               # Writes the key value to the dictionary of Gmsh configuration options.
        i += 1
    
    # Additional treatment of local parameters values:
    if GCF['InflationLayersSurfaces'][0]:
        S = GCF['InflationLayersSurfaces'][0]
        n = len(S)
        S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
              for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['InflationLayersSurfaces'][0] = S
    if GCF['LocalMeshSurfaces'][0]:
        S = GCF['LocalMeshSurfaces'][0]
        n = len(S)
        S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
              for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['LocalMeshSurfaces'][0] = S
    if GCF['LocalMeshVolumes'][0]:
        S = GCF['LocalMeshVolumes'][0]
        n = len(S)
        S = [[j for j in S[i] if j not in [jj for ii in S[:i] for jj in ii]]  \
              for i in range(n)]                                               # Removes duplicate surface labels across different lists.
        GCF['LocalMeshVolumes'][0] = S
    
    # Validation of input key values for configuration options:
    if GCF['WorkingDirectoryPath'][0] is None:
        Log.append("Error: The Gmsh Configuration File " + GCFName + " is mis"\
                   "sing an absolute path to the working directory.")          # Raises error if the Gmsh Configuration File lacks an absolute path to the working directory.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    if GCF['Name'][0] is None:
        Log.append("Error: The Gmsh Configuration File " + GCFName + " is mis"\
                   "sing a name of the model.")                                # Raises error if the Gmsh Configuration File lacks a name of the model.
        writeToLogFile.write(Log,"Log")                                        # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    if GCF['STLMaxNumCores'][0] == 0:
        Log.append("Error: Value of the key 'STLMaxNumCores' must be larger t"\
                   "han '0'.")                                                 # Raises error if the parameter 'STLMaxNumCores' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    if GCF['MeshAlgorithm'][0] not in [1,2,3,5,6,7,8,9]:
        Log.append("Error: Value of the key 'MeshAlgorithm' should be either "\
                   "'1', '2', '3', '5', '6', '7', '8', or '9'.")               # Raises error if the parameter ''MeshAlgorithm' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    if GCF['MeshAlgorithm3D'][0] not in [1,3,4,7,9,10]:
        Log.append("Error: Value of the key 'MeshAlgorithm3D' should be eithe"\
                   "r '1', '3', '4', '7', '9', or '10'.")                      # Raises error if the parameter 'MeshAlgorithm3D' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    if GCF['HighOrderOptimize'][0] not in [0,1,2,3,4]:
        Log.append("Error: Value of the key 'HighOrderOptimize' should be eit"\
                   "her '0', '1', '2', '3', or '4'.")                          # Raises error if the parameter 'HighOrderOptimize' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    
    # Launch of the respective meshing module:
    TModel = perf_counter() - TModel; TPreP = perf_counter() - TPreP
    if GCF['Format'][0] == 1:
        GmshMesh.IGS(GCF,Log,TModel,TPreP)
    elif GCF['Format'][0] == 2:
        if GCF['STLRemesh'][0] == True:
            # Log,TModel,TPreP,TMesh = GmshMesh.RemeshSTL(GCF,Log,TModel,TPreP)
            Log,TModel,TPreP,TMesh = GmshMesh.RemeshSTLSeriall(GCF,Log,TModel,TPreP)
        else:
            TMesh = 0
        GmshMesh.STL(GCF,Log,TModel,TPreP,TMesh)
    else:
        Log.append("Error: Value of the key 'Format' should be either '1' or "\
                   "'2'.")                                                     # Raises error if the parameter 'Format' contains invalid value.
        writeToLogFile.write(Log,GCF['Name'][0])                               # The saved file can be located in the working directory.
        raise Exception("Fatal error occured, see " + GCF['Name'][0] + ".log "\
                        "file for details.")
    
    #Finalization:
    return
