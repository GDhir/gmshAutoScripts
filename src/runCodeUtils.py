import os
import subprocess
import re
import meshFileUtils
import fileNameUtils
import regexUtils
from os import listdir
from os.path import isfile, join

def buildMesh( gmshfilecmd, gmshfileargs ):

    gmshbuildFolder = "/home/gaurav/gmshAutoScripts/build/"
    compilecmd = ["make", "-j", "4", gmshfilecmd]
    subprocess.run( compilecmd, cwd = gmshbuildFolder )

    gmshfileargs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshlstcmd = [ gmshbuildFolder + gmshfilecmd, gmshfileargs ]
    subprocess.run( meshlstcmd, cwd = gmshbuildFolder )

def setFinchTriangleQuadrature( quadratureValue ):

    statement = "xyw = triangle_quadrature_nodes_weights(" +  str( quadratureValue ) + ");\n"
    filename = "/home/gaurav/.julia/packages/Finch/ECEMc/src/triangle_nodes.jl"
    pattern = "triangle_quadrature_nodes_weights"

    with open( filename ) as fileHandle:

        allLines = fileHandle.readlines()

        for idx, lineVal in enumerate( allLines ):

            if re.search( pattern, lineVal ):

                allLines[idx] = statement
                break

    with open( filename, "w" ) as fileHandle:

        fileHandle.writelines( allLines )

def setFinchQuadQuadrature( quadratureValue ):

    statement = "quadratureOrder = " + str( quadratureValue ) + ";\n"
    filename = "/home/gaurav/.julia/packages/Finch/ECEMc/src/refel.jl"
    pattern = "quadratureOrder"

    with open( filename ) as fileHandle:

        allLines = fileHandle.readlines()

        for idx, lineVal in enumerate( allLines ):

            if re.search( pattern, lineVal ):

                allLines[idx] = statement
                break

    with open( filename, "w" ) as fileHandle:

        fileHandle.writelines( allLines )

def setFinchGeneratedSolveFunction( optionsParam, filename, pattern, replacement ):

    with open( filename ) as fileHandle:

        allLines = fileHandle.readlines()

        for idx, lineVal in enumerate( allLines ):

            if re.search( pattern, lineVal ):

                allLines[idx] = replacement
                break

    with open( filename, "w" ) as fileHandle:

        fileHandle.writelines( allLines )


def runFinchSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams, comparisonParam, meshPath, gaussPointsViz = False ):
    
    levelsArr = meshFileUtils.getAllLevels( meshArr )

    for paramValue in allParams[ comparisonParam ]:
        optionsParam[ comparisonParam ] = paramValue
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )

        for levelsVal in levelsArr:

            optionsParam[ "level" ] = str(levelsVal)
            
            regexCriteriaVals = [ str( levelsVal ) ]

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            print(meshFileName)

            if re.search( "triangle", optionsParam["meshRegexVal"] ):
                setFinchTriangleQuadrature( optionsParam[ "quadratureOrder" ] )
            elif re.search( "regular", optionsParam[ "meshRegexVal" ] ):
                setFinchQuadQuadrature( optionsParam[ "quadratureOrder" ] )
            elif re.search( "mesh", optionsParam[ "meshRegexVal" ] ):
                setFinchTriangleQuadrature( optionsParam[ "quadratureOrder" ] )

                filename = "/home/gaurav/Finch/src/examples/mixed2dcodein.jl"
                pattern = "function genfunction_1"
                returnStr = "return (" + optionsParam["coeff_F"] + "*pi*pi*sin(" + optionsParam[ "sin(kpix)" ] + \
                    "*pi*x)*sin(" + optionsParam[ "sin(kpix)" ] + "*pi*y));"
                replacement = "function genfunction_1(x::Union{Float64,Float64},y::Union{Float64,Float64},z::Union{Float64,Float64},t::Union{Float64,Float64},node_index::Int,face_index::Int,indices::Vector{Int}) " + returnStr + " end\n"

                setFinchGeneratedSolveFunction( optionsParam, filename, pattern, replacement )

            if gaussPointsViz:
                filename = "/home/gaurav/Finch/src/examples/mixed2dcodein.jl"
                pattern = "pythonVarName ="
                pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
                replacement = "pythonVarName = \"" + pythonVarName + "\"\n"

                setFinchGeneratedSolveFunction( optionsParam, filename, pattern, replacement )

            runFinchSimWithOptions( optionsParam, meshFileName )

def runFinchSimWithOptions( optionsParam, meshval ):

    exefilename = "/home/gaurav/Finch/src/examples/example-mixed-element-2d.jl"
    pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
    allFinchOptions = [ optionsParam["sin(kpix)"], optionsParam["coeff_F"], meshval, pythonVarName ]

    runJulia( exefilename, allFinchOptions )
    # runJulia( exefilename )

def runJulia( exefilename, options = [] ):

    # juliapath = "/home/gaurav/Downloads/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    juliapath = "/home/gaurav/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    finchPath = "/home/gaurav/Finch/src/examples/"
    julialstcmd = [juliapath, exefilename]

    for option in options:
        julialstcmd.append( option )

    subprocess.run( julialstcmd, cwd = finchPath, env = os.environ.copy() )

def removeFiles( dirval ):

    meshvals = [join(dirval, f) for f in listdir(dirval) if isfile(join(dirval, f))]

    for meshfilename in meshvals:
        cmdlist = [ "rm", meshfilename ]
        subprocess.run( cmdlist )

def buildAllMeshes( gmshFileCmdNames, gmshfileArgs ):

    for gmshFileCmd in gmshFileCmdNames:
        buildMesh( gmshFileCmd, gmshfileArgs )

def runFinchSim():

    exefilename = "example-mixed-element-2d.jl"
    runJulia( exefilename )

def runDealiiSim():

    exefilename = "step-5.debug"
    dealiiPath = "/home/gaurav/dealii-9.5.1/build/bin/"

    subprocess.run( [dealiiPath + exefilename], cwd = dealiiPath )

def setDealiiDefinition( srcFileName, patternStr, replacementStr ):

    with open( srcFileName ) as fileHandle:

        allLines = fileHandle.readlines()

        for idx, lineVal in enumerate( allLines ):

            if re.search( patternStr, lineVal ):

                allLines[idx] = replacementStr
                break

    with open( srcFileName, "w" ) as fileHandle:

        fileHandle.writelines( allLines )

def setDealiiOptions( srcFileName, optionsParam ):

    kval = optionsParam[ "sin(kpix)" ]
    patternStr = "#define k"
    replacementStr = patternStr + " " + kval + "\n"
    setDealiiDefinition( srcFileName, patternStr, replacementStr )

    cval = optionsParam[ "coeff_F" ]
    patternStr = "#define c"
    replacementStr = patternStr + " " + cval + "\n"
    setDealiiDefinition( srcFileName, patternStr, replacementStr )

def buildDealii( makeArg ):

    makeCmd = [ "make", makeArg ]
    dealiiBuildPath = "/home/gaurav/dealii-9.5.1/build/examples/"

    subprocess.run( makeCmd, cwd = dealiiBuildPath )

def runDealiiSimWithOptions( optionsParam, meshval, srcFileName, exefilename = "step-5.debug",
                            makeArg = "example_step_5_debug"):

    dealiiBinPath = "/home/gaurav/dealii-9.5.1/build/bin/"
    # exefilename = "step-5.debug"

    buildDealii( makeArg )

    pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
    subprocess.run( [ dealiiBinPath + exefilename, meshval, pythonVarName,
                      optionsParam["quadratureOrder"] ], cwd = dealiiBinPath )

def runDealiiSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams, comparisonParam, 
                                        meshPath, srcFileName, exefilename = "step-5.debug",
                                        makeArg = "example_step_5_debug" ):
    
    levelsArr = meshFileUtils.getAllLevels( meshArr )

    for paramValue in allParams[ comparisonParam ]:
        optionsParam[ comparisonParam ] = paramValue
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )

        for levelsVal in levelsArr:

            optionsParam[ "level" ] = str(levelsVal)
            regexCriteriaVals = [ str(levelsVal) ]
            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            setDealiiOptions( srcFileName, optionsParam )
            runDealiiSimWithOptions( optionsParam, meshFileName, srcFileName, exefilename, makeArg )

def runSim( simPlotFolderName, gmshFileCmdNames, regexVals,
            gmshFolderPath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/",
            removeFileOption = True, 
            buildMeshes = True ):

    if removeFileOption:
        removeFiles( gmshFolderPath )

    if buildMeshes:
        buildAllMeshes( gmshFileCmdNames, gmshFolderPath )

    # runFinchSim() 
    runDealiiSim()   

    # showFinchPlot( simPlotFolderName, regexVals )
    # showParaviewPlot( simPlotFolderName, regexVals )
    # removeFiles( gmshfileargs )
