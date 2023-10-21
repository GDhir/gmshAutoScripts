import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import re
from os import listdir
from os.path import isfile, join
from math import sin, cos, pi
import os
import subprocess
import meshio
import pdb
import time
import h5py
import folderUtils

def getSortedMeshVals( meshvals, regexVal, regexCriterias ):

    keyVals = []
    prevIndexMap = dict()

    for index, meshval in enumerate(meshvals):

        if re.search( regexVal, meshval  ):

            regexCriteriaVals = []

            for regexCriteria in regexCriterias:
                rval = re.search( regexCriteria, meshval )
                offset = rval.start() + len( regexCriteria )
                rval = re.search( "[0-9]+", meshval[offset:] )
                rval = rval.group(0)

                regexCriteriaVals.append( rval )

            keyVal = 1

            for idx, regexCriteria in enumerate( regexCriterias ):
                regexCriteriaVal = int( regexCriteriaVals[ idx ] )
                keyVal = keyVal * regexCriteriaVal

            keyVals.append( keyVal )
            prevIndexMap[ keyVal ] = index            

    keyVals = sorted( keyVals )
    sortedMeshVals = []

    for keyVal in keyVals:

        idx = prevIndexMap[keyVal]
        sortedMeshVals.append( meshvals[idx] )

    return sortedMeshVals

def getFileNameFromMeshName( meshval, folderName, varName, regexTypeVal, regexCriterias, extension ):

    regexCriteriaVals = []

    for regexCriteria in regexCriterias:
        rval = re.search( regexCriteria, meshval )

        if rval:
            offset = rval.start() + len( regexCriteria )
            rval = re.search( "[0-9]+", meshval[offset:] )
            rval = rval.group(0)

            regexCriteriaVals.append( rval )

    fileName = folderName + regexTypeVal + "_" + varName + "_"

    for idx, regexCriteria in enumerate( regexCriterias ):

        regexCriteriaVal = regexCriteriaVals[ idx ]

        fileName = fileName + regexCriteria + "=" + regexCriteriaVal

    fileName = fileName + extension

    return fileName

def getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName ):

    nLevels = len( levelsArr )
    minVals = np.ones( nLevels )*100000
    maxVals = np.ones( nLevels )*-1

    for optionIdx, option in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = option

        for idx, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )

            pythonVarName = getPythonVarName( optionsParam )
            textFile = getTextFileName( folderUtils.finchTextfoldername, pythonVarName, juliaVarName )

            errVals = getData( textFile )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )
    minMaxRangeVals = []

    for idx in range( nLevels ):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals

def getDealiiData( fileName, format = "hdf" ):

    if format == "hdf":
        f = h5py.File( fileName, 'r')
        nodes = np.array( f["nodes"] )
        solution = np.array( f["solution"] )
    elif format == "vtu":
        mesh = meshio.read( fileName )
        nodes = mesh.points
        # print(nodes)
        solution = mesh.point_data['solution']

    return (nodes, solution)

def getExactSol( xval, yval, negative = -1, pival = 2*pi ):

    return negative*sin(pival*xval)*sin(pival*yval)

def getDealiiError( nodes, solution, negative = 1, pival = 2*pi ):

    errVals = []
    exactvals = []

    for idx, node in enumerate( nodes ):

        xval = node[0]
        yval = node[1]
        # print(xval, yval)
        exactval = getExactSol( xval, yval, negative, pival )
        exactvals.append( exactval )

        if len( solution.shape ) == 2:
            solval = solution[ idx, 0 ]
        elif len( solution.shape ) == 1:
            solval = solution[ idx ]

        errorval = np.abs( exactval - solval )
        errVals.append( errorval )

    # print(errVals)

    exactvals = np.array(exactvals)

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], exactvals )

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], solution )
    # plt.colorbar()
    # plt.savefig( "solutionFile.png" )

    return np.array( errVals )

def getMeshFileName( optionsParam, criterias, criteriaVals, meshPath ):

    criteriaStr = getCriteriaValsString( criterias, criteriaVals )
    fileName = meshPath + optionsParam[ "meshRegexVal" ] + criteriaStr + ".msh"

    return fileName

def getDealiiMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, negative = 1, pival = 2*pi ):

    nLevels = len( levelsArr )
    minVals = np.ones( nLevels )*100000
    maxVals = np.ones( nLevels )*-1

    dataFileFormat = "vtu"

    for optionIdx, option in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = option

        for idx, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )

            pythonVarName = getPythonVarName( optionsParam )
            textFile = getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, "solutionvalues", dataFileFormat )

            nodes, solution = getDealiiData( textFile, dataFileFormat )
            errVals = getDealiiError( nodes, solution, negative, pival )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )

    minMaxRangeVals = []

    for idx in range( nLevels ):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals

def getFinchError( xvals, yvals, uvals, pival = 2*pi ):

    errVals = []
    exactvals = []

    for idx, xval in enumerate( xvals ):

        yval = yvals[idx]
        solval = uvals[idx]
        # print(xval, yval)
        exactval = getExactSol( xval, yval, 1, pival )
        exactvals.append( exactval )

        errorval = np.abs( exactval - solval )
        errVals.append( errorval )

    # print(errVals)

    exactvals = np.array(exactvals)

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], exactvals )

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], solution[:, 0] )

    return np.array( errVals )

def getData( filename ):

    uvals = []

    with open(filename) as fval:

        uvals = fval.readlines()

    # print(uvals)
    uvals = uvals[0].split(",")

    uvals[0] = uvals[0][1:]
    uvals[-1] = uvals[-1][:-1]

    for idx, val in enumerate(uvals):

        uvals[idx] = float( uvals[idx] )

    return uvals

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

def getRegexVal( fileName, allRegexes ):

    for regexVal in allRegexes:
        if re.search( regexVal, fileName ):

            return regexVal


def runFinchSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams ):
    
    for meshval in meshArr:

        regexVal = getRegexVal( meshval, allParams[ "meshRegexVal" ] )
        optionsParam[ "meshRegexVal" ] = regexVal

        levelVal = getLevelFromFileName( meshval )
        optionsParam[ "level" ] = levelVal

        print( regexVal, levelVal )
        runFinchSimWithOptions( optionsParam, meshval )

def runFinchSimWithOptions( optionsParam, meshval ):

    exefilename = "/home/gaurav/Finch/src/examples/example-mixed-element-2d.jl"
    setFinchTriangleQuadrature( optionsParam[ "quadratureOrder" ] )
    pythonVarName = getPythonVarName( optionsParam )
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

    setDealiiOptions( srcFileName, optionsParam )

    buildDealii( makeArg )

    pythonVarName = getPythonVarName( optionsParam )
    subprocess.run( [ dealiiBinPath + exefilename, meshval, pythonVarName,
                      optionsParam["quadratureOrder"] ], cwd = dealiiBinPath )
    
def runDealiiSimWithOptionsVariousMeshes( optionsParam, meshArr, srcFileName, exefilename = "step-5.debug",
                            makeArg = "example_step_5_debug" ):
    
    for meshval in meshArr:

        meshRegexVal = getRegexVal( meshval, allParams["meshRegexVal"] )
        levelVal = getLevelFromFileName( meshval )

        optionsParam[ "meshRegexVal" ] = meshRegexVal
        optionsParam[ "level" ] = levelVal
        runDealiiSimWithOptions( optionsParam, meshval, srcFileName, exefilename, makeArg )

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

def getRegexCriterias( regexVal ):

    # if regexVal == "hanging":
    #     regexCriterias = [ "Nx=", "Ny=" ]
    # elif regexVal == "regular" or re.search( "triangle", regexVal ):
    #     regexCriterias = [ "N=" ]
    # elif regexVal == "mesh":
    #     regexCriterias = [ "lvl" ]

    regexCriterias = ["lvl"]
    
    return regexCriterias

def getNumNodes( gmshFileName ):

    numNodes = 0

    with open( gmshFileName ) as filehandle:

        idx = 0

        for line in filehandle.readlines():

            if idx == 4:

                numNodes = int( line )

                break

            idx += 1

    return numNodes

def getDistance( vertex1, vertex2 ):

    distVal = 0

    for i, x in enumerate( vertex1 ):

        distVal += ( vertex1[ i ] - vertex2[ i ] )**2

    distVal = np.sqrt( distVal )

    return distVal

def getMaxH( gmshFileName, regexVal ):

    allNodes = getAllNodes( gmshFileName )
    elementNodeIndices = getElementNodeIndices( gmshFileName, regexVal )

    # print( len( allNodes ) )

    h = -1

    for nodeIndices in elementNodeIndices:
        # print( nodeIndices )
        for i, idx1 in enumerate( nodeIndices ):

            numIndices = len( nodeIndices )
            vertex1Coords = allNodes[ idx1 - 1 ]
            vertex2Coords = allNodes[ nodeIndices[ ( i + 1 )%numIndices ] - 1 ]

            distVal = getDistance( vertex1Coords, vertex2Coords )

            if distVal > h:
                h = distVal

    return h

def getAllNodes( gmshFileName ):

    allNodes = []

    with open( gmshFileName ) as filehandle:

        idxStart = 0
        idxEnd = 0

        idx = 0

        allLines = filehandle.readlines()

        for lineval in allLines:

            if lineval == "$Nodes\n":

                idxStart = idx + 2

            if lineval == "$EndNodes\n":

                idxEnd = idx
                break

            idx += 1
            
        nodeLines = allLines[ idxStart : idxEnd ]

        for nodeLine in nodeLines:

            nodeLine = [ float(node) for node in nodeLine.split() ]
            allNodes.append( nodeLine[ 1: ] )

    return allNodes

def getElementNodeIndices( gmshFileName, regexVal ):

    if regexVal == "regular":
        types = ["3"]

    if re.search( "triangle", regexVal ):
        types = ["2"]

    if re.search( "hanging", regexVal ):
        types = [ "2", "3" ]

    else:
        types = ["2", "3"]

    elementNodeIndices = []

    with open( gmshFileName ) as filehandle:

        idxStart = 0
        idxEnd = 0

        idx = 0

        allLines = filehandle.readlines()

        for lineval in allLines:

            if lineval == "$Elements\n":

                idxStart = idx + 2

            if lineval == "$EndElements\n":

                idxEnd = idx
                break

            idx += 1        

        elementLines = allLines[ idxStart : idxEnd ]

    for lineval in elementLines:

        lineval = lineval.split()

        if lineval[1] in types:

            ntags = lineval[2]
            indices = lineval[ 3 + int( ntags ): ]
            indices = [ int(index) for index in indices ]

            elementNodeIndices.append( indices )

    return elementNodeIndices

def get2DAreaFromNodes( coords ):

    nnodes = len( coords )

    area = 0

    for idx in range(nnodes):

        area = area + coords[ idx ][0] * coords[ ( idx + 1 ) % nnodes ][1] - coords[ (idx + 1)%nnodes ][ 0 ] * coords[ idx ][1]

    return abs( area/2 )


def getAverage2DArea( gmshFileName, regexVal ):

    allNodes = getAllNodes( gmshFileName )
    elementNodeIndices = getElementNodeIndices( gmshFileName, regexVal )

    averageArea = 0

    nElements = len( elementNodeIndices )

    for element in elementNodeIndices:

        coords = []

        for indexVal in element:
            coords.append( allNodes[ indexVal - 1 ][:2] )

        averageArea += get2DAreaFromNodes( coords )

    return averageArea/nElements


def makePlotAdjustmentsAndSave( axHandle, figHandle, plotFileName ):

    box = axHandle.get_position()
    axHandle.set_position([box.x0, box.y0 + box.height * 0.1,
                    box.width, box.height * 0.9])

    # Put a legend below current axis
    lgd = axHandle.legend(loc=9, bbox_to_anchor=(0.5, -0.05),
            fancybox=True, shadow=True, ncol=5)  
    
    text = axHandle.text(-0.2,1.05, "         ", transform=axHandle.transAxes)

    figHandle.savefig( plotFileName, bbox_extra_artists=( lgd, text ), bbox_inches = 'tight' )

def showFinchDerivPlot( simPlotFolderName, regexVals ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    # print(allSortedMeshVals)

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexVal = regexVals[idx]
        regexCriterias = getRegexCriterias( regexVal )

        for index, meshval in enumerate( sortedMeshVals ):

            numNodeVals.append( getNumNodes( meshpath + meshval ) )
            areaVals.append( getAverage2DArea( meshpath + meshval, regexVal ) )

            xvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "centroid_xvalues",
                                                     regexVal, regexCriterias, ".txt" )
            yvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "centroid_yvalues",
                                                    regexVal, regexCriterias, ".txt" )

            deriv_xvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "deriv_xvalues",
                                                          regexVal, regexCriterias, ".txt" )
            deriv_yvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "deriv_yvalues",
                                                          regexVal, regexCriterias, ".txt" )
            # errvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )

            xvals = getData( xvalsFileName )
            yvals = getData( yvalsFileName )
            deriv_xvals = getData( deriv_xvalsFileName )
            deriv_yvals = getData( deriv_yvalsFileName )

            plt.figure
            plt.tricontourf( xvals, yvals, deriv_xvals)
            plt.show()

def createMultipleFigureAxis( n ):

    axisArray = []
    figArray = []

    for idx in range(n):
        figArray.append( plt.figure() )
        axisArray.append( figArray[ idx ].add_subplot(1, 1, 1) )

    return [ figArray, axisArray ]

def makeMultiplePlotAdjustmentsAndSave( axisArray, figureArray, fileNames ):

    for idx, axVal in enumerate( axisArray ):
        makePlotAdjustmentsAndSave( axisArray[idx], figureArray[idx], fileNames[idx] )

def showFinchPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshvals, meshPath):
    
    folderUtils.checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    levelsArr = getAllLevels( meshvals )

    juliaVarName = "errorvalues"
    minMaxRangeVals = getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    allMaxErrorList = []
    allL2ErrorList = []

    figArrayNumNodes, axisArrayNumNodes = createMultipleFigureAxis( 2 )
    figArrayArea, axisArrayArea = createMultipleFigureAxis( 2 )
    figArrayDx, axisArrayDx = createMultipleFigureAxis( 2 )

    # print(allSortedMeshVals)

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexCriterias = getRegexCriterias( "lvl" )
        dxVals = []
        optionsParam[ comparisonParam ] = paramValue

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = getPythonVarName( optionsParam )
            regexCriteriaVals = [level]
            criteriaValsStr = getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, optionsParam[ "meshRegexVal" ] ) )
            areaVals.append( getAverage2DArea( meshFileName, optionsParam[ "meshRegexVal" ] ) )

            xvalsFileName = getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "xvalues" )
            yvalsFileName = getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "yvalues" )
            uvalsFileName = getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "uvalues" )
            
            xvals = getData( xvalsFileName )
            yvals = getData( yvalsFileName )
            uvals = getData( uvalsFileName )
            errvals = getFinchError( xvals, yvals, uvals, int( optionsParam[ "sin(kpix)" ] )*pi )

            # print(errvalsFileName)

            curmaxval = np.max(errvals)
            curminval = np.min(errvals)

            curMaxErrorList.append( curmaxval )
            curL2ErrorList.append( np.sqrt( np.sum( np.array(errvals)**2 )/numNodeVals[-1] ) ) 

            curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"

            if not os.path.exists( curPlotFolderName ):
                os.mkdir( curPlotFolderName )

            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
            # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
            # print( np.max(errvals) )

            minVal, maxVal = minMaxRangeVals[ index ]

            contourLevels = np.linspace( minVal, maxVal, 40 )
            plt.figure()
            plt.tricontourf( xvals, yvals, errvals, levels = contourLevels, vmin = minVal, vmax = maxVal, cmap = cm )
            plt.colorbar()
            plotVarName = "errorvaluesContour"
            plotfilename = getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            contourLevels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            plt.tricontourf( xvals, yvals, errvals, levels = contourLevels, colors = 'r')
            plt.colorbar()
            plotVarName = "large" + "errorvalues" + "Contour"
            plotfilename = getTextFileName( curPlotFolderName, pythonVarName, plotVarName, ".png" )
            plt.savefig( plotfilename )
            plt.close()
            
        allMaxErrorList.append( curMaxErrorList )
        allL2ErrorList.append( curL2ErrorList )

        plotMaxAndL2Error( axisArrayNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue )
        plotMaxAndL2Error( axisArrayArea, areaVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue ) 
        plotMaxAndL2Error( axisArrayDx, dxVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue )      

        h2Vals = [ dxVal**2 for dxVal in dxVals ]

        print( comparisonParam + "=" + paramValue )
        print( dxVals )
        print( curMaxErrorList )
        print( curL2ErrorList )

    xValsArray = [ numNodeVals, numNodeVals, areaVals, areaVals, dxVals, dxVals ]
    yValsArray = [ h2Vals ]*len(xValsArray)
    labelNameArray = [ "$h^2$" ]*len(xValsArray)

    axisArray = [ axisArrayNumNodes[0], axisArrayNumNodes[1],
                    axisArrayArea[0],
                    axisArrayArea[1],
                    axisArrayDx[0],
                    axisArrayDx[1] ]

    fileNameArray = [   simPlotFolderName + "maxErrorNumNodes" + ".png",
                        simPlotFolderName + "l2ErrorNumNodes" + ".png",
                        simPlotFolderName + "maxErrorArea" + ".png",
                        simPlotFolderName + "l2ErrorArea" + ".png",
                        simPlotFolderName + "maxErrorDx" + ".png",
                        simPlotFolderName + "l2ErrorDx" + ".png" ]

    figArray = [ figArrayNumNodes[0], figArrayNumNodes[1],
                figArrayArea[0],
                figArrayArea[1],
                figArrayDx[0],
                figArrayDx[1] ]


    plotMultipleAxis( axisArray, xValsArray, yValsArray, labelNameArray, "-x" )

    makeMultiplePlotAdjustmentsAndSave( axisArray, figArray, fileNameArray )

    # if len( regexVals ) > 1:
    #     errorDiffMax = [ allMaxErrorList[1][idx] - allMaxErrorList[0][idx] for idx in range( len(allMaxErrorList[0]) ) ]
    #     errorDiffL2 = [ allL2ErrorList[1][idx] - allL2ErrorList[0][idx] for idx in range( len(allL2ErrorList[0]) ) ]

    #     plt.figure()
    #     plt.plot( dxVals, errorDiffMax, "-o" )
    #     titleval = regexVals[1] + " - " + regexVals[0] + " Max Error Plot"
    #     plt.xlabel( "h" )
    #     plt.ylabel( "Error" )
    #     plt.title( titleval )

    #     plt.figure()
    #     titleval = regexVals[1] + " - " + regexVals[0] + " L2 Error Plot"
    #     plt.plot( dxVals, errorDiffL2, "-o" )
    #     plt.xlabel( "h" )
    #     plt.ylabel( "Error" )
    #     plt.title( titleval )

    plt.show()

    return

def plotMultipleAxis( axisArray, xvalsArray, yvalsArray, labelNameArray, markerVal ):

    for i, ax in enumerate( axisArray ):

        ax.loglog( xvalsArray[i], yvalsArray[i], markerVal, label = labelNameArray[i] )

def plotMaxAndL2Error( axisArray, xvals, maxError, l2Error, param, paramValue ):

    labelNameArray = [ param + "=" + paramValue + " $L^{\infty}$ Error_" + "h",
                       param + "=" + paramValue + " $L^{2}$ Error_" + "h" ]
    xvalsArray = [ xvals, xvals ]
    yvalsArray = [ maxError, l2Error ]

    plotMultipleAxis( axisArray, xvalsArray, yvalsArray, labelNameArray, "-o" )

def showDealiiPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshvals, meshPath, 
                   negative = 1, pival = 2*pi ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    levelsArr = getAllLevels( meshvals )

    # dealiiVarName = "errorvalues"
    # minMaxRangeVals = getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName )

    minMaxRangeVals = getDealiiMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, negative, pival )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    allMaxErrorList = []
    allL2ErrorList = []

    figArrayNumNodes, axisArrayNumNodes = createMultipleFigureAxis( 2 )
    figArrayArea, axisArrayArea = createMultipleFigureAxis( 2 )
    figArrayDx, axisArrayDx = createMultipleFigureAxis( 2 )

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexCriterias = getRegexCriterias( "lvl" )
        dxVals = []
        optionsParam[ comparisonParam ] = paramValue

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = getPythonVarName( optionsParam )
            regexCriteriaVals = [level]
            criteriaValsStr = getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, optionsParam[ "meshRegexVal" ] ) )
            areaVals.append( getAverage2DArea( meshFileName, optionsParam[ "meshRegexVal" ] ) )

            dealiiVarName = "solutionvalues"
            solvalsFileName = getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, dealiiVarName, "vtu" )
            
            (nodes, solution) = getDealiiData( solvalsFileName, "vtu" )
            xvals = nodes[:, 0]
            yvals = nodes[:, 1]
            errvals = getDealiiError( nodes, solution, negative, pival )
            # print(errvals)
            # print(solvalsFileName)

            curmaxval = np.max(errvals)
            curminval = np.min(errvals)

            curMaxErrorList.append( curmaxval )
            curL2ErrorList.append( np.sqrt( np.sum( np.array(errvals)**2 )/ numNodeVals[-1] ) ) 

            curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"

            folderUtils.checkAndCreateFolder( curPlotFolderName )

            minVal, maxVal = minMaxRangeVals[ index ]
            contourLevels = np.linspace( minVal, maxVal, 40 )
            plt.figure()
            plt.tricontourf( xvals, yvals, errvals, levels = contourLevels, vmin = minVal, vmax = maxVal, cmap = cm )
            plt.colorbar()
            plotVarName = "errorvaluesContour"
            plotfilename = getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            contourLevels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, errvals, levels = contourLevels, colors = 'r')
            plt.colorbar()
            plotVarName = "large" + "errorvalues" + "Contour"
            plotfilename = getTextFileName( curPlotFolderName, pythonVarName, plotVarName, ".png" )            
            plt.savefig( plotfilename )
            plt.close()
            
        allMaxErrorList.append( curMaxErrorList )
        allL2ErrorList.append( curL2ErrorList )

        plotMaxAndL2Error( axisArrayNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue )
        plotMaxAndL2Error( axisArrayArea, areaVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue ) 
        plotMaxAndL2Error( axisArrayDx, dxVals, curMaxErrorList, curL2ErrorList, comparisonParam, paramValue )      

        print( comparisonParam + "=" + paramValue )
        print( dxVals )
        print( curMaxErrorList )
        print( curL2ErrorList )

    h2Vals = [ dxVal**2 for dxVal in dxVals ]

    xValsArray = [ numNodeVals, numNodeVals, areaVals, areaVals, dxVals, dxVals ]
    yValsArray = [ h2Vals ]*len(xValsArray)
    labelNameArray = [ "$h^2$" ]*len(xValsArray)

    axisArray = [ axisArrayNumNodes[0], axisArrayNumNodes[1],
                    axisArrayArea[0],
                    axisArrayArea[1],
                    axisArrayDx[0],
                    axisArrayDx[1] ]

    fileNameArray = [   simPlotFolderName + "maxErrorNumNodes" + ".png",
                        simPlotFolderName + "l2ErrorNumNodes" + ".png",
                        simPlotFolderName + "maxErrorArea" + ".png",
                        simPlotFolderName + "l2ErrorArea" + ".png",
                        simPlotFolderName + "maxErrorDx" + ".png",
                        simPlotFolderName + "l2ErrorDx" + ".png" ]

    figArray = [ figArrayNumNodes[0], figArrayNumNodes[1],
                figArrayArea[0],
                figArrayArea[1],
                figArrayDx[0],
                figArrayDx[1] ]

    plotMultipleAxis( axisArray, xValsArray, yValsArray, labelNameArray, "-x" )

    makeMultiplePlotAdjustmentsAndSave( axisArray, figArray, fileNameArray )

    # if len( regexVals ) > 1:
    #     errorDiffMax = [ allMaxErrorList[1][idx] - allMaxErrorList[0][idx] for idx in range( len(allMaxErrorList[0]) ) ]
    #     errorDiffL2 = [ allL2ErrorList[1][idx] - allL2ErrorList[0][idx] for idx in range( len(allL2ErrorList[0]) ) ]

    #     plt.figure()
    #     plt.plot( dxvals, errorDiffMax, "-o" )
    #     titleval = regexVals[1] + " - " + regexVals[0] + " Max Error Plot"
    #     plt.xlabel( "h" )
    #     plt.ylabel( "Error" )
    #     plt.title( titleval )

    #     plt.figure()
    #     titleval = regexVals[1] + " - " + regexVals[0] + " L2 Error Plot"
    #     plt.plot( dxvals, errorDiffL2, "-o" )
    #     plt.xlabel( "h" )
    #     plt.ylabel( "Error" )
    #     plt.title( titleval )

    plt.show()

    return

def plotMaxAndL2ErrorComparison( axErrorHandle, xvals, maxError1,
                                 maxError2, l2Error1, l2Error2, h2vals, label1, label2 ):


    axErrorHandle[0].loglog( xvals, maxError1, "-o", label = label1 + "$L^{\infty}$ Error_" + "h" )
    axErrorHandle[0].loglog( xvals, maxError2, "-o", label = label2 + "$L^{\infty}$ Error_" + "h" )

    axErrorHandle[1].loglog( xvals, l2Error1, "-o", label = label1 + "$L^{2}$ Error_" + "h" )
    axErrorHandle[1].loglog( xvals, l2Error2, "-o", label = label2 + "$L^{2}$ Error_" + "h" )

    axErrorHandle[0].loglog( xvals, h2vals, "-x", label = "$h^2$" )
    axErrorHandle[1].loglog( xvals, h2vals, "-x", label = "$h^2$" )

def getMaxL2ErrorFigHandles():

    figMaxError = plt.figure()
    axMaxError = figMaxError.add_subplot(1, 1, 1)
    figL2Error = plt.figure()
    axL2Error = figL2Error.add_subplot(1, 1, 1)

    return [ figMaxError, figL2Error, axMaxError, axL2Error ]

def compareDealiiFinch( simPlotFolderName, regexVals, varName,
                       meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/" ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    # minMaxRangeVals = getDealiiMinMaxRange( allSortedMeshVals )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    dealiiAllMaxErrorList = []
    dealiiAllL2ErrorList = []

    finchAllMaxErrorList = []
    finchAllL2ErrorList = []

    # print(allSortedMeshVals)

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        dealiiCurMaxErrorList = []
        dealiiCurL2ErrorList = []

        regexVal = regexVals[idx]
        regexCriterias = getRegexCriterias( regexVal )

        finchCurMaxErrorList = []
        finchCurL2ErrorList = []
        numNodeVals = []
        areaVals = []
        dxVals = []

        figArrayNumNodes, axisArrayNumNodes = createMultipleFigureAxis( 2 )
        figArrayArea, axisArrayArea = createMultipleFigureAxis( 2 )
        figArrayDx, axisArrayDx = createMultipleFigureAxis( 2 )

        for index, meshval in enumerate( sortedMeshVals ):

            meshFileName = meshPath + meshval

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, regexVal ) )
            areaVals.append( getAverage2DArea( meshFileName, regexVal ) )

            dealiiSolvalsFileName = getFileNameFromMeshName( meshval, folderUtils.dealiiTextfoldername, varName,
                                                            regexVal, regexCriterias,  ".vtu" )
            
            (nodes, solution) = getDealiiData( dealiiSolvalsFileName, format = "vtu" )
            dealiiErrvals = getDealiiError( nodes, solution, negative=-1 )
            # print(solvalsFileName)

            dealiiCurmaxval = np.max(dealiiErrvals)

            dealiiCurMaxErrorList.append( dealiiCurmaxval )
            dealiiCurL2ErrorList.append( np.sqrt( np.sum( np.array( dealiiErrvals )**2 ) / numNodeVals[ -1 ] ) ) 

            # Finch data
            finchXvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "xvalues",
                                                         regexVal, regexCriterias, ".txt" )
            finchYvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "yvalues",
                                                         regexVal, regexCriterias, ".txt" )
            finchErrvalsFileName = getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "errorvalues",
                                                           regexVal, regexCriterias, ".txt" )

            finchXvals = getData( finchXvalsFileName )
            finchYvals = getData( finchYvalsFileName )
            finchErrvals = getData( finchErrvalsFileName )

            # print(errvalsFileName)

            finchCurmaxval = np.max( finchErrvals )
            finchCurMaxErrorList.append( finchCurmaxval )
            finchCurL2ErrorList.append( np.sqrt( np.sum( np.array( finchErrvals )**2 ) / numNodeVals[ -1 ] ) ) 

            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
            # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
            # print( np.max(errvals) )

        dealiiAllMaxErrorList.append( dealiiCurMaxErrorList )
        dealiiAllL2ErrorList.append( dealiiCurL2ErrorList )

        finchAllMaxErrorList.append( finchCurMaxErrorList )
        finchAllL2ErrorList.append( finchCurL2ErrorList )

        label1 = regexVals[idx] + " Dealii "
        label2 = regexVals[idx] + " Finch "

        h2Vals = [ dxVal**2 for dxVal in dxVals ]

        plotMaxAndL2ErrorComparison( axisArrayNumNodes, numNodeVals,
                                     dealiiCurMaxErrorList, finchCurMaxErrorList, dealiiCurL2ErrorList,
                                       finchCurL2ErrorList, h2Vals, label1, label2 )
        
        plotMaxAndL2ErrorComparison( axisArrayArea, areaVals, 
                                    dealiiCurMaxErrorList, finchCurMaxErrorList, dealiiCurL2ErrorList, 
                                    finchCurL2ErrorList, h2Vals, label1, label2 )
        
        plotMaxAndL2ErrorComparison( axisArrayDx, dxVals, 
                                    dealiiCurMaxErrorList, finchCurMaxErrorList, dealiiCurL2ErrorList, 
                                    finchCurL2ErrorList, h2Vals, label1, label2 )


        curPlotFolderName = simPlotFolderName + "Finch_Dealii_Comparison" + "/"
        folderUtils.checkAndCreateFolder( curPlotFolderName )

        axisArray = [ axisArrayNumNodes[0], axisArrayNumNodes[1],
                    axisArrayArea[0],
                    axisArrayArea[1],
                    axisArrayDx[0],
                    axisArrayDx[1] ]

        fileNameArray = [ curPlotFolderName + regexVal + "_maxErrorNumNodes" + ".png",
                        curPlotFolderName + regexVal + "_l2ErrorNumNodes" + ".png",
                        curPlotFolderName + regexVal + "_maxErrorArea" + ".png",
                        curPlotFolderName + regexVal + "_l2ErrorArea" + ".png",
                        curPlotFolderName + regexVal + "_maxErrorDx" + ".png",
                        curPlotFolderName + regexVal + "_l2ErrorDx" + ".png" ]

        figArray = [ figArrayNumNodes[0], figArrayNumNodes[1],
                    figArrayArea[0],
                    figArrayArea[1],
                    figArrayDx[0],
                    figArrayDx[1] ]

        makeMultiplePlotAdjustmentsAndSave( axisArray, figArray, fileNameArray )

    plt.show()

    return

def getPythonVarName( optionsParam ):

    pythonVarName = ""

    for param, paramValue in optionsParam.items():

        pythonVarName += param + "_" + paramValue

    return pythonVarName

def getTextFileName( folderName, pythonVarName, juliaVarName, extension = "txt" ):

    filename = folderName + pythonVarName + juliaVarName + "." + extension

    return filename

def getCriteriasFromFileName( meshval, regexCriterias ):

    regexCriteriaVals = ""

    for regexCriteria in regexCriterias:
        rval = re.search( regexCriteria, meshval )

        if rval:
            offset = rval.start() + len( regexCriteria )
            rval = re.search( "[0-9]+", meshval[offset:] )
            rval = rval.group(0)

            regexCriteriaVals += regexCriteria + "=" + rval

    return regexCriteriaVals[:-1]

def getCriteriaValsString( regexCriterias, criteriaVals ):

    regexCriteriaVals = ""

    for idx, regexCriteria in enumerate( regexCriterias ):
        regexCriteriaVals += regexCriteria + "=" + criteriaVals[idx]

    return regexCriteriaVals    

def getLevelFromFileName( meshval ):

    rval = re.search( "lvl", meshval )

    if rval:
        offset = rval.start() + len( "lvl" )
        rval = re.search( "[0-9]+", meshval[offset:] )
        rval = rval.group(0)

    return rval

def getAllLevels( meshFileNameArr ):

    levelsSet = set()

    for fileName in meshFileNameArr:

        level = getLevelFromFileName( fileName )
        levelsSet.add( level )

    return list( levelsSet )

def getMeshFilesFromFolder( meshpath ):
    
    meshvals = [ join( meshpath, f ) for f in listdir(meshpath) if isfile(join(meshpath, f))]
    return meshvals

if __name__ == "__main__":

    getMinMaxRangeFunc = dict()
    getMinMaxRangeFunc["Finch"] = getFinchMinMaxRange
    getMinMaxRangeFunc["Dealii"] = getDealiiMinMaxRange

    # gmshFileCmdNames = ["triangleMeshv1", "triangleMeshv2"]
    # regexVals = ["triangleMeshStruct", "triangleMeshUnstruct"]
    regexVals = ["triangleMeshUnstruct", "triangleMeshStruct", "regularMesh"]
    gmshFileCmdNames = ["triangleMeshv2", "triangleMeshv1", "regularMeshv3"]
    # gmshFileCmdNames = ["hangingMeshv8"]
    # regexVals = ["hanging"]

    allParams = dict()
    allParams["software"] = ["Finch", "Dealii"]
    allParams["meshRegexVal"] = regexVals
    allParams["gmshFileCmdName"] = gmshFileCmdNames
    allParams["quadratureOrder"] = [ "2", "3", "4" ]
    allParams["sin(kpix)"] = [ "1", "2", "4" ]
    allParams["coeff_F"] = [ "-2", "-8", "-32" ]

    optionsParam = dict()
    optionsParam["quadratureOrder"] = "2"
    optionsParam["sin(kpix)"] = "1"
    optionsParam["coeff_F"] = "-2"
    optionsParam["software"] = "Finch"
    optionsParam["meshRegexVal"] = regexVals[0]
    optionsParam["level"] = "0"

    # pythonVarName = getPythonVarName( optionsParam )
    # print(pythonVarName)
    # filename = getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, "solutionValues", "vtu" )
    # print(filename) 

    comparisonParam = "meshRegexVal"

    simPlotRootFolderName = folderUtils.gmshImageFolderName + "PlotNewSetup_2pi/"
    meshPlotRootFolderName = folderUtils.gmshImageFolderName + "MeshPlotsHangingLevel_QuadratureOrder=2_2pi/"

    # regexVals = [ "mesh" ]
    meshPath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    # showMeshes( folderUtils.meshPlotRootFolderName, regexVals )

    # runSim( simPlotRootFolderName, gmshFileCmdNames, regexVals )
    # createMeshVTU( meshPlotRootFolderName, regexVals )
    # compareDealiiFinch( simPlotRootFolderName, regexVals, "solutionvalues", meshPath )

    # simPlotFolderName = simPlotRootFolderName + "Dealii/"
    # print( "Dealii" )
    # showDealiiPlot( simPlotFolderName, regexVals, meshpath= meshPath, varName = "solutionvalues", negative = -1, pival = 2*pi )

    simPlotFolderName = simPlotRootFolderName + "Finch/"
    print( "Finch" )
    # buildAllMeshes( gmshFileCmdNames, meshPath )
    meshArr = getMeshFilesFromFolder( meshPath )
    # runFinchSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams )
    # showFinchPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath )

    srcFileName = "/home/gaurav/dealii-9.5.1/examples/step-5/step-5.cc"
    runDealiiSimWithOptionsVariousMeshes( optionsParam, meshArr, srcFileName )

    simPlotFolderName = simPlotRootFolderName + "Dealii/"
    showDealiiPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath, negative=-1, pival = pi )

    # fileName = folderUtils.textFolderNames["Dealii"] + "Mesh_solutionvalues_lvl=7.vtu"
    # getDealiiData( fileName, "vtu" )
    
    # showFinchDerivPlot( simPlotFolderName, ["regular"] )
    # showParaviewPlot( simPlotFolderName, regexVals, getMinMaxRangeFunc, meshPath )
    # showplot( simPlotFolderName, regexVals )
    # compareParaview( simPlotFolderName, regexVals, meshPath )

    # setFinchTriangleQuadrature( 2 )

