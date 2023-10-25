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
import fileNameUtils
import meshFileUtils
import regexUtils
import dataUtils
import errorUtils

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

def runFinchSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams, comparisonParam, meshPath ):
    
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
        regexCriterias = regexUtils.getRegexCriterias( regexVal )
        sortedMeshVals = meshFileUtils.getSortedMeshVals( meshvals, regexVal, regexCriterias )
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
        regexCriterias = regexUtils.getRegexCriterias( regexVal )

        for index, meshval in enumerate( sortedMeshVals ):

            numNodeVals.append( getNumNodes( meshpath + meshval ) )
            areaVals.append( getAverage2DArea( meshpath + meshval, regexVal ) )

            xvalsFileName = fileNameUtils.getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "centroid_xvalues",
                                                     regexVal, regexCriterias, ".txt" )
            yvalsFileName = fileNameUtils.getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "centroid_yvalues",
                                                    regexVal, regexCriterias, ".txt" )

            deriv_xvalsFileName = fileNameUtils.getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "deriv_xvalues",
                                                          regexVal, regexCriterias, ".txt" )
            deriv_yvalsFileName = fileNameUtils.getFileNameFromMeshName( meshval, folderUtils.finchTextfoldername, "deriv_yvalues",
                                                          regexVal, regexCriterias, ".txt" )
            # errvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )

            xvals = dataUtils.getData( xvalsFileName )
            yvals = dataUtils.getData( yvalsFileName )
            deriv_xvals = dataUtils.getData( deriv_xvalsFileName )
            deriv_yvals = dataUtils.getData( deriv_yvalsFileName )

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

    levelsArr = meshFileUtils.getAllLevels( meshvals )

    juliaVarName = "errorvalues"
    minMaxRangeVals = dataUtils.getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName )

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
    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )
        dxVals = []
        optionsParam[ comparisonParam ] = paramValue

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str( level ) ]
            criteriaValsStr = regexUtils.getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, optionsParam[ "meshRegexVal" ] ) )
            areaVals.append( getAverage2DArea( meshFileName, optionsParam[ "meshRegexVal" ] ) )

            xvalsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "xvalues" )
            yvalsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "yvalues" )
            uvalsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "uvalues" )
            
            xvals = dataUtils.getData( xvalsFileName )
            yvals = dataUtils.getData( yvalsFileName )
            uvals = dataUtils.getData( uvalsFileName )
            errvals = errorUtils.getFinchError( xvals, yvals, uvals, int( optionsParam[ "sin(kpix)" ] )*pi )

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
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            contourLevels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            plt.tricontourf( xvals, yvals, errvals, levels = contourLevels, colors = 'r')
            plt.colorbar()
            plotVarName = "large" + "errorvalues" + "Contour"
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, ".png" )
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

    levelsArr = meshFileUtils.getAllLevels( meshvals )

    # dealiiVarName = "errorvalues"
    # minMaxRangeVals = getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName )

    minMaxRangeVals = dataUtils.getDealiiMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, negative, pival )

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
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )
        dxVals = []
        optionsParam[ comparisonParam ] = paramValue

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str( level )]
            criteriaValsStr = regexUtils.getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, optionsParam[ "meshRegexVal" ] ) )
            areaVals.append( getAverage2DArea( meshFileName, optionsParam[ "meshRegexVal" ] ) )

            dealiiVarName = "solutionvalues"
            solvalsFileName = fileNameUtils.getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, dealiiVarName, "vtu" )
            
            (nodes, solution) = dataUtils.getDealiiData( solvalsFileName, "vtu" )
            xvals = nodes[:, 0]
            yvals = nodes[:, 1]
            errvals = errorUtils.getDealiiError( nodes, solution, negative, pival )
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
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )
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
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, ".png" )            
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

def compareDealiiFinch( simPlotFolderName, allParams, optionsParam, comparisonParam, meshvals, meshPath, negative = 1, pival = 2*pi ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    levelsArr = meshFileUtils.getAllLevels( meshvals )
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

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        dealiiCurMaxErrorList = []
        dealiiCurL2ErrorList = []

        regexVal = regexVals[idx]
        regexCriterias = regexUtils.getRegexCriterias( regexVal )

        finchCurMaxErrorList = []
        finchCurL2ErrorList = []
        numNodeVals = []
        areaVals = []

        regexCriterias = regexUtils.getRegexCriterias( "lvl" )
        dxVals = []
        optionsParam[ comparisonParam ] = paramValue

        figArrayNumNodes, axisArrayNumNodes = createMultipleFigureAxis( 2 )
        figArrayArea, axisArrayArea = createMultipleFigureAxis( 2 )
        figArrayDx, axisArrayDx = createMultipleFigureAxis( 2 )

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str( level )]
            criteriaValsStr = regexUtils.getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, optionsParam[ "meshRegexVal" ] ) )
            areaVals.append( getAverage2DArea( meshFileName, optionsParam[ "meshRegexVal" ] ) )

            dealiiVarName = "solutionvalues"
            solvalsFileName = fileNameUtils.getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, dealiiVarName, "vtu" )
            
            (nodes, solution) = dataUtils.getDealiiData( solvalsFileName, "vtu" )
            dealiiErrvals = errorUtils.getDealiiError( nodes, solution, negative = -1, pival = pival )
            
            # print(solvalsFileName)

            dealiiCurmaxval = np.max(dealiiErrvals)

            dealiiCurMaxErrorList.append( dealiiCurmaxval )
            dealiiCurL2ErrorList.append( np.sqrt( np.sum( np.array( dealiiErrvals )**2 ) / numNodeVals[ -1 ] ) ) 

            # Finch data
            finchErrvalsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "errorvalues" )
            
            finchErrvals = dataUtils.getData( finchErrvalsFileName )

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

        label1 = comparisonParam + "=" + paramValue + " Dealii "
        label2 = comparisonParam + "=" + paramValue + " Finch "

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

        fileNameLabel = comparisonParam + "=" + paramValue

        fileNameArray = [ curPlotFolderName + fileNameLabel + "_maxErrorNumNodes" + ".png",
                        curPlotFolderName + fileNameLabel + "_l2ErrorNumNodes" + ".png",
                        curPlotFolderName + fileNameLabel + "_maxErrorArea" + ".png",
                        curPlotFolderName + fileNameLabel + "_l2ErrorArea" + ".png",
                        curPlotFolderName + fileNameLabel + "_maxErrorDx" + ".png",
                        curPlotFolderName + fileNameLabel + "_l2ErrorDx" + ".png" ]

        figArray = [ figArrayNumNodes[0], figArrayNumNodes[1],
                    figArrayArea[0],
                    figArrayArea[1],
                    figArrayDx[0],
                    figArrayDx[1] ]

        makeMultiplePlotAdjustmentsAndSave( axisArray, figArray, fileNameArray )

    plt.show()

    return


if __name__ == "__main__":

    getMinMaxRangeFunc = dict()
    getMinMaxRangeFunc["Finch"] = dataUtils.getFinchMinMaxRange
    getMinMaxRangeFunc["Dealii"] = dataUtils.getDealiiMinMaxRange

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
    optionsParam["meshRegexVal"] = regexVals[1]
    optionsParam["level"] = "0"

    # pythonVarName = getPythonVarName( optionsParam )
    # print(pythonVarName)
    # filename = getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, "solutionValues", "vtu" )
    # print(filename) 

    comparisonParam = "quadratureOrder"

    simPlotRootFolderName = folderUtils.gmshImageFolderName + "PlotTriangleStructuredQuadrature_pi/"
    # meshPlotRootFolderName = folderUtils.gmshImageFolderName + "MeshPlotsHangingLevel_QuadratureOrder=2_pi/"

    # regexVals = [ "mesh" ]
    meshPath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    buildAllMeshes( gmshFileCmdNames, meshPath )
    meshArr = meshFileUtils.getMeshFilesFromFolder( meshPath )
    # showMeshes( folderUtils.meshPlotRootFolderName, regexVals )

    # createMeshVTU( meshPlotRootFolderName, regexVals )
 
    simPlotFolderName = simPlotRootFolderName + "Finch/"
    print( "Finch" )
    # runFinchSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams, comparisonParam, meshPath )
    showFinchPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath )

    srcFileName = "/home/gaurav/dealii-9.5.1/examples/step-5/step-5.cc"
    # runDealiiSimWithOptionsVariousMeshes( optionsParam, meshArr, allParams, comparisonParam, meshPath, srcFileName )

    simPlotFolderName = simPlotRootFolderName + "Dealii/"
    print( "dealii" )
    showDealiiPlot( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath, negative=-1, pival = pi )

    # compareDealiiFinch( simPlotRootFolderName, allParams,
                    # optionsParam, comparisonParam, meshArr, meshPath, negative =-1, pival = pi )

    # fileName = folderUtils.textFolderNames["Dealii"] + "Mesh_solutionvalues_lvl=7.vtu"
    # getDealiiData( fileName, "vtu" )
    
    # showFinchDerivPlot( simPlotFolderName, ["regular"] )
    # showParaviewPlot( simPlotFolderName, regexVals, getMinMaxRangeFunc, meshPath )
    # showplot( simPlotFolderName, regexVals )
    # compareParaview( simPlotFolderName, regexVals, meshPath )

    # setFinchTriangleQuadrature( 2 )
    # setFinchQuadQuadrature(3)
