import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import re
from os import listdir
from os.path import isfile, join
from math import sin, cos, pi
import os
import subprocess
from paraview.simple import *
import meshio
import pdb
import time
import h5py
paraview.simple._DisableFirstRenderCameraReset()

finchRootfoldername = "/media/gaurav/easystore/Finch/MixedElement/"
dealiiRootfoldername = "/media/gaurav/easystore/dealii/MixedElement/"

finchTextfoldername = finchRootfoldername + "TextFiles/"
finchPlotfoldername = finchRootfoldername + "PlotFiles/SimPlots/"

dealiiTextfoldername = dealiiRootfoldername + "TextFiles/"
dealiiPlotfoldername = dealiiRootfoldername + "PlotFiles/SimPlots/"
gmshImageFolderName = "/home/gaurav/gmshAutoScripts/Images/"

textFolderNames = dict()
textFolderNames["Finch"] = finchTextfoldername
textFolderNames["Dealii"] = dealiiTextfoldername

plotFolderNames = dict()
plotFolderNames["Finch"] = finchPlotfoldername
plotFolderNames["Dealii"] = dealiiPlotfoldername

def checkAndCreateFolder( folderName ):

    if not os.path.exists( folderName ):

        folderSplit = folderName.split( "/" )

        folderStr = ""

        for folderName in folderSplit:

            if not folderName:
                continue
            
            folderStr = folderStr + "/" + folderName

            if not os.path.exists( folderStr ):
                os.mkdir( folderStr )        


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

def getFinchMinMaxRange( sortedMeshValsArr, regexVals, varName = "errorvalues" ):

    meshValsLen = len( sortedMeshValsArr[0] )
    minVals = np.ones( meshValsLen )*100000
    maxVals = np.ones( meshValsLen )*-1

    for typeIdx, sortedMeshVals in enumerate( sortedMeshValsArr ):

        regexVal = regexVals[ typeIdx ]
        regexCriterias = getRegexCriterias( regexVal )
        for idx, meshval in enumerate( sortedMeshVals ):

            textFile = getFileNameFromMeshName( meshval, finchTextfoldername, varName, regexVal, regexCriterias, ".txt" )

            errVals = getData( textFile )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )

    minMaxRangeVals = []

    for idx in range(meshValsLen):
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


def getDealiiMinMaxRange( sortedMeshValsArr, regexVals, varName = "solutionvaluesGaussModified", format = "h5", negative = 1, pival = 2*pi ):

    meshValsLen = len( sortedMeshValsArr[0] )
    minVals = np.ones( meshValsLen )*100000
    maxVals = np.ones( meshValsLen )*-1

    for typeIdx, sortedMeshVals in enumerate( sortedMeshValsArr ):

        regexVal = regexVals[typeIdx]
        regexCriterias = getRegexCriterias( regexVal )
        for idx, meshval in enumerate( sortedMeshVals ):

            h5FileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, varName,
                                                regexVal, regexCriterias, "." + format )
            
            nodes, solution = getDealiiData( h5FileName, format )
            errVals = getDealiiError( nodes, solution, negative, pival )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )

    minMaxRangeVals = []

    for idx in range(meshValsLen):
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

    statement = "xyw = triangle_quadrature_nodes_weights(" +  quadratureValue + ");"
    filename = "/home/gaurav/.julia/packages/Finch/ECEMc/src/triangle_nodes.jl"
    pattern = "triangle_quadrature_nodes_weights"

    with open( filename ) as fileHandle:

        allLines = fileHandle.readlines()

        for idx, lineVal in enumerate( allLines ):

            if re.search( pattern, lineVal ):

                allLines[idx] = statement
                break

    with open( filename ) as fileHandle:

        fileHandle.writelines( allLines )

def runFinchSimWithOptions( optionsParam ):

    exefilename = "example-mixed-element-2d.jl"
    setFinchTriangleQuadrature( optionsParam[ "quadratureOrder" ] )
    pythonVarName = getPythonVarName( optionsParam )
    allFinchOptions = [ optionsParam["sin(k*pi*x)"], optionsParam["coeff_F"], pythonVarName ]

    runJulia( exefilename, allFinchOptions )

def runJulia( exefilename, options = [] ):

    # juliapath = "/home/gaurav/Downloads/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    juliapath = "/home/gaurav/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    finchPath = "/home/gaurav/Finch/src/examples/"
    julialstcmd = [juliapath, exefilename]

    for option in options:
        julialstcmd.append( option )

    subprocess.run( julialstcmd, cwd = finchPath )

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

def createMeshVTU( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )
    
    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:

        regexCriterias = getRegexCriterias( regexVal )

        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

def showMeshes( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )
    
    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:

        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    for typeIdx, sortedMeshVals in enumerate( allSortedMeshVals ):
        
        regexVal = regexVals[typeIdx]
        regexCriterias = getRegexCriterias( regexVal )

        for index, meshval in enumerate(sortedMeshVals):

            meshvtkName = meshVizPath + meshval[:-4] + ".vtu"

            gmshfile = OpenDataFile( meshvtkName )
            dpGmsh = GetDisplayProperties( gmshfile )
            dpGmsh.Representation = 'Wireframe'
            gmshdisplay = Show(gmshfile)
            
            myview = GetActiveView()
            myview.ViewSize = [1920, 1080]
            myview.InteractionMode = '2D'
            myview.AxesGrid = 'Grid Axes 3D Actor'
            myview.CenterOfRotation = [0.5, 0.5, 0.0]
            myview.StereoType = 'Crystal Eyes'
            myview.CameraPosition = [0.5, 0.5, 3.0349403797187358]
            myview.CameraFocalPoint = [0.5, 0.5, 0.0]
            myview.CameraFocalDisk = 1.0
            myview.CameraParallelScale = 0.7908298380174797
            myview.LegendGrid = 'Legend Grid Actor'

            Render()

            curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"  
            checkAndCreateFolder( curPlotFolderName )
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, 
                                                   "paraview_error", regexVal, 
                                                   regexCriterias, ".png" )
            SaveScreenshot( plotfilename, myview)

            Hide( gmshfile )

def showParaviewPlot( simPlotFolderName, regexVals, getMinMaxRangeFunc, 
                     meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/",
                      meshVizPath =  "/home/gaurav/Finch/src/examples/Mesh/MeshViz/", 
                      software = "Finch" ):

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getMinMaxRangeFunc[software]( allSortedMeshVals, regexVals )

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        regexVal = regexVals[idx]
        regexCriterias = getRegexCriterias( regexVal )

        for index, meshval in enumerate(sortedMeshVals):

            meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
            hangingfilename = getFileNameFromMeshName( meshval, textFolderNames[ software ],
                                                    "errorAndu", regexVal, regexCriterias, ".vtu" )

            solfile = OpenDataFile( hangingfilename )
            display = Show(solfile)
            ColorBy(display, ('POINTS', 'err'))
            minVal, maxVal = minMaxRangeVals[ index ]
            colorMap = GetColorTransferFunction('err')
            colorMap.RescaleTransferFunction( minVal, maxVal )        

            gmshfile = OpenDataFile( meshvtkName )
            dpGmsh = GetDisplayProperties( gmshfile )
            dpGmsh.Representation = 'Wireframe'
            gmshdisplay = Show(gmshfile)
            
            myview = GetActiveView()
            myview.ViewSize = [1920, 1080]
            myview.InteractionMode = '2D'
            # myview.AxesGrid = 'Grid Axes 3D Actor'
            myview.CenterOfRotation = [0.5, 0.5, 0.0]
            myview.StereoType = 'Crystal Eyes'
            myview.CameraPosition = [0.5, 0.5, 3.0349403797187358]
            myview.CameraFocalPoint = [0.5, 0.5, 0.0]
            myview.CameraFocalDisk = 1.0
            myview.CameraParallelScale = 0.7908298380174797
            # myview.LegendGrid = 'Legend Grid Actor'

            Render()
            
            dpSol = GetDisplayProperties(solfile, myview)
            # to show the color legend
            dpSol.SetScalarBarVisibility(myview, True)

            curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"   
            if not os.path.exists( curPlotFolderName ):
                os.mkdir( curPlotFolderName )

            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "paraview_error",
                                                   regexVal, regexCriterias, ".png" )
            SaveScreenshot( plotfilename, myview)

            Hide( solfile )
            Hide( gmshfile )

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

    checkAndCreateFolder( simPlotFolderName )

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

            xvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "centroid_xvalues",
                                                     regexVal, regexCriterias, ".txt" )
            yvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "centroid_yvalues",
                                                    regexVal, regexCriterias, ".txt" )

            deriv_xvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "deriv_xvalues",
                                                          regexVal, regexCriterias, ".txt" )
            deriv_yvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "deriv_yvalues",
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

def showFinchPlot( simPlotFolderName, regexVals, 
                meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/", 
                varName = "errorvalues" ):

    checkAndCreateFolder( simPlotFolderName )

    import matplotlib.ticker as ticker

    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getFinchMinMaxRange( allSortedMeshVals, regexVals, varName )

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

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexVal = regexVals[idx]
        regexCriterias = getRegexCriterias( regexVal )
        dxVals = []

        for index, meshval in enumerate( sortedMeshVals ):

            meshFileName = meshPath + meshval

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, regexVal ) )
            areaVals.append( getAverage2DArea( meshFileName, regexVal ) )

            xvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "xvalues", 
                                                    regexVal, regexCriterias, ".txt" )
            yvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "yvalues", 
                                                    regexVal, regexCriterias, ".txt" )
            # errvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )
            uvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "uvalues",
                                                    regexVal, regexCriterias, ".txt" )

            xvals = getData( xvalsFileName )
            yvals = getData( yvalsFileName )
            uvals = getData( uvalsFileName )
            errvals = getFinchError( xvals, yvals, uvals )

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

            levels = np.linspace( minVal, maxVal, 40 )
            plt.figure()
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = levelsmin[-1], vmax = levelsmax[-1] )
            plt.tricontourf( xvals, yvals, errvals, levels = levels, vmin = minVal, vmax = maxVal, cmap = cm )
            plt.colorbar()
            # plt.scatter( xvalsHanging, yvalsHanging )
            # plt.show()
            fileVarName = varName + "Contour"
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, fileVarName, 
                                                   regexVal, regexCriterias, ".png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            levels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, errvals, levels = levels, colors = 'r')
            plt.colorbar()
            fileVarName = "large" + varName + "Contour"
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, fileVarName,
                                                    regexVal, regexCriterias, ".png" )
            plt.savefig( plotfilename )
            plt.close()
            # plt.figure()
            # plt.tricontourf( xvals, yvals, uvals, levels = 20 )
            # plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            # plt.savefig( plotfoldername + "ucontour_index=" + str(index) + ".png" )

            # plt.figure()
            # plt.tricontourf( xvals, yvals, uexactvals, levels = 20 )
            # plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            # plt.savefig( plotfoldername + "uexactcontour_index=" + str(index) + ".png" )

        allMaxErrorList.append( curMaxErrorList )
        allL2ErrorList.append( curL2ErrorList )

        plotMaxAndL2Error( axisArrayNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, regexVal )
        plotMaxAndL2Error( axisArrayArea, areaVals, curMaxErrorList, curL2ErrorList, regexVal ) 
        plotMaxAndL2Error( axisArrayDx, dxVals, curMaxErrorList, curL2ErrorList, regexVal )      

        h2Vals = [ dxVal**2 for dxVal in dxVals ]

        print( regexVal )
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

    if len( regexVals ) > 1:
        errorDiffMax = [ allMaxErrorList[1][idx] - allMaxErrorList[0][idx] for idx in range( len(allMaxErrorList[0]) ) ]
        errorDiffL2 = [ allL2ErrorList[1][idx] - allL2ErrorList[0][idx] for idx in range( len(allL2ErrorList[0]) ) ]

        plt.figure()
        plt.plot( dxVals, errorDiffMax, "-o" )
        titleval = regexVals[1] + " - " + regexVals[0] + " Max Error Plot"
        plt.xlabel( "h" )
        plt.ylabel( "Error" )
        plt.title( titleval )

        plt.figure()
        titleval = regexVals[1] + " - " + regexVals[0] + " L2 Error Plot"
        plt.plot( dxVals, errorDiffL2, "-o" )
        plt.xlabel( "h" )
        plt.ylabel( "Error" )
        plt.title( titleval )

    plt.show()

    return

def plotMultipleAxis( axisArray, xvalsArray, yvalsArray, labelNameArray, markerVal ):

    for i, ax in enumerate( axisArray ):

        ax.loglog( xvalsArray[i], yvalsArray[i], markerVal, label = labelNameArray[i] )

def plotMaxAndL2Error( axisArray, xvals, maxError, l2Error, regexVal ):

    labelNameArray = [ regexVal + " $L^{\infty}$ Error_" + "h", regexVal + " $L^{2}$ Error_" + "h" ]
    xvalsArray = [ xvals, xvals ]
    yvalsArray = [ maxError, l2Error ]

    plotMultipleAxis( axisArray, xvalsArray, yvalsArray, labelNameArray, "-o" )

    # labelName = regexVal + " $L^{\infty}$ Error_" + "2h"
    # axMaxErrorHandle.loglog( xvals[:-1], maxError[:-1], "-o", label = labelName )

    # labelName = regexVal + " $L^{\infty}$ Error_" + "h"
    # axMaxErrorHandle.loglog( xvals[:-1], maxError[1:], "-o", label = labelName )

    # labelName = regexVal + " $L^{2}$ Error_" + "2h"
    # axL2ErrorHandle.loglog( xvals[:-1], l2Error[:-1], "-o", label = labelName )

    # labelName = regexVal + " $L^{2}$ Error_" + "h"
    # axL2ErrorHandle.loglog( xvals[:-1], l2Error[1:], "-o", label = labelName )

def showDealiiPlot( simPlotFolderName, regexVals, varName, 
                   meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/", 
                   negative = 1, pival = 2*pi ):

    checkAndCreateFolder( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfiletriangleunstruct.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getDealiiMinMaxRange( allSortedMeshVals, regexVals, varName,
                                        format = "vtu", negative = negative, pival =  pival )

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

    print(allSortedMeshVals)

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        dxVals = []
        regexVal = regexVals[idx]
        regexCriterias = getRegexCriterias( regexVal )

        for index, meshval in enumerate( sortedMeshVals ):

            meshFileName = meshPath + meshval

            numNodeVals.append( getNumNodes( meshFileName ) )
            dxVals.append( getMaxH( meshFileName, regexVal ) )
            areaVals.append( getAverage2DArea( meshFileName, regexVal ) )

            solvalsFileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, varName,
                                                      regexVal, regexCriterias, ".vtu" )
            
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

            checkAndCreateFolder( curPlotFolderName )

            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
            # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
            # print( np.max(errvals) )

            minVal, maxVal = minMaxRangeVals[ index ]

            levels = np.linspace( minVal, maxVal, 40 )
            plt.figure()
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = levelsmin[-1], vmax = levelsmax[-1] )
            plt.tricontourf( xvals, yvals, errvals, levels = levels, vmin = minVal, vmax = maxVal, cmap = cm )
            plt.colorbar()
            # plt.scatter( xvalsHanging, yvalsHanging )
            # plt.show()
            fileVarName = varName + "Contour"
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, fileVarName, regexVal, regexCriterias, ".png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            levels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, errvals, levels = levels, colors = 'r')
            plt.colorbar()
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "largeErrorContour_", regexVal, regexCriterias, ".png" )
            plt.savefig( plotfilename )
            plt.close()
            # plt.figure()
            # plt.tricontourf( xvals, yvals, uvals, levels = 20 )
            # plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            # plt.savefig( plotfoldername + "ucontour_index=" + str(index) + ".png" )

            # plt.figure()
            # plt.tricontourf( xvals, yvals, uexactvals, levels = 20 )
            # plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            # plt.savefig( plotfoldername + "uexactcontour_index=" + str(index) + ".png" )

        allMaxErrorList.append( curMaxErrorList )
        allL2ErrorList.append( curL2ErrorList )

        plotMaxAndL2Error( axisArrayNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, regexVal )
        plotMaxAndL2Error( axisArrayArea, areaVals, curMaxErrorList, curL2ErrorList, regexVal ) 
        plotMaxAndL2Error( axisArrayDx, dxVals, curMaxErrorList, curL2ErrorList, regexVal )

        print( regexVal )
        print( dxVals )
        print( curMaxErrorList )
        print( curL2ErrorList )

    h2Vals = [ dxVal**2 for dxVal in dxVals ]

    # print( areaVals )
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

    if len( regexVals ) > 1:
        errorDiffMax = [ allMaxErrorList[1][idx] - allMaxErrorList[0][idx] for idx in range( len(allMaxErrorList[0]) ) ]
        errorDiffL2 = [ allL2ErrorList[1][idx] - allL2ErrorList[0][idx] for idx in range( len(allL2ErrorList[0]) ) ]

        plt.figure()
        plt.plot( dxvals, errorDiffMax, "-o" )
        titleval = regexVals[1] + " - " + regexVals[0] + " Max Error Plot"
        plt.xlabel( "h" )
        plt.ylabel( "Error" )
        plt.title( titleval )

        plt.figure()
        titleval = regexVals[1] + " - " + regexVals[0] + " L2 Error Plot"
        plt.plot( dxvals, errorDiffL2, "-o" )
        plt.xlabel( "h" )
        plt.ylabel( "Error" )
        plt.title( titleval )

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

    checkAndCreateFolder( simPlotFolderName )

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

            dealiiSolvalsFileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, varName,
                                                            regexVal, regexCriterias,  ".vtu" )
            
            (nodes, solution) = getDealiiData( dealiiSolvalsFileName, format = "vtu" )
            dealiiErrvals = getDealiiError( nodes, solution, negative=-1 )
            # print(solvalsFileName)

            dealiiCurmaxval = np.max(dealiiErrvals)

            dealiiCurMaxErrorList.append( dealiiCurmaxval )
            dealiiCurL2ErrorList.append( np.sqrt( np.sum( np.array( dealiiErrvals )**2 ) / numNodeVals[ -1 ] ) ) 

            # Finch data
            finchXvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "xvalues",
                                                         regexVal, regexCriterias, ".txt" )
            finchYvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "yvalues",
                                                         regexVal, regexCriterias, ".txt" )
            finchErrvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues",
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
        checkAndCreateFolder( curPlotFolderName )

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

def compareParaview( simPlotFolderName, regexVals, getMinMaxRangeFunc,
                    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/", 
                     software = "Finch" ):

    assert( len(regexVals) == 2 )

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        regexCriterias = getRegexCriterias( regexVal )
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getMinMaxRangeFunc[software]( allSortedMeshVals )

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    meshValsLen = len( allSortedMeshVals[0] )

    allViews = dict()
    for regexVal in regexVals:
        allViews[regexVal] = CreateRenderView()
        renderView = allViews[regexVal]
        renderView.ViewSize = [701, 784]
        renderView.InteractionMode = '2D'
        renderView.AxesGrid = 'Grid Axes 3D Actor'
        renderView.CenterOfRotation = [0.5, 0.5, 0.0]
        renderView.StereoType = 'Crystal Eyes'
        renderView.CameraPosition = [0.5, 0.5, 3.0349403797187358]
        renderView.CameraFocalPoint = [0.5, 0.5, 0.0]
        renderView.CameraFocalDisk = 1.0
        renderView.CameraParallelScale = 0.7908298380174797
        renderView.LegendGrid = 'Legend Grid Actor'

    layout = CreateLayout(name='Layout #1')
    layout.SplitHorizontal(0, 0.500000)

    for idx, regexVal in enumerate(regexVals):
        layout.AssignView( idx + 1, allViews[regexVal] )

    layout.SetSize(1403, 784)
    # layout.SetSize(1920, 1080)

    for idx in range(meshValsLen):

        curPlotFolderName = simPlotFolderName + "Plot" + str(idx) + "/"   
        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        solfiles = []
        gmshfiles = []

        for regexIdx, regexVal in enumerate( regexVals ):
            
            SetActiveView( allViews[regexVal] )
            myview = GetActiveView()

            regexCriterias = getRegexCriterias[ regexVal ]

            meshval = allSortedMeshVals[regexIdx][idx]
            meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
            filename = getFileNameFromMeshName( meshval, textFolderNames[software],
                                             "errorAndu", regexVal, regexCriterias, ".vtu" )

            solfile = OpenDataFile( filename )
            solfiles.append(solfile)
            display = Show(solfile)
            ColorBy(display, ('POINTS', 'err'))
            minVal, maxVal = minMaxRangeVals[ idx ]
            colorMap = GetColorTransferFunction('err')
            colorMap.RescaleTransferFunction( minVal, maxVal )        

            gmshfile = OpenDataFile( meshvtkName )
            gmshfiles.append(gmshfile)
            dpGmsh = GetDisplayProperties( gmshfile )
            dpGmsh.Representation = 'Wireframe'
            gmshdisplay = Show(gmshfile)
            
            dpSol = GetDisplayProperties(solfile, myview)
            # # to show the color legend
            dpSol.SetScalarBarVisibility(myview, True)
            myview.Update()

        Render()
        plotfilename = curPlotFolderName + "paraview_error_comparison.png" 
        SaveScreenshot( plotfilename, layout)        

        for regexIdx, regexVal in enumerate( regexVals ):

            SetActiveView( allViews[regexVal] )
            Hide( solfiles[regexIdx] )
            Hide( gmshfiles[regexIdx] )

    return

def getPythonVarName( optionsParam ):

    pythonVarName = ""

    for param, paramValue in optionsParam.items():

        pythonVarName += param + "_" + paramValue

    return pythonVarName

def getFileNameFromMeshName( folderName, meshFileName, pythonVarName, juliaVarName, criteria = "lvl", extension = "txt" ):

    regexCriterias = getRegexCriterias( criteria )
    criteriaVals = getCriteriasFromFileName( meshFileName, regexCriterias )

    filename = getTextFileName( folderName, pythonVarName, juliaVarName, criteriaVals, criteria, extension )

    return filename

def getTextFileName( folderName, pythonVarName, juliaVarName, criteriaVals, criteria = "lvl", extension = "txt" ):

    filename = folderName + pythonVarName + "_" + juliaVarName + "_" + \
        criteria + "=" + criteriaVals + "." + extension

    return filename

def getCriteriasFromFileName( meshval, regexCriterias ):

    regexCriteriaVals = ""

    for regexCriteria in regexCriterias:
        rval = re.search( regexCriteria, meshval )

        if rval:
            offset = rval.start() + len( regexCriteria )
            rval = re.search( "[0-9]+", meshval[offset:] )
            rval = rval.group(0)

            regexCriteriaVals += rval + "_"

    return regexCriteriaVals[:-1]

if __name__ == "__main__":

    getMinMaxRangeFunc = dict()
    getMinMaxRangeFunc["Finch"] = getFinchMinMaxRange
    getMinMaxRangeFunc["Dealii"] = getDealiiMinMaxRange

    # gmshFileCmdNames = ["triangleMeshv1", "triangleMeshv2"]
    # regexVals = ["triangleMeshStruct", "triangleMeshUnstruct"]
    regexVals = ["triangleMeshUnstruct", "triangleMeshStruct", "regular"]
    gmshFileCmdNames = ["triangleMeshv2", "triangleMeshv1", "regularMeshv3"]
    # gmshFileCmdNames = ["hangingMeshv8"]
    # regexVals = ["hanging"]

    allParams = dict()
    allParams["softwares"] = ["Finch", "Dealii"]
    allParams["meshRegexVals"] = regexVals
    allParams["gmshFileCmdNames"] = gmshFileCmdNames
    allParams["quadratureOrders"] = [ 2, 3, 4, 5, 6, 7 ]
    allParams["sin(k*pi*x)s"] = [ 1, 2, 4 ]
    allParams["coeff_Fs"] = [ -2, -8, -32 ]

    optionsParam = dict()
    optionsParam["quadratureOrder"] = 2
    optionsParam["sin(k*pi*x)"] = 1
    optionsParam["coeff_F"] = -2
    optionsParam["software"] = "Finch"
    optionsParam["meshRegexVal"] = regexVals[0]

    comparisonParam = "meshRegexVals"

    simPlotRootFolderName = gmshImageFolderName + "PlotRegularMeshFinchTriPoints3_2pi/"
    meshPlotRootFolderName = gmshImageFolderName + "MeshPlotsHangingLevel_QuadratureOrder=2_2pi/"

    # regexVals = [ "mesh" ]
    meshPath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    # showMeshes( meshPlotRootFolderName, regexVals )

    # runSim( simPlotRootFolderName, gmshFileCmdNames, regexVals )
    # createMeshVTU( meshPlotRootFolderName, regexVals )
    # compareDealiiFinch( simPlotRootFolderName, regexVals, "solutionvalues", meshPath )

    simPlotFolderName = simPlotRootFolderName + "Dealii/"
    print( "Dealii" )
    # showDealiiPlot( simPlotFolderName, regexVals, meshpath= meshPath, varName = "solutionvalues", negative = -1, pival = 2*pi )

    simPlotFolderName = simPlotRootFolderName + "Finch/"
    print( "Finch" )
    showFinchPlot( simPlotFolderName, regexVals, meshPath )

    # fileName = textFolderNames["Dealii"] + "Mesh_solutionvalues_lvl=7.vtu"
    # getDealiiData( fileName, "vtu" )
    
    # showFinchDerivPlot( simPlotFolderName, ["regular"] )
    # showParaviewPlot( simPlotFolderName, regexVals, getMinMaxRangeFunc, meshPath )
    # showplot( simPlotFolderName, regexVals )
    # compareParaview( simPlotFolderName, regexVals, meshPath )

