import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import re
from os import listdir
from os.path import isfile, join
from math import sin, cos, pi
import os
import subprocess
# from paraview.simple import *
import meshio
import pdb
import time
import h5py
# paraview.simple._DisableFirstRenderCameraReset()

finchRootfoldername = "/media/gaurav/easystore/Finch/MixedElement/"
dealiiRootfoldername = "/media/gaurav/easystore/dealii/MixedElement/"

finchTextfoldername = finchRootfoldername + "TextFiles/"
finchPlotfoldername = finchRootfoldername + "PlotFiles/SimPlots/"

dealiiTextfoldername = dealiiRootfoldername + "TextFiles/"
dealiiPlotfoldername = dealiiRootfoldername + "PlotFiles/SimPlots/"
gmshImageFolderName = "/home/gaurav/gmshAutoScripts/Images/"

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


def getSortedMeshVals( meshvals, regexVal ):

    keyvals = []
    prevIndexMap = dict()

    if regexVal == "hanging":

        for index, meshval in enumerate(meshvals):

            if re.search( regexVal, meshval  ):
                Nxval = re.search( "Nx=", meshval )
                offset = Nxval.start() + 3
                Nxval = re.search( "[0-9]+", meshval[offset:] )
                Nxval = int( Nxval.group(0) )

                Nyval = re.search( "Ny=", meshval )
                offset = Nyval.start() + 3
                Nyval = re.search( "[0-9]+", meshval[offset:] )
                Nyval = int( Nyval.group(0) )

                keyval = Nxval*Nyval
                keyvals.append( keyval )

                prevIndexMap[ keyval ] = index

    if regexVal == "regular" or regexVal == "triangle" or regexVal == "triangleMeshUnstruct" or regexVal == "triangleMeshStruct":

        for index, meshval in enumerate(meshvals):

            if re.search( regexVal, meshval ):

                regularNval = re.search( "[0-9]+"  , meshval  )
                regularNval = regularNval.group(0)
                
                keyval = int(regularNval)
                keyvals.append( keyval )

                prevIndexMap[ keyval ] = index

    keyvals = sorted( keyvals )
    sortedMeshVals = []

    for keyval in keyvals:

        idx = prevIndexMap[keyval]
        sortedMeshVals.append( meshvals[idx] )

    return sortedMeshVals

def getFileNameFromMeshName( meshval, folderName, varName, extension ):

    if re.search( "hanging", meshval  ):
        Nxval = re.search( "Nx=", meshval )
        offset = Nxval.start() + 3
        Nxval = re.search( "[0-9]+", meshval[offset:] )
        Nxval = Nxval.group(0)

        Nyval = re.search( "Ny=", meshval )
        offset = Nyval.start() + 3
        Nyval = re.search( "[0-9]+", meshval[offset:] )
        Nyval = Nyval.group(0)

        fileName = folderName + "mixed_" + varName + "Nx=" + Nxval + "Ny=" + Nyval + extension

    if re.search( "regular", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        fileName = folderName + "regular_" + varName + "N=" + regularNval + extension

    if re.search( "triangle", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        fileName = folderName + "triangle_" + varName + "N=" + regularNval + extension

    if re.search( "triangleMeshStruct", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        fileName = folderName + "triangleMeshStruct_" + varName + "N=" + regularNval + extension

    if re.search( "triangleMeshUnstruct", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        fileName = folderName + "triangleMeshUnstruct_" + varName + "N=" + regularNval + extension

    return fileName

def getFinchMinMaxRange( sortedMeshValsArr ):

    meshValsLen = len( sortedMeshValsArr[0] )
    minVals = np.ones( meshValsLen )*100000
    maxVals = np.ones( meshValsLen )*-1

    for sortedMeshVals in sortedMeshValsArr:
        for idx, meshval in enumerate( sortedMeshVals ):

            textFile = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )

            errVals = getData( textFile )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )

    minMaxRangeVals = []

    for idx in range(meshValsLen):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals

def getDealiiData( h5FileName ):

    f = h5py.File( h5FileName, 'r')
    nodes = np.array( f["nodes"] )
    solution = np.array( f["solution"] )

    return (nodes, solution)

def getExactSol( xval, yval ):

    return -sin(2*pi*xval)*sin(2*pi*yval)

def getDealiiError( nodes, solution ):

    errVals = []
    exactvals = []

    for idx, node in enumerate( nodes ):

        xval = node[0]
        yval = node[1]
        # print(xval, yval)
        exactval = getExactSol( xval, yval )
        exactvals.append( exactval )
        solval = solution[ idx, 0 ]

        errorval = np.abs( exactval - solval )
        errVals.append( errorval )

    # print(errVals)

    exactvals = np.array(exactvals)

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], exactvals )

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], solution[:, 0] )

    return np.array( errVals )


def getDealiiMinMaxRange( sortedMeshValsArr ):

    meshValsLen = len( sortedMeshValsArr[0] )
    minVals = np.ones( meshValsLen )*100000
    maxVals = np.ones( meshValsLen )*-1

    for sortedMeshVals in sortedMeshValsArr:
        for idx, meshval in enumerate( sortedMeshVals ):

            h5FileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, "solutionvaluesGaussModified_", ".h5" )
            nodes, solution = getDealiiData( h5FileName )
            errVals = getDealiiError( nodes, solution )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )

    minMaxRangeVals = []

    for idx in range(meshValsLen):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals


def getData( filename ):

    uvals = []

    with open(filename) as fval:

        uvals = fval.readlines()

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

def runJulia( exefilename ):

    # juliapath = "/home/gaurav/Downloads/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    juliapath = "/home/gaurav/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    finchPath = "/home/gaurav/Finch/src/examples/"
    julialstcmd = [juliapath, exefilename]

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

def runSim( simPlotFolderName, gmshFileCmdNames, regexVals ):

    gmshfileArgs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    removeFiles( gmshfileArgs )

    buildAllMeshes( gmshFileCmdNames, gmshfileArgs )
    runFinchSim() 
    runDealiiSim()   

    # showFinchPlot( simPlotFolderName, regexVals )
    # showParaviewPlot( simPlotFolderName, regexVals )
    # removeFiles( gmshfileargs )

def createMeshVTU( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )
    
    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
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
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
        allSortedMeshVals.append( sortedMeshVals )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    for sortedMeshVals in allSortedMeshVals:
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
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "paraview_error", ".png" )
            SaveScreenshot( plotfilename, myview)

            Hide( gmshfile )

# def showParaviewPlot( simPlotFolderName, regexVals ):

#     if not os.path.exists(simPlotFolderName):
#         os.mkdir( simPlotFolderName )

#     meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
#     meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"
#     meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

#     allSortedMeshVals = []

#     for regexVal in regexVals:
#         sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
#         allSortedMeshVals.append( sortedMeshVals )

#     minMaxRangeVals = getMinMaxRange( allSortedMeshVals )

#     for sortedMeshVals in allSortedMeshVals:
#         for index, meshval in enumerate( sortedMeshVals ):
#             mesh = meshio.read( meshpath + meshval )
#             meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

#     for sortedMeshVals in allSortedMeshVals:
#         for index, meshval in enumerate(sortedMeshVals):

#             meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
#             hangingfilename = getFileNameFromMeshName( meshval, textfoldername, "errorAndu", ".vtu" )

#             solfile = OpenDataFile( hangingfilename )
#             display = Show(solfile)
#             ColorBy(display, ('POINTS', 'err'))
#             minVal, maxVal = minMaxRangeVals[ index ]
#             colorMap = GetColorTransferFunction('err')
#             colorMap.RescaleTransferFunction( minVal, maxVal )        

#             gmshfile = OpenDataFile( meshvtkName )
#             dpGmsh = GetDisplayProperties( gmshfile )
#             dpGmsh.Representation = 'Wireframe'
#             gmshdisplay = Show(gmshfile)
            
#             myview = GetActiveView()
#             myview.ViewSize = [1920, 1080]
#             myview.InteractionMode = '2D'
#             myview.AxesGrid = 'Grid Axes 3D Actor'
#             myview.CenterOfRotation = [0.5, 0.5, 0.0]
#             myview.StereoType = 'Crystal Eyes'
#             myview.CameraPosition = [0.5, 0.5, 3.0349403797187358]
#             myview.CameraFocalPoint = [0.5, 0.5, 0.0]
#             myview.CameraFocalDisk = 1.0
#             myview.CameraParallelScale = 0.7908298380174797
#             myview.LegendGrid = 'Legend Grid Actor'

#             Render()
            
#             dpSol = GetDisplayProperties(solfile, myview)
#             # to show the color legend
#             dpSol.SetScalarBarVisibility(myview, True)

#             curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"   
#             if not os.path.exists( curPlotFolderName ):
#                 os.mkdir( curPlotFolderName )

#             plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "paraview_error", ".png" )
#             SaveScreenshot( plotfilename, myview)

#             Hide( solfile )
#             Hide( gmshfile )

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

def showFinchPlot( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfiletriangleunstruct.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getFinchMinMaxRange( allSortedMeshVals )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    allMaxErrorList = []
    allL2ErrorList = []

    figMaxErrorNumNodes = plt.figure()
    axMaxErrorNumNodes = figMaxErrorNumNodes.add_subplot(1, 1, 1)
    figL2ErrorNumNodes = plt.figure()
    axL2ErrorNumNodes = figL2ErrorNumNodes.add_subplot(1, 1, 1)

    figMaxErrorArea = plt.figure()
    axMaxErrorArea = figMaxErrorArea.add_subplot(1, 1, 1)
    figL2ErrorArea = plt.figure()
    axL2ErrorArea = figL2ErrorArea.add_subplot(1, 1, 1)

    # print(allSortedMeshVals)

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexVal = regexVals[idx]

        for index, meshval in enumerate( sortedMeshVals ):

            numNodeVals.append( getNumNodes( meshpath + meshval ) )
            areaVals.append( getAverage2DArea( meshpath + meshval, regexVal ) )

            xvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "xvalues_", ".txt" )
            yvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "yvalues_", ".txt" )
            errvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )

            xvals = getData( xvalsFileName )
            yvals = getData( yvalsFileName )
            errvals = getData( errvalsFileName )

            # print(errvalsFileName)

            curmaxval = np.max(errvals)
            curminval = np.min(errvals)

            curMaxErrorList.append( curmaxval )
            curL2ErrorList.append( np.sum( np.array(errvals)**2 ) ) 

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
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "errorContour_", ".png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            levels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, errvals, levels = levels, colors = 'r')
            plt.colorbar()
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "largeErrorContour_", ".png" )
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

        labelName = regexVals[idx] + " $L^{\infty}$ Error_" + "2h"
        plotMaxAndL2Error( axMaxErrorNumNodes, axL2ErrorNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, regexVal )
        plotMaxAndL2Error( axMaxErrorArea, axL2ErrorArea, areaVals, curMaxErrorList, curL2ErrorList, regexVal )      

    h2vals = dxvals**2

    axMaxErrorNumNodes.loglog( numNodeVals[:-1], h2vals[:-1], "-x", label = "$h^2$" )
    axL2ErrorNumNodes.loglog( numNodeVals[:-1], h2vals[:-1], "-x", label = "$h^2$" )

    axMaxErrorArea.loglog( areaVals[:-1], h2vals[:-1], "-x", label = "$h^2$" )
    axL2ErrorArea.loglog( areaVals[:-1], h2vals[:-1], "-x", label = "$h^2$" )

    makePlotAdjustmentsAndSave( axMaxErrorNumNodes, figMaxErrorNumNodes, simPlotFolderName + "maxErrorNumNodes" + ".png" )
    makePlotAdjustmentsAndSave( axL2ErrorNumNodes, figL2ErrorNumNodes, simPlotFolderName + "l2ErrorNumNodes" + ".png" )

    makePlotAdjustmentsAndSave( axMaxErrorArea, figMaxErrorArea, simPlotFolderName + "maxErrorArea" + ".png" )
    makePlotAdjustmentsAndSave( axL2ErrorArea, figL2ErrorArea, simPlotFolderName + "l2ErrorArea" + ".png" )

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

def plotMaxAndL2Error( axMaxErrorHandle, axL2ErrorHandle, xvals, maxError, l2Error, regexVal ):

    labelName = regexVal + " $L^{\infty}$ Error_" + "2h"
    axMaxErrorHandle.loglog( xvals[:-1], maxError[:-1], "-o", label = labelName )

    labelName = regexVal + " $L^{\infty}$ Error_" + "h"
    axMaxErrorHandle.loglog( xvals[:-1], maxError[1:], "-o", label = labelName )

    labelName = regexVal + " $L^{2}$ Error_" + "2h"
    axL2ErrorHandle.loglog( xvals[:-1], l2Error[:-1], "-o", label = labelName )

    labelName = regexVal + " $L^{2}$ Error_" + "h"
    axL2ErrorHandle.loglog( xvals[:-1], l2Error[1:], "-o", label = labelName )

def showDealiiPlot( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfiletriangleunstruct.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
        allSortedMeshVals.append( sortedMeshVals )

    minMaxRangeVals = getDealiiMinMaxRange( allSortedMeshVals )

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    allMaxErrorList = []
    allL2ErrorList = []

    figMaxErrorNumNodes = plt.figure()
    axMaxErrorNumNodes = figMaxErrorNumNodes.add_subplot(1, 1, 1)
    figL2ErrorNumNodes = plt.figure()
    axL2ErrorNumNodes = figL2ErrorNumNodes.add_subplot(1, 1, 1)

    figMaxErrorArea = plt.figure()
    axMaxErrorArea = figMaxErrorArea.add_subplot(1, 1, 1)
    figL2ErrorArea = plt.figure()
    axL2ErrorArea = figL2ErrorArea.add_subplot(1, 1, 1)

    # print(allSortedMeshVals)

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        curMaxErrorList = []
        curL2ErrorList = []
        numNodeVals = []
        areaVals = []
        regexVal = regexVals[idx]

        for index, meshval in enumerate( sortedMeshVals ):

            numNodeVals.append( getNumNodes( meshpath + meshval ) )
            areaVals.append( getAverage2DArea( meshpath + meshval, regexVal ) )

            solvalsFileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, "solutionvaluesGaussModified_", ".h5" )
            
            (nodes, solution) = getDealiiData( solvalsFileName )
            xvals = nodes[:, 0]
            yvals = nodes[:, 1]
            errvals = getDealiiError( nodes, solution )
            # print(solvalsFileName)

            curmaxval = np.max(errvals)
            curminval = np.min(errvals)

            curMaxErrorList.append( curmaxval )
            curL2ErrorList.append( np.sum( np.array(errvals)**2 ) ) 

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
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "errorContour_", ".png" )
            plt.savefig( plotfilename )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            levels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, errvals, levels = levels, colors = 'r')
            plt.colorbar()
            plotfilename = getFileNameFromMeshName( meshval, curPlotFolderName, "largeErrorContour_", ".png" )
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

        plotMaxAndL2Error( axMaxErrorNumNodes, axL2ErrorNumNodes, numNodeVals, curMaxErrorList, curL2ErrorList, regexVal )
        plotMaxAndL2Error( axMaxErrorArea, axL2ErrorArea, areaVals, curMaxErrorList, curL2ErrorList, regexVal ) 

    h2vals = dxvals**2

    print( areaVals )
    makePlotAdjustmentsAndSave( axMaxErrorNumNodes, figMaxErrorNumNodes, simPlotFolderName + "maxErrorNumNodes" + ".png" )
    makePlotAdjustmentsAndSave( axL2ErrorNumNodes, figL2ErrorNumNodes, simPlotFolderName + "l2ErrorNumNodes" + ".png" )

    makePlotAdjustmentsAndSave( axMaxErrorArea, figMaxErrorArea, simPlotFolderName + "maxErrorArea" + ".png" )
    makePlotAdjustmentsAndSave( axL2ErrorArea, figL2ErrorArea, simPlotFolderName + "l2ErrorArea" + ".png" )

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

def plotMaxAndL2ErrorComparison( axMaxErrorHandle, axL2ErrorHandle, xvals, maxError1,
                                 maxError2, l2Error1, l2Error2, h2vals, label1, label2 ):


    axMaxErrorHandle.loglog( xvals, maxError1, "-o", label = label1 + "$L^{\infty}$ Error_" + "h" )
    axMaxErrorHandle.loglog( xvals, maxError2, "-o", label = label2 + "$L^{\infty}$ Error_" + "h" )

    axL2ErrorHandle.loglog( xvals, l2Error1, "-o", label = label1 + "$L^{2}$ Error_" + "h" )
    axL2ErrorHandle.loglog( xvals, l2Error2, "-o", label = label2 + "$L^{2}$ Error_" + "h" )

    axMaxErrorHandle.loglog( xvals, h2vals, "-x", label = "$h^2$" )
    axL2ErrorHandle.loglog( xvals, h2vals, "-x", label = "$h^2$" )

def getMaxL2ErrorFigHandles():

    figMaxError = plt.figure()
    axMaxError = figMaxError.add_subplot(1, 1, 1)
    figL2Error = plt.figure()
    axL2Error = figL2Error.add_subplot(1, 1, 1)

    return [ figMaxError, figL2Error, axMaxError, axL2Error ]


def compareDealiiFinch( simPlotFolderName, regexVals ):

    checkAndCreateFolder( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfiletriangleunstruct.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:
        sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
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

    h2vals = dxvals**2

    for idx, sortedMeshVals in enumerate( allSortedMeshVals ):

        dealiiCurMaxErrorList = []
        dealiiCurL2ErrorList = []

        regexVal = regexVals[idx]

        finchCurMaxErrorList = []
        finchCurL2ErrorList = []
        numNodeVals = []
        areaVals = []

        [ figMaxErrorNumNodes, figL2ErrorNumNodes, axMaxErrorNumNodes, axL2ErrorNumNodes ] = \
            getMaxL2ErrorFigHandles()

        [ figMaxErrorArea, figL2ErrorArea, axMaxErrorArea, axL2ErrorArea ] = \
            getMaxL2ErrorFigHandles()

        for index, meshval in enumerate( sortedMeshVals ):

            numNodeVals.append( getNumNodes( meshpath + meshval ) )
            areaVals.append( getAverage2DArea( meshpath + meshval, regexVal ) )

            dealiiSolvalsFileName = getFileNameFromMeshName( meshval, dealiiTextfoldername, "solutionvaluesGaussModified_", ".h5" )
            
            (nodes, solution) = getDealiiData( dealiiSolvalsFileName )
            dealiiErrvals = getDealiiError( nodes, solution )
            # print(solvalsFileName)

            dealiiCurmaxval = np.max(dealiiErrvals)

            dealiiCurMaxErrorList.append( dealiiCurmaxval )
            dealiiCurL2ErrorList.append( np.sum( np.array( dealiiErrvals )**2 ) ) 

            # Finch data
            finchXvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "xvalues_", ".txt" )
            finchYvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "yvalues_", ".txt" )
            finchErrvalsFileName = getFileNameFromMeshName( meshval, finchTextfoldername, "errorvalues_", ".txt" )

            finchXvals = getData( finchXvalsFileName )
            finchYvals = getData( finchYvalsFileName )
            finchErrvals = getData( finchErrvalsFileName )

            # print(errvalsFileName)

            finchCurmaxval = np.max( finchErrvals )
            finchCurMaxErrorList.append( finchCurmaxval )
            finchCurL2ErrorList.append( np.sum( np.array( finchErrvals )**2 ) ) 

            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
            # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
            # print( np.max(errvals) )

        dealiiAllMaxErrorList.append( dealiiCurMaxErrorList )
        dealiiAllL2ErrorList.append( dealiiCurL2ErrorList )

        finchAllMaxErrorList.append( finchCurMaxErrorList )
        finchAllL2ErrorList.append( finchCurL2ErrorList )

        label1 = regexVals[idx] + " Dealii "
        label2 = regexVals[idx] + " Finch "

        plotMaxAndL2ErrorComparison( axMaxErrorNumNodes, axL2ErrorNumNodes, numNodeVals,
                                     dealiiCurMaxErrorList, finchCurMaxErrorList, dealiiCurL2ErrorList,
                                       finchCurL2ErrorList, h2vals, label1, label2 )
        
        plotMaxAndL2ErrorComparison( axMaxErrorArea, axL2ErrorArea, areaVals, 
                                    dealiiCurMaxErrorList, finchCurMaxErrorList, dealiiCurL2ErrorList, 
                                    finchCurL2ErrorList, h2vals, label1, label2 )


        curPlotFolderName = simPlotFolderName + "Finch_Dealii_Comparison" + "/"
        checkAndCreateFolder( curPlotFolderName )

        makePlotAdjustmentsAndSave( axMaxErrorNumNodes, figMaxErrorNumNodes,
                                    curPlotFolderName + regexVal + "_maxErrorNumNodes" + ".png" )
        makePlotAdjustmentsAndSave( axMaxErrorArea, figMaxErrorArea,
                                    curPlotFolderName + regexVal + "_l2Error" + ".png" )
        
    plt.show()

    return

# def compareParaview( simPlotFolderName, regexVals ):

#     assert( len(regexVals) == 2 )

#     if not os.path.exists(simPlotFolderName):
#         os.mkdir( simPlotFolderName )

#     meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
#     meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"
#     meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

#     allSortedMeshVals = []

#     for regexVal in regexVals:
#         sortedMeshVals = getSortedMeshVals( meshvals, regexVal )
#         allSortedMeshVals.append( sortedMeshVals )

#     minMaxRangeVals = getMinMaxRange( allSortedMeshVals )

#     for sortedMeshVals in allSortedMeshVals:
#         for index, meshval in enumerate( sortedMeshVals ):
#             mesh = meshio.read( meshpath + meshval )
#             meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

#     meshValsLen = len( allSortedMeshVals[0] )

#     allViews = dict()
#     for regexVal in regexVals:
#         allViews[regexVal] = CreateRenderView()
#         renderView = allViews[regexVal]
#         renderView.ViewSize = [701, 784]
#         renderView.InteractionMode = '2D'
#         renderView.AxesGrid = 'Grid Axes 3D Actor'
#         renderView.CenterOfRotation = [0.5, 0.5, 0.0]
#         renderView.StereoType = 'Crystal Eyes'
#         renderView.CameraPosition = [0.5, 0.5, 3.0349403797187358]
#         renderView.CameraFocalPoint = [0.5, 0.5, 0.0]
#         renderView.CameraFocalDisk = 1.0
#         renderView.CameraParallelScale = 0.7908298380174797
#         renderView.LegendGrid = 'Legend Grid Actor'

#     layout = CreateLayout(name='Layout #1')
#     layout.SplitHorizontal(0, 0.500000)

#     for idx, regexVal in enumerate(regexVals):
#         layout.AssignView( idx + 1, allViews[regexVal] )

#     layout.SetSize(1403, 784)
#     # layout.SetSize(1920, 1080)

#     for idx in range(meshValsLen):

#         curPlotFolderName = simPlotFolderName + "Plot" + str(idx) + "/"   
#         if not os.path.exists( curPlotFolderName ):
#             os.mkdir( curPlotFolderName )

#         solfiles = []
#         gmshfiles = []

#         for regexIdx, regexVal in enumerate( regexVals ):
            
#             SetActiveView( allViews[regexVal] )
#             myview = GetActiveView()

#             meshval = allSortedMeshVals[regexIdx][idx]
#             meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
#             filename = getFileNameFromMeshName( meshval, textfoldername, "errorAndu", ".vtu" )

#             solfile = OpenDataFile( filename )
#             solfiles.append(solfile)
#             display = Show(solfile)
#             ColorBy(display, ('POINTS', 'err'))
#             minVal, maxVal = minMaxRangeVals[ idx ]
#             colorMap = GetColorTransferFunction('err')
#             colorMap.RescaleTransferFunction( minVal, maxVal )        

#             gmshfile = OpenDataFile( meshvtkName )
#             gmshfiles.append(gmshfile)
#             dpGmsh = GetDisplayProperties( gmshfile )
#             dpGmsh.Representation = 'Wireframe'
#             gmshdisplay = Show(gmshfile)
            
#             dpSol = GetDisplayProperties(solfile, myview)
#             # # to show the color legend
#             dpSol.SetScalarBarVisibility(myview, True)
#             myview.Update()

#         Render()
#         plotfilename = curPlotFolderName + "paraview_error_comparison.png" 
#         SaveScreenshot( plotfilename, layout)        

#         for regexIdx, regexVal in enumerate( regexVals ):

#             SetActiveView( allViews[regexVal] )
#             Hide( solfiles[regexIdx] )
#             Hide( gmshfiles[regexIdx] )

#     return

if __name__ == "__main__":

    # gmshFileCmdNames = ["triangleMeshv1", "triangleMeshv2"]
    # regexVals = ["triangleMeshStruct", "triangleMeshUnstruct"]
    regexVals = ["triangleMeshUnstruct", "regular"]
    gmshFileCmdNames = ["triangleMeshv2", "regularMeshv3"]
    # gmshFileCmdNames = ["hangingMeshv8"]
    # regexVals = ["hanging"]

    simPlotRootFolderName = gmshImageFolderName + "PlotNumNodes_2pi/"
    meshPlotRootFolderName = gmshImageFolderName + "MeshPlots28_2pi/"

    # showMeshes( meshPlotRootFolderName, regexVals )

    runSim( simPlotRootFolderName, gmshFileCmdNames, regexVals )
    createMeshVTU( meshPlotRootFolderName, regexVals )
    compareDealiiFinch( simPlotRootFolderName, regexVals )

    simPlotFolderName = simPlotRootFolderName + "Dealii/"
    showDealiiPlot( simPlotFolderName, regexVals )

    simPlotFolderName = simPlotRootFolderName + "Finch/"
    showFinchPlot( simPlotFolderName, regexVals )    

    # showParaviewPlot( simPlotFolderName, regexVals )
    # showplot( simPlotFolderName, regexVals )
    # compareParaview( simPlotFolderName, regexVals )

