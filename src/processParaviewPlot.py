from os import listdir
from os.path import isfile, join
from math import sin, cos, pi
import os
from paraview.simple import *
import meshio
import folderUtils
paraview.simple._DisableFirstRenderCameraReset()


def showMeshes( simPlotFolderName, regexVals ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )
    
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
            folderUtils.checkAndCreateFolder( curPlotFolderName )
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
            hangingfilename = getFileNameFromMeshName( meshval, folderUtils.textFolderNames[ software ],
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
            filename = getFileNameFromMeshName( meshval, folderUtils.textFolderNames[software],
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

def createMeshVTU( simPlotFolderName, regexVals ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )
    
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

