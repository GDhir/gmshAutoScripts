from os import listdir
from os.path import isfile, join
import subprocess
from math import sin, cos, pi
import os
from paraview.simple import *
import meshio
import folderUtils
import meshFileUtils
import fileNameUtils
import regexUtils
import dataUtils
import vtk
import errorUtils
paraview.simple._DisableFirstRenderCameraReset()


def showMeshes( simPlotFolderName, regexVals ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )
    
    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    allSortedMeshVals = []

    for regexVal in regexVals:

        regexCriterias = regexUtils.getRegexCriterias( regexVal )
        sortedMeshVals = meshFileUtils.getSortedMeshVals( meshvals, regexVal, regexCriterias )
        allSortedMeshVals.append( sortedMeshVals )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    for typeIdx, sortedMeshVals in enumerate( allSortedMeshVals ):
        
        regexVal = regexVals[typeIdx]
        regexCriterias = regexUtils.getRegexCriterias( regexVal )

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
            plotfilename = fileNameUtils.getFileNameFromMeshName( meshval, curPlotFolderName, 
                                                   "paraview_error", regexVal, 
                                                   regexCriterias, ".png" )
            SaveScreenshot( plotfilename, myview)

            Hide( gmshfile )

def showParaviewPlotFinch( simPlotFolderName, allParams, optionsParam, comparisonParam, meshvals, meshPath ):

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    levelsArr = meshFileUtils.getAllLevels( meshvals )
    levelsArr = levelsArr[:-2]

    juliaVarName = "errorvalues"
    minMaxRangeVals = dataUtils.getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = paramValue
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str(level)]
            criteriaValsStr = regexUtils.getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, "" )

            mesh = meshio.read( meshPath + meshFileName )
            meshio.write( meshVizPath + meshFileName[:-4] + ".vtu", mesh )

            meshvtkName = meshVizPath + meshFileName[:-4] + ".vtu"

            solValsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "errorAndu", "vtu" )

            solfile = OpenDataFile( solValsFileName )
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

            plotVarName = "paraview_error"
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )

            SaveScreenshot( plotfilename, myview)

            Hide( solfile )
            Hide( gmshfile )

def showParaviewPlotDealii( simPlotFolderName, allParams, optionsParam,
                        comparisonParam, meshvals, meshPath, negative = -1, pival = 2*pi ):

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    levelsArr = meshFileUtils.getAllLevels( meshvals )

    minMaxRangeVals = dataUtils.getDealiiMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, negative, pival )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    for idx, paramValue in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = paramValue
        regexCriterias = regexUtils.getRegexCriterias( "lvl" )

        for index, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )
            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str(level)]

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, "" )

            mesh = meshio.read( meshPath + meshFileName )
            meshio.write( meshVizPath + meshFileName[:-4] + ".vtu", mesh )

            meshvtkName = meshVizPath + meshFileName[:-4] + ".vtu"

            dealiiVarName = "solutionvalues"
            solvalsFileName = fileNameUtils.getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, dealiiVarName, "vtu" ) 
            # print(nodes.shape)

            dealiiVarName = "errorvalues"
            errorvalsFileName = fileNameUtils.getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, dealiiVarName, "vtu" )
            createDataVTK( solvalsFileName, errorvalsFileName, negative, pival )

            solfile = OpenDataFile( errorvalsFileName )
            display = Show(solfile)
            ColorBy(display, ('POINTS', 'solution'))
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

            plotVarName = "paraview_error"
            plotfilename = fileNameUtils.getTextFileName( curPlotFolderName, pythonVarName, plotVarName, "png" )

            SaveScreenshot( plotfilename, myview)

            Hide( solfile )
            Hide( gmshfile )

def compareParaview( simPlotFolderName, allParams, optionsParam, comparisonParam,
                    meshvals, software = "Finch" ):

    assert( len(allParams[ comparisonParam ]) == 2 )

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"

    levelsArr = meshFileUtils.getAllLevels( meshvals )
    minMaxRangeVals = dataUtils.getFinchMinMaxRange( levelsArr, allParams, optionsParam,
                                                     comparisonParam, "errorvalues" )

    for sortedMeshVals in allSortedMeshVals:
        for index, meshval in enumerate( sortedMeshVals ):
            mesh = meshio.read( meshpath + meshval )
            meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    meshValsLen = len( levelsArr )

    allViews = dict()
    for paramVal in allParams[ comparisonParam ]:
        allViews[ paramVal ] = CreateRenderView()
        renderView = allViews[ paramVal ]
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

    for idx, paramVal in enumerate( allParams[ comparisonParam ] ):
        layout.AssignView( idx + 1, allViews[ paramVal ] )

    layout.SetSize(1403, 784)
    # layout.SetSize(1920, 1080)

    for idx in range( meshValsLen ):

        curPlotFolderName = simPlotFolderName + "Plot" + str(idx) + "/"   
        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        solfiles = []
        gmshfiles = []
        level = levelsArr[idx]
        optionsParam[ "level" ] = str( level )

        for paramIdx, paramVal in enumerate( allParams[ comparisonParam ] ):
            
            SetActiveView( allViews[ paramVal ] )
            myview = GetActiveView()
            optionsParam[ comparisonParam ] = paramVal

            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            regexCriteriaVals = [str(level)]
            criteriaValsStr = regexUtils.getCriteriaValsString( regexCriterias, regexCriteriaVals )

            meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, "" )

            mesh = meshio.read( meshPath + meshFileName )
            meshio.write( meshVizPath + meshFileName[:-4] + ".vtu", mesh )

            meshvtkName = meshVizPath + meshFileName[:-4] + ".vtu"
            solValsFileName = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, "errorAndu", "vtu" )

            solfile = OpenDataFile( solValsFileName )
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

        for paramIdx, paramVal in enumerate( allParams[comparisonParam] ):

            SetActiveView( allViews[paramVal] )
            Hide( solfiles[paramIdx] )
            Hide( gmshfiles[paramIdx] )

    return

def createMeshVTU( simPlotFolderName, regexVals ):

    folderUtils.checkAndCreateFolder( simPlotFolderName )
    
    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"
    
    for index, meshval in enumerate( meshvals ):
        mesh = meshio.read( meshpath + meshval )
        meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

def createDataVTK( solvalsFileName, errvalsFileName, negative = -1, pival = 2*pi ):

    mesh = meshio.read( solvalsFileName )
    nodes = mesh.points
    # cells = mesh.cells[0]
    # print(cells)
    solution = mesh.point_data['solution']
    errvals = errorUtils.getDealiiError( nodes, solution, negative, pival )

    mesh.point_data['solution'] = errvals
    meshio.write( errvalsFileName, mesh )

    return

def buildAllMeshes( gmshFileCmdNames, gmshfileArgs ):

    for gmshFileCmd in gmshFileCmdNames:
        buildMesh( gmshFileCmd, gmshfileArgs )

def buildMesh( gmshfilecmd, gmshfileargs ):

    gmshbuildFolder = "/home/gaurav/gmshAutoScripts/build/"
    compilecmd = ["make", "-j", "4", gmshfilecmd]
    subprocess.run( compilecmd, cwd = gmshbuildFolder )

    gmshfileargs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshlstcmd = [ gmshbuildFolder + gmshfilecmd, gmshfileargs ]
    subprocess.run( meshlstcmd, cwd = gmshbuildFolder )



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
    # regexVals = ["mesh"]

    allParams = dict()
    allParams["software"] = ["Finch", "Dealii"]
    allParams["meshRegexVal"] = regexVals
    allParams["gmshFileCmdName"] = gmshFileCmdNames
    allParams["quadratureOrder"] = [ "2", "3", "4" ]
    allParams["sin(kpix)"] = [ "1", "2", "4" ]
    allParams["coeff_F"] = [ "-2", "-8", "-32" ]

    optionsParam = dict()
    optionsParam["quadratureOrder"] = "2"
    optionsParam["sin(kpix)"] = "2"
    optionsParam["coeff_F"] = "-8"
    optionsParam["software"] = "Finch"
    optionsParam["meshRegexVal"] = regexVals[0]
    optionsParam["level"] = "0"

    # pythonVarName = getPythonVarName( optionsParam )
    # print(pythonVarName)
    # filename = getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, "solutionValues", "vtu" )
    # print(filename) 

    comparisonParam = "meshRegexVal"

    simPlotRootFolderName = folderUtils.gmshImageFolderName + "PlotGaussQuad_2pi/"
    folderUtils.checkAndCreateFolder(simPlotRootFolderName)
    meshPlotRootFolderName = folderUtils.gmshImageFolderName + "MeshPlotsHangingLevel_QuadratureOrder=4_pi/"

    # regexVals = [ "mesh" ]
    meshPath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    # showMeshes( folderUtils.meshPlotRootFolderName, regexVals )

    # createMeshVTU( meshPlotRootFolderName, regexVals )

    # simPlotFolderName = simPlotRootFolderName + "Dealii/"
    # print( "Dealii" )
    # showDealiiPlot( simPlotFolderName, regexVals, meshpath= meshPath, varName = "solutionvalues", negative = -1, pival = 2*pi )

    simPlotFolderName = simPlotRootFolderName + "Finch/"
    print( "Finch" )
    buildAllMeshes( gmshFileCmdNames, meshPath )
    meshArr = meshFileUtils.getMeshFilesFromFolder( meshPath )

    # showParaviewPlotFinch( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath )
    # compareParaview( simPlotFolderName, regexVals, meshPath )

    # setFinchTriangleQuadrature( 2 )

    level = 3
    pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
    regexCriterias = regexUtils.getRegexCriterias( "lvl" )
    regexCriteriaVals = [str( level )]

    meshFileName = fileNameUtils.getMeshFileName( optionsParam, regexCriterias, regexCriteriaVals, meshPath )

    simPlotFolderName = simPlotRootFolderName + "Dealii/"
    showParaviewPlotDealii( simPlotFolderName, allParams, optionsParam, comparisonParam, meshArr, meshPath, -1, 2*pi )
    