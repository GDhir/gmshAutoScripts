import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
sys.path.append( "/home/gaurav/.local/lib/python3.9/site-packages/" )
import gmsh
import math
import os
import sys
import gmshUtils
import meshio
from paraview.simple import *
import vtk
import paraview

# state file generated using paraview version 5.11.1
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

def createGMSHPNG( folderName, vtuFileName, plotFileName ):

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # ----------------------------------------------------------------
    # setup views used in the visualization
    # ----------------------------------------------------------------

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # Create a new 'Render View'
    renderView1 = CreateView('RenderView')
    renderView1.ViewSize = [1016, 789]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraPosition = [0.6140295491416724, 1.9877760738188548, 3.4949410743454217]
    renderView1.CameraFocalPoint = [0.4999999999999997, 0.49999999999999983, 0.49999999999999983]
    renderView1.CameraViewUp = [0.014716823347220036, 0.895263407023515, -0.4452941130929253]
    renderView1.CameraFocalDisk = 1.0
    renderView1.CameraParallelScale = 0.8660254037844386
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.OSPRayMaterialLibrary = materialLibrary1

    SetActiveView(None)

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.AssignView(0, renderView1)
    layout1.SetSize(1016, 789)

    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(renderView1)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'XML Unstructured Grid Reader'
    nodeConf3vtu = XMLUnstructuredGridReader(registrationName=vtuFileName, FileName=[folderName + vtuFileName])
    nodeConf3vtu.CellArrayStatus = ['gmsh:physical', 'gmsh:geometrical']
    nodeConf3vtu.TimeArray = 'None'

    # create a new 'Extract Cells By Type'
    extractCellsByType1 = ExtractCellsByType(registrationName='ExtractCellsByType1', Input=nodeConf3vtu)
    extractCellsByType1.CellTypes = ['Pyramid']

    # create a new 'Extract Edges'
    extractEdges1 = ExtractEdges(registrationName='ExtractEdges1', Input=nodeConf3vtu)

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from nodeConf3vtu
    nodeConf3vtuDisplay = Show(nodeConf3vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    nodeConf3vtuDisplay.Representation = 'Point Gaussian'
    nodeConf3vtuDisplay.ColorArrayName = [None, '']
    nodeConf3vtuDisplay.SelectTCoordArray = 'None'
    nodeConf3vtuDisplay.SelectNormalArray = 'None'
    nodeConf3vtuDisplay.SelectTangentArray = 'None'
    nodeConf3vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    nodeConf3vtuDisplay.SelectOrientationVectors = 'None'
    nodeConf3vtuDisplay.ScaleFactor = 0.1
    nodeConf3vtuDisplay.SelectScaleArray = 'None'
    nodeConf3vtuDisplay.GlyphType = 'Arrow'
    nodeConf3vtuDisplay.GlyphTableIndexArray = 'None'
    nodeConf3vtuDisplay.GaussianRadius = 0.02
    nodeConf3vtuDisplay.SetScaleArray = [None, '']
    nodeConf3vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    nodeConf3vtuDisplay.OpacityArray = [None, '']
    nodeConf3vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    nodeConf3vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
    nodeConf3vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
    nodeConf3vtuDisplay.ScalarOpacityUnitDistance = 0.36950037695844756
    nodeConf3vtuDisplay.OpacityArrayName = ['CELLS', 'gmsh:geometrical']
    nodeConf3vtuDisplay.SelectInputVectors = [None, '']
    nodeConf3vtuDisplay.WriteLog = ''

    # show data from extractEdges1
    extractEdges1Display = Show(extractEdges1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractEdges1Display.Representation = 'Surface'
    extractEdges1Display.AmbientColor = [0.0, 0.6666666666666666, 0.0]
    extractEdges1Display.ColorArrayName = [None, '']
    extractEdges1Display.DiffuseColor = [0.0, 0.6666666666666666, 0.0]
    extractEdges1Display.SelectTCoordArray = 'None'
    extractEdges1Display.SelectNormalArray = 'None'
    extractEdges1Display.SelectTangentArray = 'None'
    extractEdges1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractEdges1Display.SelectOrientationVectors = 'None'
    extractEdges1Display.ScaleFactor = 0.1
    extractEdges1Display.SelectScaleArray = 'None'
    extractEdges1Display.GlyphType = 'Arrow'
    extractEdges1Display.GlyphTableIndexArray = 'None'
    extractEdges1Display.GaussianRadius = 0.005
    extractEdges1Display.SetScaleArray = [None, '']
    extractEdges1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractEdges1Display.OpacityArray = [None, '']
    extractEdges1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractEdges1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractEdges1Display.PolarAxes = 'PolarAxesRepresentation'
    extractEdges1Display.SelectInputVectors = [None, '']
    extractEdges1Display.WriteLog = ''

    SaveScreenshot( folderName + plotFileName + "AllEdges.png", renderView1)

    # show data from extractCellsByType1
    extractCellsByType1Display = Show(extractCellsByType1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractCellsByType1Display.Representation = 'Surface'
    extractCellsByType1Display.ColorArrayName = [None, '']
    extractCellsByType1Display.SelectTCoordArray = 'None'
    extractCellsByType1Display.SelectNormalArray = 'None'
    extractCellsByType1Display.SelectTangentArray = 'None'
    extractCellsByType1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractCellsByType1Display.SelectOrientationVectors = 'None'
    extractCellsByType1Display.ScaleFactor = 0.1
    extractCellsByType1Display.SelectScaleArray = 'None'
    extractCellsByType1Display.GlyphType = 'Arrow'
    extractCellsByType1Display.GlyphTableIndexArray = 'None'
    extractCellsByType1Display.GaussianRadius = 0.005
    extractCellsByType1Display.SetScaleArray = [None, '']
    extractCellsByType1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractCellsByType1Display.OpacityArray = [None, '']
    extractCellsByType1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractCellsByType1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractCellsByType1Display.PolarAxes = 'PolarAxesRepresentation'
    extractCellsByType1Display.ScalarOpacityUnitDistance = 1.0129094569634634
    extractCellsByType1Display.OpacityArrayName = ['CELLS', 'gmsh:geometrical']
    extractCellsByType1Display.SelectInputVectors = [None, '']
    extractCellsByType1Display.WriteLog = ''

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(nodeConf3vtu)
    # ----------------------------------------------------------------
    SaveScreenshot( folderName + plotFileName + "AllPyramids.png", renderView1)


if __name__ == "__main__":

    folderName = "/home/gaurav/gmshAutoScripts/build/"
    plotFileName = "nodeConf1"

    createGMSHPNG( folderName, "nodeConf1.vtu", plotFileName )
    # SaveExtracts(ExtractsOutputDirectory='extracts')

