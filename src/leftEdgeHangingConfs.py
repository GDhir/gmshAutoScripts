# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
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
renderView1.CameraPosition = [0.2347476171737969, 1.4818594258074875, 3.687749310103378]
renderView1.CameraFocalPoint = [0.49999999999999956, 0.4999999999999992, 0.49999999999999956]
renderView1.CameraViewUp = [-0.1411465107891691, 0.9427609737891935, -0.30212482319745654]
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
nodeConf1SingleEdgevtu = XMLUnstructuredGridReader(registrationName='nodeConf1SingleEdge.vtu', FileName=['/home/gaurav/gmshAutoScripts/build/nodeConf1SingleEdge.vtu'])
nodeConf1SingleEdgevtu.CellArrayStatus = ['gmsh:physical', 'gmsh:geometrical']
nodeConf1SingleEdgevtu.TimeArray = 'None'

# create a new 'XML Unstructured Grid Reader'
nodeConf2SingleEdgevtu = XMLUnstructuredGridReader(registrationName='nodeConf2SingleEdge.vtu', FileName=['/home/gaurav/gmshAutoScripts/build/nodeConf2SingleEdge.vtu'])
nodeConf2SingleEdgevtu.CellArrayStatus = ['gmsh:physical', 'gmsh:geometrical']
nodeConf2SingleEdgevtu.TimeArray = 'None'

# create a new 'Extract Cells By Type'
extractCellsByType1 = ExtractCellsByType(registrationName='ExtractCellsByType1', Input=nodeConf1SingleEdgevtu)
extractCellsByType1.CellTypes = ['Pyramid']

# create a new 'Extract Cells By Type'
extractCellsByType2 = ExtractCellsByType(registrationName='ExtractCellsByType2', Input=nodeConf2SingleEdgevtu)
extractCellsByType2.CellTypes = ['Pyramid']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from nodeConf2SingleEdgevtu
nodeConf2SingleEdgevtuDisplay = Show(nodeConf2SingleEdgevtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
nodeConf2SingleEdgevtuDisplay.Representation = 'Point Gaussian'
nodeConf2SingleEdgevtuDisplay.ColorArrayName = [None, '']
nodeConf2SingleEdgevtuDisplay.SelectTCoordArray = 'None'
nodeConf2SingleEdgevtuDisplay.SelectNormalArray = 'None'
nodeConf2SingleEdgevtuDisplay.SelectTangentArray = 'None'
nodeConf2SingleEdgevtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
nodeConf2SingleEdgevtuDisplay.SelectOrientationVectors = 'None'
nodeConf2SingleEdgevtuDisplay.ScaleFactor = 0.1
nodeConf2SingleEdgevtuDisplay.SelectScaleArray = 'None'
nodeConf2SingleEdgevtuDisplay.GlyphType = 'Arrow'
nodeConf2SingleEdgevtuDisplay.GlyphTableIndexArray = 'None'
nodeConf2SingleEdgevtuDisplay.GaussianRadius = 0.02
nodeConf2SingleEdgevtuDisplay.SetScaleArray = [None, '']
nodeConf2SingleEdgevtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
nodeConf2SingleEdgevtuDisplay.OpacityArray = [None, '']
nodeConf2SingleEdgevtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
nodeConf2SingleEdgevtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
nodeConf2SingleEdgevtuDisplay.PolarAxes = 'PolarAxesRepresentation'
nodeConf2SingleEdgevtuDisplay.ScalarOpacityUnitDistance = 0.45544897693044056
nodeConf2SingleEdgevtuDisplay.OpacityArrayName = ['CELLS', 'gmsh:geometrical']
nodeConf2SingleEdgevtuDisplay.SelectInputVectors = [None, '']
nodeConf2SingleEdgevtuDisplay.WriteLog = ''

# show data from extractCellsByType2
extractCellsByType2Display = Show(extractCellsByType2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
extractCellsByType2Display.Representation = 'Wireframe'
extractCellsByType2Display.ColorArrayName = ['POINTS', '']
extractCellsByType2Display.LineWidth = 2.0
extractCellsByType2Display.SelectTCoordArray = 'None'
extractCellsByType2Display.SelectNormalArray = 'None'
extractCellsByType2Display.SelectTangentArray = 'None'
extractCellsByType2Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractCellsByType2Display.SelectOrientationVectors = 'None'
extractCellsByType2Display.ScaleFactor = 0.1
extractCellsByType2Display.SelectScaleArray = 'None'
extractCellsByType2Display.GlyphType = 'Arrow'
extractCellsByType2Display.GlyphTableIndexArray = 'None'
extractCellsByType2Display.GaussianRadius = 0.005
extractCellsByType2Display.SetScaleArray = [None, '']
extractCellsByType2Display.ScaleTransferFunction = 'PiecewiseFunction'
extractCellsByType2Display.OpacityArray = [None, '']
extractCellsByType2Display.OpacityTransferFunction = 'PiecewiseFunction'
extractCellsByType2Display.DataAxesGrid = 'GridAxesRepresentation'
extractCellsByType2Display.PolarAxes = 'PolarAxesRepresentation'
extractCellsByType2Display.ScalarOpacityUnitDistance = 1.0911236359717214
extractCellsByType2Display.OpacityArrayName = ['CELLS', 'gmsh:geometrical']
extractCellsByType2Display.SelectInputVectors = [None, '']
extractCellsByType2Display.WriteLog = ''

# ----------------------------------------------------------------
# restore active source
SetActiveSource(extractCellsByType2)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')