import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh
import math
import os
import sys
import gmshUtils
import meshio
import numpy as np
import miscUtils
import pandas as pd
import meshDataUtils
import rotateReflect
import gmshConfGenerate
import meshFileUtils
from matplotlib import pyplot as plt

def runConfs( folderName ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    linesVal = []
    ptsVal = []

    lc = 0.5

    allpts = [ (0, 0, 0, lc), (0, 0, 0.5, lc), (0, 0, 1, lc), (0, 0.5, 1, lc), (0, 1, 1, lc), (0, 1, 0.5, lc), 
                (0, 1, 0, lc),(0, 0.5, 0, lc), (0, 0.5, 0.5, lc), (1, 0, 0, 2 * lc), (1, 0, 1, 2 * lc), 
                (1, 1, 1, 2 * lc), (1, 1, 0, 2 * lc) ]

    lineIdxVals = [ (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0), (1, 8), (8, 5),\
                (7, 8), (8, 3), (0, 9), (2, 10), (6, 12), (4, 11), (9, 10), (10, 11), (11, 12), (12, 9) ]

    for ptval in allpts:
        ptsVal.append( gmsh.model.occ.addPoint( ptval[0], ptval[1], ptval[2], ptval[3] ) )

    for lineIdxVal in lineIdxVals:
        linesVal.append( gmsh.model.occ.addLine( ptsVal[ lineIdxVal[0] ], ptsVal[ lineIdxVal[1] ] ) )

    surfaceCurveLoopVals = [ (0, 8, 10, 7), (1, 2, 11, 8), (11, 3, 4, 9), (10, 9, 5, 6), (16, 17, 18, 19), (12, 7, 6, 14, 19),\
                            (13, 2, 3, 15, 17), (12, 16, 13, 1, 0), (14, 18, 15, 4, 5) ]

    orientationVals = [ (1, 1, -1, 1), (1, 1, -1, -1), (1, 1, 1, -1), (1, 1, 1, 1), (1, 1, 1, 1), (-1, -1, -1, 1, 1),\
                    (-1, 1, 1, 1, -1), (1, 1, -1, -1, -1), (1, -1, -1, 1, 1) ]

    surfaceVals = []

    for idx, clVals in enumerate( surfaceCurveLoopVals ):

        oval = orientationVals[idx]

        lineIdxVals = []
        for lineIdx, lineIdxVal in enumerate( clVals ):

            lineIdxVals.append( oval[ lineIdx] * linesVal[ lineIdxVal ] )
        
        cl = gmsh.model.occ.addCurveLoop( lineIdxVals )
        surfaceVals.append( gmsh.model.occ.addPlaneSurface( [ cl ] ) )

    gmsh.model.occ.synchronize()
    gmshUtils.setTransfiniteCurves( [linesVal], 2, occ = True )

    cornerPtsAllSurfaces = [ (0, 1, 8, 7), (1, 2, 3, 8), (8, 3, 4, 5), (7, 8, 5, 6), (9, 10, 11, 12), (0, 6, 12, 9),\
                            (2, 4, 11, 10), (0, 9, 10, 2), (4, 6, 12, 11) ]

    for sIdx, surfaceVal in enumerate( surfaceVals[:5] ):

        cornerPtVals = []

        for cornerPtIdx in cornerPtsAllSurfaces[sIdx]:
            cornerPtVals.append( ptsVal[ cornerPtIdx ] )

        # print(cornerPtVals)
        gmsh.model.occ.synchronize()
        gmshUtils.setTransfiniteSurfaces( [surfaceVal], cornerPtVals, occ = True )

    sl = gmsh.model.occ.addSurfaceLoop( [ surfaceVals[0], surfaceVals[1], surfaceVals[2], surfaceVals[3],\
                                        surfaceVals[4], surfaceVals[5], surfaceVals[6], surfaceVals[7], surfaceVals[8] ] )
                
    vl = gmsh.model.occ.addVolume( [sl] )
    # gmsh.model.geo.rotate([(3, vl)], 0, 0.3, 0, 0, 0, 1, -math.pi / 4)

    # gmsh.model.occ.synchronize()
    # gmsh.model.mesh.setTransfiniteVolume( vl )
    # gmshUtils.setTransfiniteSurfaces( surfaceVals )

    # gmsh.model.geo.addPoint(0.2, 0.2, 0.25, lc, 150)
    # gmsh.model.mesh.setTransfiniteAutomatic()
    gmsh.model.occ.synchronize()

    # gmsh.model.mesh.embed(0, [150], 2, 1)
    # gmsh.model.mesh.setCompound(1, [linesVal[6], linesVal[7], linesVal[12], linesVal[19], linesVal[14], \
                                # linesVal[13], linesVal[17], linesVal[15], linesVal[3], linesVal[2] ])
    # gmsh.model.mesh.setCompound( 2, [surfaceVals[5], surfaceVals[6] ] )
    # gmsh.model.mesh.setCompound( 2, [surfaceVals[7], surfaceVals[8] ] )
    # gmsh.model.geo.mesh.setSmoothing( 3, vl, 4 )

    # gmsh.model.mesh.setCompound( 2, [surfaceVals[0], surfaceVals[1], surfaceVals[2], surfaceVals[3] ] )

    gmshUtils.recombineSurfaces( surfaceVals[:5] )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
    gmsh.option.setNumber("Mesh.Algorithm", 9)
    # gmsh.model.mesh.refine()
    # gmsh.view.option.setNumber("DrawPoints", 1)
    t2 = gmsh.view.add("Hanging Node Conf")

    gmsh.option.setNumber("Geometry.Points", 1)
    gmsh.option.setNumber("Geometry.PointSize", 20)

    gmsh.option.setNumber("Mesh.Points", 1)
    gmsh.option.setNumber("Mesh.PointSize", 10)

    gmsh.view.option.setNumber(t2, "PointType", 2)
    gmsh.view.option.setNumber(t2, "PointSize", 10)

    # gmsh.option.setNumber("General.RotationX",
    #                           30)
    # gmsh.option.setNumber("General.RotationY",
    #                         gmsh.option.getNumber("General.RotationX") / 2)
    # gmsh.option.setNumber("General.RotationZ",
    #                         gmsh.option.getNumber("General.RotationZ") + 0.1)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = folderName + "nodeConf2" + ".msh"
    gmsh.write( regMeshFileName )

    # gmsh.graphics.draw()
    # gmsh.model.setColor( [(0, ptVal) for ptVal in ptsVal], 2, 2, 127)  # Gray50

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.fltk.initialize()

    pngFileName = folderName + "nodeConf2.png"
    gmsh.write( pngFileName )

    vtuFileName = "nodeConf2"
    meshFileUtils.createGMSHVTU( folderName, vtuFileName )

    gmsh.finalize()


def runConfs2( folderName ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    linesVal = []
    ptsVal = []

    lc = 0.2

    allpts = [ (0, 0, 0, lc), (0, 0, 0.5, lc), (0, 0, 1, lc), (0, 0.5, 1, lc), (0, 1, 1, lc), (0, 1, 0.5, lc), 
                (0, 1, 0, lc),(0, 0.5, 0, lc), (0, 0.5, 0.5, lc), (1, 0, 0, 2 * lc), (1, 0, 1, 2 * lc), 
                (1, 1, 1, 2 * lc), (1, 1, 0, 2 * lc) ]

    lineIdxVals = [ (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0), (1, 8), (8, 5),\
                (7, 8), (8, 3), (0, 9), (2, 10), (6, 12), (4, 11), (9, 10), (10, 11), (11, 12), (12, 9) ]

    for ptval in allpts:
        ptsVal.append( gmsh.model.geo.addPoint( ptval[0], ptval[1], ptval[2], ptval[3] ) )

    for lineIdxVal in lineIdxVals:
        linesVal.append( gmsh.model.geo.addLine( ptsVal[ lineIdxVal[0] ], ptsVal[ lineIdxVal[1] ] ) )

    surfaceCurveLoopVals = [ (0, 1, 2, 3, 4, 5, 6, 7), (16, 17, 18, 19), (12, 7, 6, 14, 19),\
                            (13, 2, 3, 15, 17), (12, 16, 13, 1, 0), (14, 18, 15, 4, 5) ]

    orientationVals = [ (1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1), (-1, -1, -1, 1, 1),\
                    (-1, 1, 1, 1, -1), (1, 1, -1, -1, -1), (1, -1, -1, 1, 1) ]

    surfaceVals = []

    for idx, clVals in enumerate( surfaceCurveLoopVals ):

        oval = orientationVals[idx]

        lineIdxVals = []
        for lineIdx, lineIdxVal in enumerate( clVals ):

            lineIdxVals.append( oval[ lineIdx] * linesVal[ lineIdxVal ] )
        
        cl = gmsh.model.geo.addCurveLoop( lineIdxVals )
        surfaceVals.append( gmsh.model.geo.addPlaneSurface( [ cl ] ) )

    gmsh.model.occ.synchronize()
    gmshUtils.setTransfiniteCurves( [linesVal], 2 )

    cornerPtsAllSurfaces = [ (0, 2, 4, 6), (9, 10, 11, 12), (0, 6, 12, 9),\
                            (2, 4, 11, 10), (0, 9, 10, 2), (4, 6, 12, 11) ]

    t2 = gmsh.view.add("Hanging Node Conf")

    for sIdx, surfaceVal in enumerate( surfaceVals[:2] ):

        cornerPtVals = []

        for cornerPtIdx in cornerPtsAllSurfaces[sIdx]:
            cornerPtVals.append( ptsVal[ cornerPtIdx ] )

        print(cornerPtVals)
        gmshUtils.setTransfiniteSurfaces( [surfaceVal], cornerPtVals )

    sl = gmsh.model.geo.addSurfaceLoop( [ surfaceVals[0], surfaceVals[1], surfaceVals[2], surfaceVals[3],\
                                        surfaceVals[4], surfaceVals[5]] )
                
    vl = gmsh.model.geo.addVolume( [sl] )
    # gmsh.model.geo.mesh.setTransfiniteVolume( vl )
    # gmshUtils.setTransfiniteSurfaces( surfaceVals )

    gmsh.model.geo.synchronize()

    gmshUtils.recombineSurfaces( surfaceVals[:2] )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

    gmsh.option.setNumber("Geometry.Points", 1)
    gmsh.option.setNumber("Geometry.PointSize", 20)

    gmsh.option.setNumber("Mesh.Points", 1)
    gmsh.option.setNumber("Mesh.PointSize", 10)

    gmsh.view.option.setNumber(t2, "PointType", 2)
    gmsh.view.option.setNumber(t2, "PointSize", 10)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = folderName + "nodeConf3" + ".msh"
    gmsh.write( regMeshFileName )

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.fltk.initialize()

    pngFileName = folderName + "nodeConf2.png"
    gmsh.write( pngFileName )

    vtuFileName = "nodeConf3"
    meshFileUtils.createGMSHVTU( folderName, vtuFileName )

    gmsh.finalize()

def runConfsEdgeHanging( folderName, fileNamePrefix, algNumber = 3 ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    linesVal = []
    ptsVal = []

    lc = 0.5

    allpts = [ (0, 0, 0, lc), (0, 0, 0.5, lc), (0, 0, 1, lc), (0, 1, 1, lc), 
                (0, 1, 0, lc), (1, 0, 0, 2 * lc), (1, 0, 1, 2 * lc), 
                (1, 1, 1, 2 * lc), (1, 1, 0, 2 * lc) ]

    lineIdxVals = [ (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),\
                (0, 5), (2, 6), (4, 8), (3, 7), (5, 6), (6, 7), (7, 8), (8, 5) ]

    for ptval in allpts:
        ptsVal.append( gmsh.model.occ.addPoint( ptval[0], ptval[1], ptval[2], ptval[3] ) )

    for lineIdxVal in lineIdxVals:
        linesVal.append( gmsh.model.occ.addLine( ptsVal[ lineIdxVal[0] ], ptsVal[ lineIdxVal[1] ] ) )

    surfaceCurveLoopVals = [ (0, 1, 2, 3, 4), (9, 10, 11, 12), (5, 12, 7, 4), (6, 10, 8, 2), (5, 9, 6, 1, 0), (7, 11, 8, 3) ]

    orientationVals = [ (1, 1, 1, 1, 1), (1, 1, 1, 1), (1, 1, -1, 1), (1, 1, -1, 1), (1, 1, -1, -1, -1), (1, -1, -1, 1) ]

    surfaceVals = []

    for idx, clVals in enumerate( surfaceCurveLoopVals ):

        oval = orientationVals[idx]

        lineIdxVals = []
        for lineIdx, lineIdxVal in enumerate( clVals ):

            lineIdxVals.append( oval[ lineIdx] * linesVal[ lineIdxVal ] )
        
        cl = gmsh.model.occ.addCurveLoop( lineIdxVals )
        surfaceVals.append( gmsh.model.occ.addPlaneSurface( [ cl ] ) )

    gmsh.model.occ.synchronize()
    gmshUtils.setTransfiniteCurves( [linesVal], 2, occ = True )

    cornerPtsAllSurfaces = [ (0, 2, 3, 4), (5, 6, 7, 8), (0, 5, 8, 4), (2, 6, 7, 3), (0, 5, 6, 2), (4, 8, 7, 3) ]

    transfiniteSurfaceIdxVals = [ 1, 2, 3, 5 ]
    transfiniteSurfaceVals = []

    for sIdx in transfiniteSurfaceIdxVals:

        surfaceVal = surfaceVals[ sIdx ]
        cornerPtVals = []
        transfiniteSurfaceVals.append( surfaceVal )

        for cornerPtIdx in cornerPtsAllSurfaces[sIdx]:
            cornerPtVals.append( ptsVal[ cornerPtIdx ] )

        # print(cornerPtVals)
        gmsh.model.occ.synchronize()
        gmshUtils.setTransfiniteSurfaces( [surfaceVal], cornerPtVals, occ = True )

    sl = gmsh.model.occ.addSurfaceLoop( [ surfaceVals[0], surfaceVals[1], surfaceVals[2], surfaceVals[3], 
                                    surfaceVals[4], surfaceVals[5] ] )
                
    vl = gmsh.model.occ.addVolume( [sl] )

    gmsh.model.occ.synchronize()

    gmshUtils.recombineSurfaces( transfiniteSurfaceVals )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
    gmsh.option.setNumber("Mesh.Algorithm", algNumber)
    # gmsh.model.mesh.refine()
    # gmsh.view.option.setNumber("DrawPoints", 1)
    t2 = gmsh.view.add("Hanging Node Conf")

    gmsh.option.setNumber("Geometry.Points", 1)
    gmsh.option.setNumber("Geometry.PointSize", 20)

    gmsh.option.setNumber("Mesh.Points", 1)
    gmsh.option.setNumber("Mesh.PointSize", 10)

    gmsh.view.option.setNumber(t2, "PointType", 2)
    gmsh.view.option.setNumber(t2, "PointSize", 10)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = folderName + fileNamePrefix + ".msh"
    gmsh.write( regMeshFileName )

    # gmsh.graphics.draw()
    # gmsh.model.setColor( [(0, ptVal) for ptVal in ptsVal], 2, 2, 127)  # Gray50

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.fltk.initialize()

    pngFileName = folderName + fileNamePrefix + ".png"
    gmsh.write( pngFileName )

    vtuFileName = fileNamePrefix
    meshFileUtils.createGMSHVTU( folderName, vtuFileName )

    gmsh.finalize()


def runConfsTwoFacesHanging( folderName, fileNamePrefix, algNumber = 3 ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    linesVal = []
    ptsVal = []

    lc = 0.5

    allpts = [ (0, 0, 0, lc), (0.5, 0, 0, lc), (1, 0, 0, lc), (0, 0.5, 0, lc), (0, 1, 0, lc), (1, 1, 0, 2 * lc), 
                (0, 0, 0.5, lc),(0.5, 0, 0.5, lc), (1, 0, 0.5, lc), (0, 0.5, 0.5, lc), (0, 1, 0.5, lc), (0, 0, 1, lc), 
                (0.5, 0, 1, lc), (1, 0, 1, lc), (0, 0.5, 1, lc), (0, 1, 1, lc), (1, 1, 1, 2 * lc) ]

    lineIdxVals = [ (0, 1), (1, 2), (4, 5), (6, 7), (7, 8),\
                (11, 12), (12, 13), (15, 16), (0, 3), (3, 4),\
                (2, 5), (6, 9), (9, 10), (11, 14), (14, 15),\
                (13, 16), (0, 6), (6, 11), (1, 7), (7, 12),\
                (2, 8), (8, 13), (3, 9), (9, 14), (4, 10), (10, 15), (5, 16) ]

    for ptval in allpts:
        ptsVal.append( gmsh.model.occ.addPoint( ptval[0], ptval[1], ptval[2], ptval[3] ) )

    for lineIdxVal in lineIdxVals:
        linesVal.append( gmsh.model.occ.addLine( ptsVal[ lineIdxVal[0] ], ptsVal[ lineIdxVal[1] ] ) )

    surfaceCurveLoopVals = [ (16, 11, -22, -8), (22, 12, -24, -9), (17, 13, -23, -11), (23, 14, -25, -12),\
                        (20, 21, 15, -26, -10), (0, 18, -3, -16), (1, 20, -4, -18), (3, 19, -5, -17), (4, 21, -6, -19),\
                        (2, 26, -7, -25, -24), (0, 1, 10, -2, -9, -8), (5, 6, 15, -7, -14, -13)]

    # orientationVals = [ (1, 1, -1, -1), (1, 1, -1, -1), (1, 1, -1, -1), (1, 1, -1, -1), (1, 1, -1, -1, -1), (1, -1, -1, 1) ]

    surfaceVals = []

    for idx, clVals in enumerate( surfaceCurveLoopVals ):

        lineIdxVals = []
        for lineIdx, lineIdxVal in enumerate( clVals ):

            oval = 1

            if lineIdxVal < 0:
                oval = -1

            lineIdxVals.append( oval * linesVal[ abs( lineIdxVal ) ] )
        
        cl = gmsh.model.occ.addCurveLoop( lineIdxVals )
        surfaceVals.append( gmsh.model.occ.addPlaneSurface( [ cl ] ) )

    gmsh.model.occ.synchronize()
    gmshUtils.setTransfiniteCurves( [linesVal], 2, occ = True )

    cornerPtsAllSurfaces = [ (0, 6, 9, 3), (3, 9, 10, 4), (6, 11, 14, 9), (9, 14, 15, 10),\
                        (2, 13, 16, 5), (0, 1, 7, 6), (1, 2, 8, 7), (6, 7, 12, 11) , (7, 8, 13, 12),\
                        (4, 5, 16, 15), (0, 2, 5, 4), (11, 13, 16, 15)]

    transfiniteSurfaceIdxVals = [ 0, 1, 2, 3, 5, 6, 7, 8 ]
    transfiniteSurfaceVals = []

    for sIdx in transfiniteSurfaceIdxVals:

        surfaceVal = surfaceVals[ sIdx ]
        cornerPtVals = []
        transfiniteSurfaceVals.append( surfaceVal )

        for cornerPtIdx in cornerPtsAllSurfaces[sIdx]:
            cornerPtVals.append( ptsVal[ cornerPtIdx ] )

        # print(cornerPtVals)
        gmsh.model.occ.synchronize()
        gmshUtils.setTransfiniteSurfaces( [surfaceVal], cornerPtVals, occ = True )

    sl = gmsh.model.occ.addSurfaceLoop( [ surfaceVals[0], surfaceVals[1], surfaceVals[2], surfaceVals[3], 
                                    surfaceVals[4], surfaceVals[5] ] )
                
    vl = gmsh.model.occ.addVolume( [sl] )

    gmsh.model.occ.synchronize()

    gmshUtils.recombineSurfaces( transfiniteSurfaceVals )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
    gmsh.option.setNumber("Mesh.Algorithm", algNumber)
    # gmsh.model.mesh.refine()
    # gmsh.view.option.setNumber("DrawPoints", 1)
    t2 = gmsh.view.add("Hanging Node Conf")

    # gmsh.option.setNumber("Geometry.Points", 1)
    # gmsh.option.setNumber("Geometry.PointSize", 20)

    # gmsh.option.setNumber("Mesh.Points", 1)
    # gmsh.option.setNumber("Mesh.PointSize", 10)

    gmsh.view.option.setNumber(t2, "PointType", 2)
    gmsh.view.option.setNumber(t2, "PointSize", 10)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = folderName + fileNamePrefix + ".msh"
    gmsh.write( regMeshFileName )

    # gmsh.graphics.draw()
    # gmsh.model.setColor( [(0, ptVal) for ptVal in ptsVal], 2, 2, 127)  # Gray50

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.fltk.initialize()

    pngFileName = folderName + fileNamePrefix + ".png"
    gmsh.write( pngFileName )

    vtuFileName = fileNamePrefix
    meshFileUtils.createGMSHVTU( folderName, vtuFileName )

    gmsh.finalize()

def plotDFStats( plotFolderName, minVals, maxVals, aveVals, tagVals, formatVal ):

    for idx, tagVal in enumerate( tagVals ):

        plt.figure()
        plt.plot( minVals[idx], ".", label = "Minimum Values of " + tagVal )
        plt.legend()
        plt.xlabel( "Hanging Node Cases" )
        plt.ylabel( tagVal )
        plt.tight_layout()
        plt.grid(linestyle='dotted')
        plotfilename = plotFolderName + "Minimum" + tagVal + "." + formatVal
        plt.savefig( plotfilename, format = formatVal )

        plt.figure()
        plt.plot( maxVals[idx], ".", label = "Maximum Values of " + tagVal )
        plt.legend()
        plt.xlabel( "Hanging Node Cases" )
        plt.ylabel( tagVal )
        plt.tight_layout()
        plt.grid(linestyle = 'dotted')
        plotfilename = plotFolderName + "Maximum" + tagVal + "." + formatVal
        plt.savefig( plotfilename, format = formatVal )

        plt.figure()
        plt.plot( aveVals[idx], ".", label = "Average Values of " + tagVal )
        plt.legend()
        plt.xlabel( "Hanging Node Cases" )
        plt.ylabel( tagVal )
        plt.tight_layout()
        plt.grid(linestyle='dotted')
        plotfilename = plotFolderName + "Average" + tagVal + "." + formatVal
        plt.savefig( plotfilename, format = formatVal )

    plt.show()

def runAllPermutations(folderName, plotFolderName, nodeConfVal, algNumberDict, lcVals, *, showMesh = False ):

    allPermutes = []
    currentPermute = []
    allIsPresent = [[5, 7], [0, 1, 4, 5, 7], [5, 7], [0, 1, 4, 5, 7], [0, 1, 4, 5], [0, 1, 4, 5, 7], [5, 7], [0, 1, 4, 5, 7], [5, 7]]
    miscUtils.getPermutations(allIsPresent, 0, currentPermute, allPermutes)

    nIdx = 1
    minVals = []
    maxVals = []
    aveVals = []
    tagVals = []

    allEdgeConfs = set()

    for possiblePermute in allPermutes:

        if miscUtils.checkValidPermutation( possiblePermute ):

            if rotateReflect.checkRotationAndReflection( allEdgeConfs, possiblePermute ):
            # if True:

                isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in possiblePermute ] )
                allEdgeConfs.add( isPresentStr )
                print( nIdx, "Permute: ", possiblePermute )

                fileNameVal = "nodeConf" + str(nodeConfVal) + "_" + isPresentStr
                dfStats = gmshConfGenerate.runConfsAuto( folderName, fileNameVal, possiblePermute, lcVals, algNumberDict[ nodeConfVal ], 
                    printStatsFile = True, showMesh = showMesh, createVTU = True, saveMeshFile = True )

                if nIdx == 1:
                    tagVals = list( dfStats.columns )

                    for tagVal in tagVals:

                        minVals.append( [] )
                        maxVals.append( [] )
                        aveVals.append( [] )

                else:
                    for idx, tagVal in enumerate( tagVals ):

                        minVals[idx].append( dfStats.loc[ "min" ][ tagVal ] )
                        maxVals[idx].append( dfStats.loc[ "max" ][ tagVal ] )
                        aveVals[idx].append( dfStats.loc[ "mean" ][ tagVal ] )

                nIdx += 1

    formatVal = "png"
    print( nIdx )
    plotDFStats( plotFolderName, minVals, maxVals, aveVals, tagVals, formatVal )

if __name__ == "__main__":

    # folderName = "/home/gaurav/gmshAutoScripts/Images/HangingNodeConfs/"
    folderName = "/media/gaurav/easystore/HangingNodeConfs/"
    plotFolderName = "/media/gaurav/easystore/HangingNodePlots/"
    # runConfsTwoFacesHanging( folderName, fileNameVal, 3 )

    fileNameVal = "nodeConfAuto"

    lcVals = [0.5] * 27
    # lcVals = [1] * 27
    # lcVals[ 0 : 3 ] = [0.5] * 3
    # lcVals[ 24 : 27 ] = [0.5] * 3
    # lcVals[ 2:27:3 ] = [1]*9
    # lcVals[ 26 ] = 1

    isPresent = [ 5, 1, 5, 1, 1, 1, 5, 1, 5 ] # One Face Hanging
    # isPresent = [ 7, 1, 5, 7, 1, 1, 7, 1, 5 ] # Two Faces hanging
    # isPresent = [ 7, 5, 5, 7, 5, 5, 7, 5, 5] # Three faces hanging
    # isPresent = [ 7, 7, 7, 7, 1, 1, 7, 7, 7 ] # Four Faces hanging
    # isPresent = [ 7, 7, 7, 7, 1, 7, 7, 7, 7 ] # Five Faces hanging

    # isPresent = [5, 0, 5, 1, 0, 0, 5, 0, 5] # One Edge Hanging
    # isPresent = [ 7, 0, 5, 0, 0, 0, 5, 0, 7 ] # Two Edges hanging
    # isPresent = [ 7, 0, 7, 0, 0, 0, 5, 0, 7 ] # Three Edges hanging

    # isPresent = [ 7, 7, 7, 0, 0, 0, 5, 4, 5 ] # One face and one edge hanging
    # isPresent = [ 7, 7, 7, 0, 0, 0, 5, 5, 5 ] # One face and two edges hanging
    # isPresent = [ 7, 7, 7, 1, 0, 0, 5, 5, 5 ] # One face and three edges hanging

    nodeConfVal = 2
    algNumberDict = dict( [(1, 5), (2, 3)] )
    # isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in isPresent ] )
    # fileNameVal = "nodeConf" + str(nodeConfVal) + "_" + isPresentStr
    # gmshConfGenerate.runConfsAuto( folderName, fileNameVal, isPresent, lcVals, algNumberDict[ nodeConfVal ] )

    runAllPermutations(folderName, plotFolderName, nodeConfVal, algNumberDict, lcVals, showMesh = False )

    # performRotation( folderName, isPresent, lcVals, nodeConfVal, algNumberDict, "z" )
    # checkNonDuplicateEdgeConfOnRotation( [], [5, 1, 5, 1, 1, 1, 7, 5, 5], printConfs = True )
    
    # performReflection( [5, 1, 5, 1, 1, 1, 7, 5, 5], lcVals = lcVals, folderName = folderName, nodeConfVal = 2,\
                    #    algNumberDict = algNumberDict, dxn = "x", createConfs = True )
    