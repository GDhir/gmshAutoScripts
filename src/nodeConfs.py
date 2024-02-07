import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh
import math
import os
import sys
import gmshUtils
import meshio
import numpy as np

def createGMSHVTU( folderName, mshFileName ):

    mesh = meshio.read( folderName + mshFileName + ".msh" )
    meshio.write( folderName + mshFileName + ".vtu", mesh )

    return 0

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
    createGMSHVTU( folderName, vtuFileName )

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
    createGMSHVTU( folderName, vtuFileName )

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
    createGMSHVTU( folderName, vtuFileName )

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
    createGMSHVTU( folderName, vtuFileName )

    gmsh.finalize()

def bitReprBase( number, baseVal = 2 ):

    output = [0, 0, 0]

    reprVal = np.base_repr( number, baseVal )

    for idx, x in enumerate( reprVal[::-1] ):
        output[idx] = int(x)

    # print(output)
    return output

def getCoords( nodeBitVals ):

    nodeCoords = [0, 0, 0]

    for idx, bitVal in enumerate( nodeBitVals ):
        nodeCoords[ idx ] = bitVal * 0.5

    # print( nodeBitVals )
    # print( nodeCoords )

    return nodeCoords

def addLine( nodeBitVals, isPointPresent, curIdx, linesIdx, isPresentBitVals, allLines, allPts, ptEdgeMap, dxn ):

    offset = 3**dxn

    if nodeBitVals[ dxn ] != 2:
        if isPointPresent:    

            nextPoint = curIdx + offset
            
            if isPresentBitVals[ nextPoint ]:
                
                allLines[ linesIdx ] = gmsh.model.occ.addLine( allPts[ curIdx ], allPts[ nextPoint ] )
                
                ptEdgeMap[ (curIdx, dxn, 0) ] = allLines[ linesIdx ]
                ptEdgeMap[ (nextPoint, dxn, 1) ] = -allLines[ linesIdx ]

            elif nodeBitVals[ dxn ] == 0 and isPresentBitVals[ nextPoint + offset ] and 1 not in nodeBitVals:
                allLines[ linesIdx ] = gmsh.model.occ.addLine( allPts[ curIdx ], allPts[ nextPoint + offset ] )

                ptEdgeMap[ (curIdx, dxn, 0) ] = allLines[ linesIdx ]
                ptEdgeMap[ (nextPoint + offset, dxn, 1) ] = -allLines[ linesIdx ]

        return 1
    
    else:
        return 0

def runConfsAuto( folderName, fileNamePrefix, isPresent, lcVals, algNumber = 3 ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    isPresentBitVals = [0]*27

    curIdx = 0
    allPts = [-1] * 27

    for edgePointConf in isPresent:

        pointConf = bitReprBase( edgePointConf )
        print(pointConf)

        for isPointPresent in pointConf:

            # print(curIdx)
            # nodeBitVals = bitReprBase( curIdx, 3 )

            # print(nodeBitVals)

            if isPointPresent:

                isPresentBitVals[curIdx] = 1
                nodeBitVals = bitReprBase( curIdx, 3 )

                # print(nodeBitVals)

                nodeCoords = getCoords( nodeBitVals )

                allPts[ curIdx ] = gmsh.model.occ.addPoint( nodeCoords[0], nodeCoords[1], nodeCoords[2], lcVals[ curIdx ] )

            curIdx = curIdx + 1

    curIdx = 0
    allLinesX = [-1] * 18
    allLinesY = [-1] * 18
    allLinesZ = [-1] * 18

    linesXIdx = 0
    linesYIdx = 0
    linesZIdx = 0

    ptEdgeMap = dict()

    for edgePointConf in isPresent:

        pointConf = bitReprBase( edgePointConf )

        for isPointPresent in pointConf:

            nodeBitVals = bitReprBase( curIdx, 3 )

            linesXIdx = linesXIdx + addLine( nodeBitVals, isPointPresent, curIdx, linesXIdx, isPresentBitVals, allLinesZ, allPts, ptEdgeMap, 0 )
            linesYIdx = linesYIdx + addLine( nodeBitVals, isPointPresent, curIdx, linesYIdx, isPresentBitVals, allLinesY, allPts, ptEdgeMap, 1 )
            linesZIdx = linesZIdx + addLine( nodeBitVals, isPointPresent, curIdx, linesZIdx, isPresentBitVals, allLinesX, allPts, ptEdgeMap, 2 )

            curIdx = curIdx + 1

    facePtStarts = [0, 2]
    surfaceVals = []

    # gmsh.model.occ.synchronize()
    
    transfiniteSurfaces = []

    for dxn in range(3):

        surfaceDxns = [ ( dxn + 1 ) % 3, ( dxn + 2 ) % 3 ]
        offsetVal = 3 ** ( ( dxn + 1 ) % 3 )

        for faceNum in range( 2 ):

            startPt = facePtStarts[ faceNum ] * ( 3 ** dxn )      
            isFull = True

            for ptIdx in range(9):

                curPtIdx = ( startPt + ptIdx * offsetVal ) % 26

                if not isPresentBitVals[ curPtIdx ]:

                    isFull = False
                    break

            if isFull:

                startPts = ( np.array( [0, 1, 3, 4], dtype = int ) * offsetVal + startPt ) % 26

                for semiFaceNum in range(4):

                    curStartPt = startPts[ semiFaceNum ]
                    # curStartPtVal = allPts[ curStartPt ]
                    curPtIdx = int( curStartPt )
                    lineIdx = 0

                    surfaceEdges = []

                    while 1:

                        # curPtVal = allPts[ curPtIdx ]
                        isNegativeDxn = lineIdx // 2
                        surfaceDxnIdx = lineIdx % 2
                        surfaceDxn = surfaceDxns[ surfaceDxnIdx ]

                        surfaceEdges.append( ptEdgeMap[ ( curPtIdx, surfaceDxn, isNegativeDxn ) ] )
                        
                        curoffsetVal =  3 ** ( surfaceDxn )
                        curPtIdx = int ( ( curPtIdx + ( 1 - 2 * isNegativeDxn ) * curoffsetVal ) )
                        lineIdx = lineIdx + 1

                        if curPtIdx == curStartPt:
                            break

                    cl = gmsh.model.occ.addCurveLoop( surfaceEdges )
                    surfaceVals.append( gmsh.model.occ.addPlaneSurface( [ cl ] ) )

                    gmsh.model.occ.synchronize()
                    # gmshUtils.setTransfiniteSurfaces( [surfaceVals[-1]], [], occ = True )

                    transfiniteSurfaces.append( surfaceVals[-1] )

            else:

                curPtIdx = int( startPt )
                lineIdx = 0
                surfaceEdges = []                

                while 1:

                    isNegativeDxn = lineIdx // 4
                    surfaceDxnIdx = ( lineIdx // 2 ) % 2
                    surfaceDxn = surfaceDxns[ surfaceDxnIdx ]

                    curoffsetVal =  3 ** ( surfaceDxn )

                    if isPresentBitVals[ curPtIdx ]:
                        surfaceEdges.append( ptEdgeMap[ ( curPtIdx, surfaceDxn, isNegativeDxn ) ] )

                    curPtIdx = int ( ( curPtIdx + ( 1 - 2 * isNegativeDxn ) * curoffsetVal ) )
                    lineIdx = lineIdx + 1

                    if curPtIdx == startPt:
                        break

                cl = gmsh.model.occ.addCurveLoop( surfaceEdges )
                surfaceVals.append( gmsh.model.occ.addPlaneSurface( [ cl ] ) )

                if len( surfaceEdges ) == 4:
                    transfiniteSurfaces.append( surfaceVals[-1] )

    gmsh.model.occ.synchronize()
    gmshUtils.setTransfiniteCurves( [allLinesX, allLinesY, allLinesZ], 2, occ = True )
    # transfiniteSurfaces.append(5)
    gmshUtils.setTransfiniteSurfaces( transfiniteSurfaces, [], occ = True )
    gmshUtils.recombineSurfaces( transfiniteSurfaces )

    sl = gmsh.model.occ.addSurfaceLoop( surfaceVals )
    vl = gmsh.model.occ.addVolume( [sl] )
    gmsh.model.occ.synchronize()

    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
    gmsh.option.setNumber("Mesh.Algorithm", algNumber)
    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = folderName + fileNamePrefix + ".msh"
    gmsh.write( regMeshFileName )

    gmsh.graphics.draw()
    gmsh.model.setColor( [(0, ptVal) for ptVal in allPts], 2, 2, 127)  # Gray50
    
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    vtuFileName = fileNamePrefix
    createGMSHVTU( folderName, vtuFileName )


if __name__ == "__main__":

    folderName = "/home/gaurav/gmshAutoScripts/Images/HangingNodeConfs/"
    # runConfsTwoFacesHanging( folderName, fileNameVal, 3 )

    fileNameVal = "nodeConfAuto"

    lcVals = [0.5] * 27
    # lcVals = [1] * 27
    # lcVals[ 0 : 3 ] = [0.5] * 3
    # lcVals[ 24 : 27 ] = [0.5] * 3
    # lcVals[ 2:27:3 ] = [1]*9
    # lcVals[ 26 ] = 1

    # isPresent = [ 5, 1, 5, 1, 1, 1, 5, 1, 5 ] # One Face Hanging
    # isPresent = [ 7, 1, 5, 7, 1, 1, 7, 1, 5 ] # Two Faces hanging
    # isPresent = [ 7, 5, 5, 7, 5, 5, 7, 5, 5] # Three faces hanging
    # isPresent = [ 7, 7, 7, 7, 1, 1, 7, 7, 7 ] # Four Faces hanging
    # isPresent = [ 7, 7, 7, 7, 1, 7, 7, 7, 7 ] # Five Faces hanging

    # isPresent = [ 7, 0, 5, 0, 0, 0, 5, 0, 7 ] # Two Edges hanging
    # isPresent = [ 7, 0, 7, 0, 0, 0, 5, 0, 7 ] # Three Edges hanging
    
    isPresent = [ 7, 7, 7, 0, 0, 0, 5, 4, 5 ] # One face and one edge hanging

    isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in isPresent ] )

    nodeConfVal = 2
    algNumberDict = dict( [(1, 5), (2, 3)] )

    fileNameVal = "nodeConf" + str(nodeConfVal) + "_" + isPresentStr
    runConfsAuto( folderName, fileNameVal, isPresent, lcVals, algNumberDict[ nodeConfVal ] )