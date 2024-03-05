import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh
import miscUtils
import numpy as np
import gmshUtils
import meshFileUtils
import pandas as pd

def runConfsAuto( folderName, fileNamePrefix, isPresent, lcVals, algNumber = 3, *, printStatsFile = False, showMesh = False,
                createVTU = False, saveMeshFile = False ):

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    isPresentBitVals = [0]*27

    curIdx = 0
    allPts = [-1] * 27

    for edgePointConf in isPresent:

        pointConf = miscUtils.bitReprBase( edgePointConf )
        print(pointConf)

        for isPointPresent in pointConf:

            # print(curIdx)
            # nodeBitVals = miscUtils.bitReprBase( curIdx, 3 )

            # print(nodeBitVals)

            if isPointPresent:

                isPresentBitVals[curIdx] = 1
                nodeBitVals = miscUtils.bitReprBase( curIdx, 3 )

                # print(nodeBitVals)

                nodeCoords = miscUtils.getCoords( nodeBitVals )

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

        pointConf = miscUtils.bitReprBase( edgePointConf )

        for isPointPresent in pointConf:

            nodeBitVals = miscUtils.bitReprBase( curIdx, 3 )

            linesXIdx = linesXIdx + gmshUtils.addLine( nodeBitVals, isPointPresent, curIdx, linesXIdx, isPresentBitVals, allLinesZ, allPts, ptEdgeMap, 0 )
            linesYIdx = linesYIdx + gmshUtils.addLine( nodeBitVals, isPointPresent, curIdx, linesYIdx, isPresentBitVals, allLinesY, allPts, ptEdgeMap, 1 )
            linesZIdx = linesZIdx + gmshUtils.addLine( nodeBitVals, isPointPresent, curIdx, linesZIdx, isPresentBitVals, allLinesX, allPts, ptEdgeMap, 2 )

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

    if saveMeshFile:
        regMeshFileName = folderName + fileNamePrefix + ".msh"
        gmsh.write( regMeshFileName )

    gmsh.graphics.draw()
    gmsh.model.setColor( [(0, ptVal) for ptVal in allPts], 2, 2, 127)  # Gray50
    
    if showMesh:
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

    if createVTU:
        vtuFileName = fileNamePrefix
        meshFileUtils.createGMSHVTU( folderName, vtuFileName )

    if printStatsFile:
        _, eleTags , _ = gmsh.model.mesh.getElements(dim=3)
        gammaVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "gamma") )
        sicnVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "minSICN") )
        sigeVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "minSIGE") )

        df = pd.DataFrame( 
            {
                "Gamma": pd.Series( gammaVals ),
                "SICN": pd.Series( sicnVals ),
                "SIGE": pd.Series( sigeVals )
            }
        )

        dfStats = df.describe()
        statsFileName = folderName + fileNamePrefix + ".csv"
        dfStats.to_csv( statsFileName )

        return dfStats
    
    # angleVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "angleShape") )
    # minAngle = min( angleVals )
    # aveAngle = np.average( angleVals )
    # maxAngle = max( angleVals )    

    # minEdgeVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "minEdge") )
    # maxEdgeVals = np.array( gmsh.model.mesh.getElementQualities(eleTags[0], "maxEdge") )

    # edgeRatioVals = maxEdgeVals / minEdgeVals
    # minEdgeRatio = min( edgeRatioVals )
    # aveEdgeRatio = np.average( edgeRatioVals )
    # maxEdgeRatio = max( edgeRatioVals )

    # print(minEdgeRatio)
    # print(maxEdgeRatio)

    # resultVals = zip( eleTags[0], q )
    # print( list( resultVals ) )

    # gmsh.plugin.setNumber("AnalyseMeshQuality", "ICNMeasure", 1.)
    # gmsh.plugin.setNumber("AnalyseMeshQuality", "HidingThreshold", 0.3)
    # gmsh.plugin.setNumber("AnalyseMeshQuality", "ThresholdGreater",  1)

    # gmsh.plugin.setNumber("AnalyseMeshQuality", "CreateView", 1.)
    # t = gmsh.plugin.run("AnalyseMeshQuality")
    # dataType, tags, data, time, numComp = gmsh.view.getModelData(t, 0)
    # print('ICN for element {0} = {1}'.format(tags[0], data[0]))

