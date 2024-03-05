import miscUtils
import gmshConfGenerate

def getRotatedIdx( base3Idx, dxn = "x" ):

    if dxn == "x":
        temp = base3Idx[1]
        base3Idx[1] = 2 - base3Idx[2]
        base3Idx[2] = temp

    elif dxn == "y":
        temp = base3Idx[0]
        base3Idx[0] = 2 - base3Idx[2]
        base3Idx[2] = temp     

    elif dxn == "z":
        temp = base3Idx[0]
        base3Idx[0] = 2 - base3Idx[1]
        base3Idx[1] = temp

    return base3Idx

def checkNonDuplicateEdgeConfOnRotation( allEdgeConfs, curEdgeConf, *, printConfs = False ):

    rotationIdx = 0

    for xRotation in range(4):

        if xRotation == 0:
            prevXEdgeConf = curEdgeConf
        else:
            prevXEdgeConf = xEdgeConf

        xEdgeConf = performRotation( prevXEdgeConf, dxn = "x" )
        rotationIdx += 1

        if printConfs:
            print( rotationIdx, xEdgeConf )

        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in xEdgeConf ] )
        if isPresentStr in allEdgeConfs:
            return False

        for yRotation in range(4):

            if yRotation == 0:
                prevYEdgeConf = xEdgeConf
            else:
                prevYEdgeConf = yEdgeConf

            yEdgeConf = performRotation( prevYEdgeConf, dxn = "y" )
            rotationIdx += 1

            if printConfs:
                print( rotationIdx, yEdgeConf )

            isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in yEdgeConf ] )
            if isPresentStr in allEdgeConfs:
                return False

            for zRotation in range(4):

                if zRotation == 0:
                    prevZEdgeConf = yEdgeConf
                else:
                    prevZEdgeConf = zEdgeConf

                zEdgeConf = performRotation( prevZEdgeConf, dxn = "z" )
                rotationIdx += 1

                if printConfs:
                    print( rotationIdx, zEdgeConf )

                isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in zEdgeConf ] )
                if isPresentStr in allEdgeConfs:
                    return False

    return True

def performRotation( edgePointConf, *, lcVals = [], folderName = "", nodeConfVal = 2, algNumberDict = dict(),\
                    dxn = "x", createConfs = False ):

    if createConfs:
        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in edgePointConf ] )
        fileNameVal = "originalNodeConf" + str(nodeConfVal) + "_" + isPresentStr
        gmshConfGenerate.runConfsAuto( folderName, fileNameVal, edgePointConf, lcVals, algNumberDict[ nodeConfVal ], \
                    printStatsFile = False, showMesh = True, createVTU = False, saveMeshFile = False  )

    isPresentBitVals = convertEdgePointConfToIsPresentBitVals( edgePointConf )
    rotatedIsPresentBitVals = getRotatedIsPresentBitVals( isPresentBitVals, dxn )
    newEdgePointConf = []

    for edgePointStartIdx in range( 0, 27, 3 ):

        edgeConf = rotatedIsPresentBitVals[ edgePointStartIdx ] + rotatedIsPresentBitVals[ edgePointStartIdx + 1 ] * ( 2 ** 1 ) + \
             rotatedIsPresentBitVals[ edgePointStartIdx + 2 ] * ( 2 ** 2 )
        
        newEdgePointConf.append( edgeConf )
    # print(rotatedIsPresentBitVals)

    if createConfs:
        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in newEdgePointConf ] )
        fileNameVal = "rotatedNodeConf" + str(nodeConfVal) + "_" + isPresentStr
        gmshConfGenerate.runConfsAuto( folderName, fileNameVal, newEdgePointConf, lcVals, algNumberDict[ nodeConfVal ], \
                    printStatsFile = False, showMesh = True, createVTU = False, saveMeshFile = False )
        
    return newEdgePointConf

def getRotatedIsPresentBitVals( isPresentBitVals, dxn ):

    rotatedIsPresentBitVals = [0]*27

    for idx, isPresentVal in enumerate( isPresentBitVals ):

        base3Idx = miscUtils.bitReprBase( idx, 3 )
        base3Idx = getRotatedIdx( base3Idx, dxn )

        newIdxVal = base3Idx[2] * (3 ** 2) + base3Idx[1] * (3 ** 1) + base3Idx[0]
        rotatedIsPresentBitVals[newIdxVal] = isPresentVal

    return rotatedIsPresentBitVals

def swapBitVals( isPresentBitVals, idxVal1, idxVal2 ):

    temp = isPresentBitVals[ idxVal1 ]
    isPresentBitVals[ idxVal1 ] = isPresentBitVals[ idxVal2 ]
    isPresentBitVals[ idxVal2 ] = temp

def convertEdgePointConfToIsPresentBitVals( edgePointConf ):

    isPresentBitVals = [0]*27
    curIdx = 0

    for edgeConf in edgePointConf:

        pointConf = miscUtils.bitReprBase( edgeConf )

        for isPointPresent in pointConf:

            if isPointPresent:
                isPresentBitVals[curIdx] = 1

            curIdx = curIdx + 1

    return isPresentBitVals

def performReflection( edgePointConf, *, lcVals = [], folderName = "", nodeConfVal = 2, algNumberDict = dict(),\
                    dxn = "x", createConfs = False ):

    if createConfs:
        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in edgePointConf ] )
        fileNameVal = "originalNodeConf" + str(nodeConfVal) + "_" + isPresentStr
        gmshConfGenerate.runConfsAuto( folderName, fileNameVal, edgePointConf, lcVals, algNumberDict[ nodeConfVal ], \
                    printStatsFile = False, showMesh = True, createVTU = False, saveMeshFile = False  )

    isPresentBitVals = convertEdgePointConfToIsPresentBitVals( edgePointConf )
    reflectedIsPresentBitVals = getReflectedIsPresentBitVals( isPresentBitVals, dxn )
    newEdgePointConf = []

    for edgePointStartIdx in range( 0, 27, 3 ):

        edgeConf = reflectedIsPresentBitVals[ edgePointStartIdx ] + reflectedIsPresentBitVals[ edgePointStartIdx + 1 ] * ( 2 ** 1 ) + \
             reflectedIsPresentBitVals[ edgePointStartIdx + 2 ] * ( 2 ** 2 )
        
        newEdgePointConf.append( edgeConf )
    # print(reflectedIsPresentBitVals)

    if createConfs:
        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in newEdgePointConf ] )
        fileNameVal = "reflectedNodeConf" + str(nodeConfVal) + "_" + isPresentStr
        gmshConfGenerate.runConfsAuto( folderName, fileNameVal, newEdgePointConf, lcVals, algNumberDict[ nodeConfVal ], \
                    printStatsFile = False, showMesh = True, createVTU = False, saveMeshFile = False )
        
    return newEdgePointConf

def getReflectedIsPresentBitVals( isPresentBitVals, dxn ):

    expDict = dict( [ ("x", [1, 2, 0]), ("y", [0, 2, 1]), ("z", [0, 1, 2]) ] )

    for idx1 in range(3):
        for idx2 in range(3):

            exp1 = 3 ** expDict[dxn][0]
            exp2 = 3 ** expDict[dxn][1]
            exp3 = 3 ** expDict[dxn][2]             

            idxVal1 = idx1 * ( exp1 ) + idx2 * ( exp2 )
            idxVal2 = idx1 * ( exp1 ) + idx2 * ( exp2 ) + 2 * ( exp3 )
            swapBitVals( isPresentBitVals, idxVal1, idxVal2 )     

    return isPresentBitVals

def checkRotationAndReflection( allEdgeConfs, curEdgeConf, *, printConfs = False ):

    reflectionIdx = 0

    for xReflection in range(2):

        if xReflection == 0:
            prevXEdgeConf = curEdgeConf
        else:
            prevXEdgeConf = xEdgeConf

        xEdgeConf = performReflection( prevXEdgeConf, dxn = "x" )
        reflectionIdx += 1

        if not checkNonDuplicateEdgeConfOnRotation( allEdgeConfs, xEdgeConf, printConfs = printConfs ):
            return False

        isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in xEdgeConf ] )
        if isPresentStr in allEdgeConfs:
            return False

        for yReflection in range(2):

            if yReflection == 0:
                prevYEdgeConf = xEdgeConf
            else:
                prevYEdgeConf = yEdgeConf

            yEdgeConf = performReflection( prevYEdgeConf, dxn = "y" )
            reflectionIdx += 1

            if not checkNonDuplicateEdgeConfOnRotation( allEdgeConfs, yEdgeConf, printConfs = printConfs ):
                return False
            
            isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in yEdgeConf ] )
            if isPresentStr in allEdgeConfs:
                return False

            for zReflection in range(2):

                if zReflection == 0:
                    prevZEdgeConf = yEdgeConf
                else:
                    prevZEdgeConf = zEdgeConf

                zEdgeConf = performReflection( prevZEdgeConf, dxn = "z" )
                reflectionIdx += 1

                if not checkNonDuplicateEdgeConfOnRotation( allEdgeConfs, zEdgeConf, printConfs = printConfs ):
                    return False

                isPresentStr = ''.join( [ str( isPresentVal ) for isPresentVal in zEdgeConf ] )
                if isPresentStr in allEdgeConfs:
                    return False

    return True