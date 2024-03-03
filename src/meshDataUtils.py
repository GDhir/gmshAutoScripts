import numpy as np

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

def getElementNodeIndices( gmshFileName, types ):

    # if regexVal == "regular":
    #     types = ["3"]

    # if re.search( "triangle", regexVal ):
    #     types = ["2"]

    # if re.search( "hanging", regexVal ):
    #     types = [ "2", "3" ]

    # else:
    #     types = ["2", "3"]

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

def getTetrahedronEdgeLength( elementNodeIndices, allNodes ):

    edgeNodePairs = []

    pass

def getPyramidEdgeLength( elementNodeIndices, allNodes ):

    pass

def getElementEdgeLengths( elementNodalMap, allNodes, elementTypes ):

    pass