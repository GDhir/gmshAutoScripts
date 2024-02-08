def getPermutations( allVals, idx, currentPermute, allPermutes):

    if idx < len( allVals ) - 1:

        for curVal in allVals[idx]:

            currentPermute.append( curVal )
            getPermutations( allVals, idx + 1, currentPermute, allPermutes )
            currentPermute.pop()

    else:
        for curVal in allVals[ idx ]:

            currentPermute.append( curVal )
            allPermutes.append( [ curNewVal for curNewVal in currentPermute ] )
            currentPermute.pop()

def checkValidPermutation( permuteVals ):

    faceIdxVals = [1, 3, 5, 7]
    edgesToCheckVals = [ (0, 2), (0, 6), (2, 8), (6, 8) ]

    edgeChecksFor4 = [1, 3, 5, 7]

    for faceIdxVal, edgesToCheck in zip( faceIdxVals, edgesToCheckVals ):
        
        if permuteVals[ faceIdxVal ] == 7:
            if permuteVals[ edgesToCheck[0] ] != 7 or permuteVals[ edgesToCheck[1] ] != 7:

                return False
            
    if permuteVals[4] % 2 == 1:
        for edgeVal in edgeChecksFor4:

            if permuteVals[edgeVal] %2 != 1:
                return False

    if ( permuteVals[4] >> 2 ) == 1:
        for edgeVal in edgeChecksFor4:

            if ( permuteVals[edgeVal] >> 2 )%2 != 1:
                return False

    return True            
