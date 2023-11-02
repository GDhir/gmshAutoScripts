import numpy as np
from scipy.linalg import solve

def runconf1():

    A = np.zeros( [5, 5] )

    hangingNodeLocs = [ (0, -1), (1, 0), (0, 1), (-1, 0) ]
    basisvals = [ "zetavals[idx]**2", "etavals[idx]**2",  "( zetavals[idx]**2 )*( etavals[idx] )",  "zetavals[idx]*( etavals[idx]**2 )" ]

    for locIdx in range(4):

        zetavals = [ -1, 1, 1, -1 ]
        etavals = [ -1, -1, 1, 1 ]

        zetavals.insert( locIdx + 1, hangingNodeLocs[locIdx][0] )
        etavals.insert( locIdx + 1, hangingNodeLocs[locIdx][1] )

        print( "Starting Hanging Node Location " + str( hangingNodeLocs[locIdx] ) + "\n" )
        print( "**********************************************************************" )

        for basisval in basisvals:

            print("Starting Basis " + basisval + "\n" )

            for idx in range(5):
                A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], eval( basisval ) ]

            for idx in range(5):

                b = np.zeros( [5, 1] )
                b[idx, 0] = 1

                try:
                    c = solve( A, b )
                except np.linalg.LinAlgError:

                    print( "Configuration with hanging node loc " + str( hangingNodeLocs[locIdx] ) + " and basis value " + basisval + " is singular \n" )
                    break

                print(np.transpose(c))

            print("Ending Basis " + basisval + "\n" )

        print( "**********************************************************************" )
        print( "Ending Hanging Node Location " + str( hangingNodeLocs[locIdx] ) + "\n" )

    # A[ 0, : ] = [ 1, -1, -1, 1, 1 ]
    # A[ 1, : ] = [ 1, 0, -1, 0, 0 ]
    # A[ 2, : ] = [ 1, 1, -1, -1, 1 ]
    # A[ 3, : ] = [ 1, 1, 1, 1, 1 ]
    # A[ 4, : ] = [ 1, -1, 1, -1, 1 ]
    

def runconf2():

    A = np.zeros( [5, 5] )

    hangingNodeLocs = [ (0, -1), (1, 0), (0, 1), (-1, 0) ]
    basisvals = [ "zetavals[idx]**2", "etavals[idx]**2",  "( zetavals[idx]**2 )*( etavals[idx] )",  "zetavals[idx]*( etavals[idx]**2 )" ]

    for locIdx1 in range(4):
        for locIdx2 in range( locIdx1 + 1, 4 ):

            zetavals = [ -1, 1, 1, -1 ]
            etavals = [ -1, -1, 1, 1 ]

            zetavals.insert( locIdx + 1, hangingNodeLocs[locIdx][0] )
            etavals.insert( locIdx + 1, hangingNodeLocs[locIdx][1] )

            print( "Starting Hanging Node Location " + str( hangingNodeLocs[locIdx] ) + "\n" )
            print( "**********************************************************************" )

            for basisval in basisvals:

                print("Starting Basis " + basisval + "\n" )

                for idx in range(5):
                    A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], eval( basisval ) ]

                for idx in range(5):

                    b = np.zeros( [5, 1] )
                    b[idx, 0] = 1

                    try:
                        c = solve( A, b )
                    except np.linalg.LinAlgError:

                        print( "Configuration with hanging node loc " + str( hangingNodeLocs[locIdx] ) + " and basis value " + basisval + " is singular \n" )
                        break

                    print(np.transpose(c))

                print("Ending Basis " + basisval + "\n" )

            print( "**********************************************************************" )
            print( "Ending Hanging Node Location " + str( hangingNodeLocs[locIdx] ) + "\n" )

        # A[ 0, : ] = [ 1, -1, -1, 1, 1 ]
        # A[ 1, : ] = [ 1, 0, -1, 0, 0 ]
        # A[ 2, : ] = [ 1, 1, -1, -1, 1 ]
        # A[ 3, : ] = [ 1, 1, 1, 1, 1 ]
        # A[ 4, : ] = [ 1, -1, 1, -1, 1 ]
    
def runconfLinear( quad = True ):

    allvals = []

    if quad:

        A = np.zeros( [4, 4] )

        zetavals = [ -1, 1, -1, 1 ]
        etavals = [ -1, -1, 1, 1 ]

        print( "**********************************************************************" )

        for idx in range(4):
            A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx] ]

        for idx in range(4):

            b = np.zeros( [4, 1] )
            b[idx, 0] = 1

            try:
                c = solve( A, b )
            except np.linalg.LinAlgError:

                print( "Configuration is singular \n" )
                break

            print(np.transpose(c))
            c = np.transpose(c)[0]
            allvals.append( c )

            print( "**********************************************************************" )
    else:

        A = np.zeros( [3, 3] )

        zetavals = [ -1, 1, -1 ]
        etavals = [ -1, -1, 1 ]

        print( "**********************************************************************" )

        for idx in range(3):
            A[idx, :] = [ 1, zetavals[idx], etavals[idx] ]

        for idx in range(3):

            b = np.zeros( [3, 1] )
            b[idx, 0] = 1

            try:
                c = solve( A, b )
            except np.linalg.LinAlgError:

                print( "Configuration is singular \n" )
                break

            print(np.transpose(c))
            c = np.transpose(c)[0]
            allvals.append( c )

            print( "**********************************************************************" )

    return np.array( allvals )

def getValueQuad( coords, coeffMat ):

    allBasisVals = np.zeros( (4, 4) )

    for cidx, coord in enumerate( coords ):
        
        allBasisVals[ 0, cidx ] = 1
        allBasisVals[ 1, cidx ] = coord[0]
        allBasisVals[ 2, cidx ] = coord[1]
        allBasisVals[ 3, cidx ] = coord[0]*coord[1]

    evaluations = np.matmul( coeffMat, allBasisVals )

    return evaluations

def getValueTri( coords, coeffMat ):

    allBasisVals = np.zeros( (3, len( coords )) )

    for cidx, coord in enumerate( coords ):
        
        allBasisVals[ 0, cidx ] = 1
        allBasisVals[ 1, cidx ] = coord[0]
        allBasisVals[ 2, cidx ] = coord[1]

    evaluations = np.matmul( coeffMat, allBasisVals )

    return evaluations

def getDerivativeQuad( coords, coeffMat ):

    allBasisVals = np.zeros( (2, 4) )  

    for cidx, coord in enumerate( coords ):
        
        allBasisVals[ 0, cidx ] = 1
        allBasisVals[ 1, cidx ] = coord[1]

    coeffVals = np.transpose( np.array( [ coeffMat[ :, 1 ], coeffMat[ :, 3 ] ] ) ) 
    derivZeta = np.transpose( np.matmul( coeffVals, allBasisVals ) )

    for cidx, coord in enumerate( coords ):
        
        allBasisVals[ 0, cidx ] = 1
        allBasisVals[ 1, cidx ] = coord[0]

    coeffVals = np.transpose( np.array( [ coeffMat[ :, 2 ], coeffMat[ :, 3 ] ] ) ) 
    derivEta = np.transpose( np.matmul( coeffVals, allBasisVals ) )

    return [derivZeta, derivEta] 

if __name__ == "__main__":

    # coeffMat = runconfLinear()

    # coordsQuadrature = [ [ -0.5773, -0.5773 ], [ 0.5773, -0.5773 ], [ -0.5773, 0.5773 ], [ 0.5773, 0.5773 ] ]
    # coordsVertices = [ [ -1, -1 ], [ 1, -1 ], [ -1, 1 ], [ 1, 1 ] ]
    # evaluations = getValueQuad( coordsVertices, coeffMat )
    # print(evaluations)
    # derivZeta, derivEta = getDerivativeQuad( coordsVertices, coeffMat )

    # print( derivZeta )

    # print( derivEta )

    # coeffMatTri = runconfLinear( False )
    # # coordsQuadrature = [ [ -0.6667, -0.6667 ], [ 0.3333, -0.6667 ], [ -0.6667, 0.3333 ] ]
    # # coordsQuadrature = [ [ -0.3333, -0.3333 ], [ -0.6, -0.6 ], [ -0.6, 0.2 ], [0.2, -0.6] ]
    # coordsQuadrature = [ [-0.8168, -0.8168], [0.63369, -0.8168], [-0.8168, 0.63369],
    #                      [-0.1081, -0.1081], [-0.78379, -0.1081], [-0.1081, -0.78379] ]
    # coordsVertices = [ [ -1, -1 ], [ 1, -1 ], [ -1, 1 ] ]
    
    # print( coeffMatTri )
    # evaluations = getValueTri( coordsQuadrature, coeffMatTri )
    # print( evaluations )

    pts = np.array( [
      (0.333333333333334,0.333333333333334),
      (0.470142064105115,0.470142064105115),
      (0.470142064105115,0.059715871789770),
      (0.059715871789770,0.470142064105115),
      (0.101286507323456,0.101286507323456),
      (0.101286507323456,0.797426985353087),
      (0.797426985353087,0.101286507323456)] )

    weights = np.array( [ 0.450000000000000/2, 0.264788305577012/2, 0.264788305577012/2,
               0.264788305577012/2, 0.251878361089654/2, 0.251878361089654/2, 0.251878361089654/2 ] )
    
    ptsref = 2*pts - 1

    print(ptsref)
    print(weights)