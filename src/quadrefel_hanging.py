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
    


if __name__ == "__main__":

    runconf()