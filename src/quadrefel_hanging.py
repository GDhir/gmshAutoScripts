import numpy as np
from scipy.linalg import solve

def runconf1():

    A = np.zeros( [5, 5] )

    zetavals = [ -1, 0, 1, 1, -1 ]
    etavals = [ -1, -1, -1, 1, 1 ]

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]**2 ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], etavals[idx]**2 ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    print("end")

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], ( zetavals[idx]**2 )*( etavals[idx] ) ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]*( etavals[idx]**2 ) ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    # A[ 0, : ] = [ 1, -1, -1, 1, 1 ]
    # A[ 1, : ] = [ 1, 0, -1, 0, 0 ]
    # A[ 2, : ] = [ 1, 1, -1, -1, 1 ]
    # A[ 3, : ] = [ 1, 1, 1, 1, 1 ]
    # A[ 4, : ] = [ 1, -1, 1, -1, 1 ]
    

def runconf2():

    A = np.zeros( [5, 5] )

    zetavals = [ -1, 1, 1, 1, -1 ]
    etavals = [ -1, -1, 0, 1, 1 ]

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]**2 ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], etavals[idx]**2 ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    print("end")

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], ( zetavals[idx]**2 )*( etavals[idx] ) ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]*( etavals[idx]**2 ) ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

def runconf3():

    A = np.zeros( [5, 5] )

    zetavals = [ -1, 1, 1, 0, -1 ]
    etavals = [ -1, -1, 1, 1, 1 ]

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]**2 ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], etavals[idx]**2 ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    print("end")

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], ( zetavals[idx]**2 )*( etavals[idx] ) ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]*( etavals[idx]**2 ) ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

def runconf4():
    A = np.zeros( [5, 5] )

    zetavals = [ -1, 1, 1, -1, -1 ]
    etavals = [ -1, -1, 1, 1, 0 ]

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]**2 ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], etavals[idx]**2 ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))

    print("end")

    # for idx in range(5):
    #     A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], ( zetavals[idx]**2 )*( etavals[idx] ) ]

    # for idx in range(5):

    #     b = np.zeros( [5, 1] )
    #     b[idx, 0] = 1
    #     c = solve( A, b )

    #     print(np.transpose(c))

    for idx in range(5):
        A[idx, :] = [ 1, zetavals[idx], etavals[idx], zetavals[idx]*etavals[idx], zetavals[idx]*( etavals[idx]**2 ) ]

    for idx in range(5):

        b = np.zeros( [5, 1] )
        b[idx, 0] = 1
        c = solve( A, b )

        print(np.transpose(c))




if __name__ == "__main__":

    runconf1()
    print("end")
    runconf2()
    print("end")
    runconf3()
    print("end")
    runconf4()