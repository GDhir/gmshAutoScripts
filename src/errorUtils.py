import numpy as np
from math import sin, cos, pi

def getExactSol( xval, yval, negative = -1, pival = 2*pi ):

    return negative*sin(pival*xval)*sin(pival*yval)

def getDealiiError( nodes, solution, negative = 1, pival = 2*pi ):

    errVals = []
    exactvals = []

    for idx, node in enumerate( nodes ):

        xval = node[0]
        yval = node[1]
        # print(xval, yval)
        exactval = getExactSol( xval, yval, negative, pival )
        # exactvals.append( exactval )

        if len( solution.shape ) == 2:
            solval = solution[ idx, 0 ]
        elif len( solution.shape ) == 1:
            solval = solution[ idx ]

        errorval = np.abs( exactval - solval )
        errVals.append( errorval )

    # print(errVals)

    # exactvals = np.array(exactvals)

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], exactvals )

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], solution )
    # plt.colorbar()
    # plt.savefig( "solutionFile.png" )

    return np.array( errVals )

def getFinchError( xvals, yvals, uvals, pival = 2*pi ):

    errVals = []
    exactvals = []

    for idx, xval in enumerate( xvals ):

        yval = yvals[idx]
        solval = uvals[idx]
        # print(xval, yval)
        exactval = getExactSol( xval, yval, 1, pival )
        exactvals.append( exactval )

        errorval = np.abs( exactval - solval )
        errVals.append( errorval )

    # print(errVals)

    exactvals = np.array(exactvals)

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], exactvals )

    # plt.figure()
    # plt.tricontourf( nodes[:, 0], nodes[:, 1], solution[:, 0] )

    return np.array( errVals )