import numpy as np
from math import sin, cos, pi

def getExactSol( xval, yval, zval, negative = -1, pival = 2*pi, dim = 2 ):

    if dim == 1:
        return negative*sin(pival*xval)
    elif dim == 2:
        return negative*sin(pival*xval) * sin(pival*yval)
    elif dim == 3:
        return negative * sin(pival*xval) * sin(pival*yval) * sin(pival*zval)

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

def getFinchError( xvals, yvals = [], zvals = [], uvals = [], pival = 2*pi, dim = 2 ):

    errVals = []
    exactvals = []

    for idx, xval in enumerate( xvals ):

        yval = 0
        zval = 0
        
        if dim > 1:
            yval = yvals[idx]
        
        if dim > 2:
            zval = zvals[idx]
            
        solval = uvals[idx]
        # print(xval, yval)
        exactval = getExactSol( xval, yval, zval, 1, pival, dim )
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