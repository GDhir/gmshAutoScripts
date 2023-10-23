import numpy as np
import fileNameUtils
import folderUtils   
import meshio
import h5py
from math import sin, cos, pi
import errorUtils

def getDealiiData( fileName, format = "hdf" ):

    if format == "hdf":
        f = h5py.File( fileName, 'r')
        nodes = np.array( f["nodes"] )
        solution = np.array( f["solution"] )
    elif format == "vtu":
        mesh = meshio.read( fileName )
        nodes = mesh.points
        # print(nodes)
        solution = mesh.point_data['solution']

    return (nodes, solution)

def getFinchMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, juliaVarName ):

    nLevels = len( levelsArr )
    minVals = np.ones( nLevels )*100000
    maxVals = np.ones( nLevels )*-1

    for optionIdx, option in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = option

        for idx, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )

            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            textFile = fileNameUtils.getTextFileName( folderUtils.finchTextfoldername, pythonVarName, juliaVarName )

            errVals = getData( textFile )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )
            # minMaxRangeVals.append( (minVal, maxVal) )
    minMaxRangeVals = []

    for idx in range( nLevels ):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals

def getDealiiMinMaxRange( levelsArr, allParams, optionsParam, comparisonParam, negative = 1, pival = 2*pi ):

    nLevels = len( levelsArr )
    minVals = np.ones( nLevels )*100000
    maxVals = np.ones( nLevels )*-1

    dataFileFormat = "vtu"

    for optionIdx, option in enumerate( allParams[ comparisonParam ] ):

        optionsParam[ comparisonParam ] = option

        for idx, level in enumerate( levelsArr ):

            optionsParam[ "level" ] = str( level )

            pythonVarName = fileNameUtils.getPythonVarName( optionsParam )
            textFile = fileNameUtils.getTextFileName( folderUtils.dealiiTextfoldername, pythonVarName, "solutionvalues", dataFileFormat )

            nodes, solution = getDealiiData( textFile, dataFileFormat )
            errVals = errorUtils.getDealiiError( nodes, solution, negative, pival )

            minVals[idx] = np.min( [minVals[idx], np.min( errVals )] )
            maxVals[idx] = np.max( [maxVals[idx], np.max( errVals )] )

    minMaxRangeVals = []

    for idx in range( nLevels ):
        minMaxRangeVals.append( ( minVals[idx], maxVals[idx] ) )

    return minMaxRangeVals

def getData( filename ):

    uvals = []

    with open(filename) as fval:

        uvals = fval.readlines()

    # print(uvals)
    uvals = uvals[0].split(",")

    uvals[0] = uvals[0][1:]
    uvals[-1] = uvals[-1][:-1]

    for idx, val in enumerate(uvals):

        uvals[idx] = float( uvals[idx] )

    return uvals