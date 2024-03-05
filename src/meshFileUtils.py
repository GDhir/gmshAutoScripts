import re
from os import listdir
from os.path import isfile, join
import meshio

def getLevelFromFileName( meshval ):

    rval = re.search( "lvl", meshval )

    if rval:
        offset = rval.start() + len( "lvl" )
        rval = re.search( "[0-9]+", meshval[offset:] )
        rval = rval.group(0)

    return rval

def getAllLevels( meshFileNameArr ):

    levelsSet = set()

    for fileName in meshFileNameArr:

        level = getLevelFromFileName( fileName )
        levelsSet.add( int(level) )

    return list( levelsSet )

def getMeshFilesFromFolder( meshpath ):
    
    meshvals = [ join( meshpath, f ) for f in listdir(meshpath) if isfile(join(meshpath, f))]
    return meshvals

def getSortedMeshVals( meshvals, regexVal, regexCriterias ):

    keyVals = []
    prevIndexMap = dict()

    for index, meshval in enumerate(meshvals):

        if re.search( regexVal, meshval  ):

            regexCriteriaVals = []

            for regexCriteria in regexCriterias:
                rval = re.search( regexCriteria, meshval )
                offset = rval.start() + len( regexCriteria )
                rval = re.search( "[0-9]+", meshval[offset:] )
                rval = rval.group(0)

                regexCriteriaVals.append( rval )

            keyVal = 1

            for idx, regexCriteria in enumerate( regexCriterias ):
                regexCriteriaVal = int( regexCriteriaVals[ idx ] )
                keyVal = keyVal * regexCriteriaVal

            keyVals.append( keyVal )
            prevIndexMap[ keyVal ] = index            

    keyVals = sorted( keyVals )
    sortedMeshVals = []

    for keyVal in keyVals:

        idx = prevIndexMap[keyVal]
        sortedMeshVals.append( meshvals[idx] )

    return sortedMeshVals

def createGMSHVTU( folderName, mshFileName ):

    mesh = meshio.read( folderName + mshFileName + ".msh" )
    meshio.write( folderName + mshFileName + ".vtu", mesh )

    return 0