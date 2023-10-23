import re
import regexUtils

def getTextFileName( folderName, pythonVarName, juliaVarName, extension = "txt" ):

    filename = folderName + pythonVarName + juliaVarName + "." + extension

    return filename

def getCriteriasFromFileName( meshval, regexCriterias ):

    regexCriteriaVals = ""

    for regexCriteria in regexCriterias:
        rval = re.search( regexCriteria, meshval )

        if rval:
            offset = rval.start() + len( regexCriteria )
            rval = re.search( "[0-9]+", meshval[offset:] )
            rval = rval.group(0)

            regexCriteriaVals += regexCriteria + "=" + rval

    return regexCriteriaVals[:-1]

def getMeshFileName( optionsParam, criterias, criteriaVals, meshPath ):

    criteriaStr = regexUtils.getCriteriaValsString( criterias, criteriaVals )
    fileName = meshPath + optionsParam[ "meshRegexVal" ] + criteriaStr + ".msh"

    return fileName

def getFileNameFromMeshName( meshval, folderName, varName, regexTypeVal, regexCriterias, extension ):

    regexCriteriaVals = []

    for regexCriteria in regexCriterias:
        rval = re.search( regexCriteria, meshval )

        if rval:
            offset = rval.start() + len( regexCriteria )
            rval = re.search( "[0-9]+", meshval[offset:] )
            rval = rval.group(0)

            regexCriteriaVals.append( rval )

    fileName = folderName + regexTypeVal + "_" + varName + "_"

    for idx, regexCriteria in enumerate( regexCriterias ):

        regexCriteriaVal = regexCriteriaVals[ idx ]

        fileName = fileName + regexCriteria + "=" + regexCriteriaVal

    fileName = fileName + extension

    return fileName

def getPythonVarName( optionsParam ):

    pythonVarName = ""

    for param, paramValue in optionsParam.items():

        pythonVarName += param + "_" + paramValue

    return pythonVarName