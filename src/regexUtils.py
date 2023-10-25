import re

def getCriteriaValsString( regexCriterias, criteriaVals ):

    regexCriteriaVals = ""

    for idx, regexCriteria in enumerate( regexCriterias ):
        regexCriteriaVals += regexCriteria + criteriaVals[idx]

    return regexCriteriaVals 

def getRegexCriterias( regexVal ):

    # if regexVal == "hanging":
    #     regexCriterias = [ "Nx=", "Ny=" ]
    # elif regexVal == "regular" or re.search( "triangle", regexVal ):
    #     regexCriterias = [ "N=" ]
    # elif regexVal == "mesh":
    #     regexCriterias = [ "lvl" ]

    regexCriterias = ["lvl"]
    
    return regexCriterias

def getRegexVal( fileName, allRegexes ):

    for regexVal in allRegexes:
        if re.search( regexVal, fileName ):

            return regexVal