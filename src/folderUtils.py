import os

finchRootfoldername = "/media/gaurav/easystore/Finch/MixedElement/"
dealiiRootfoldername = "/media/gaurav/easystore/dealii/MixedElement/"

finchTextfoldername = finchRootfoldername + "TextFiles/"


finchPlotfoldername = finchRootfoldername + "PlotFiles/SimPlots/"

dealiiTextfoldername = dealiiRootfoldername + "TextFiles/"
dealiiPlotfoldername = dealiiRootfoldername + "PlotFiles/SimPlots/"
gmshImageFolderName = "/home/gaurav/gmshAutoScripts/Images/"

textFolderNames = dict()
textFolderNames["Finch"] = finchTextfoldername
textFolderNames["Dealii"] = dealiiTextfoldername

plotFolderNames = dict()
plotFolderNames["Finch"] = finchPlotfoldername
plotFolderNames["Dealii"] = dealiiPlotfoldername

def checkAndCreateFolder( folderName ):

    if not os.path.exists( folderName ):

        folderSplit = folderName.split( "/" )

        folderStr = ""

        for folderName in folderSplit:

            if not folderName:
                continue
            
            folderStr = folderStr + "/" + folderName

            if not os.path.exists( folderStr ):
                os.mkdir( folderStr )        
