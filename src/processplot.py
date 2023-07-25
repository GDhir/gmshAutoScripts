import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import re
from os import listdir
from os.path import isfile, join
import os
import subprocess
from paraview.simple import *
import meshio

rootfoldername = "/media/gaurav/easystore/Finch/MixedElement/"
textfoldername = rootfoldername + "TextFiles/"
plotfoldername = rootfoldername + "PlotFiles/SimPlots/"
gmshImageFolderName = "/home/gaurav/gmshAutoScripts/Images/"

def getSortedMeshVals( meshvals, regexVal ):

    keyvals = []
    prevIndexMap = dict()

    if regexVal == "hanging":

        for index, meshval in enumerate(meshvals):

            if re.search( regexVal, meshval  ):
                Nxval = re.search( "Nx=", meshval )
                offset = Nxval.start() + 3
                Nxval = re.search( "[0-9]+", meshval[offset:] )
                Nxval = int( Nxval.group(0) )

                Nyval = re.search( "Ny=", meshval )
                offset = Nyval.start() + 3
                Nyval = re.search( "[0-9]+", meshval[offset:] )
                Nyval = int( Nyval.group(0) )

                keyval = Nxval*Nyval
                keyvals.append( keyval )

                prevIndexMap[ keyval ] = index

    if regexVal == "regular" or regexVal == "triangle":

        for index, meshval in enumerate(meshvals):

            if re.search( regexVal, meshval ):

                regularNval = re.search( "[0-9]+"  , meshval  )
                regularNval = regularNval.group(0)
                
                keyval = int(regularNval)
                keyvals.append( keyval )

                prevIndexMap[ keyval ] = index

    keyvals = sorted( keyvals )
    sortedMeshVals = []

    for keyval in keyvals:

        idx = prevIndexMap[keyval]
        sortedMeshVals.append( meshvals[idx] )

    return sortedMeshVals

def getTextFileName( meshval, varName ):

    if re.search( "hanging", meshval  ):
        Nxval = re.search( "Nx=", meshval )
        offset = Nxval.start() + 3
        Nxval = re.search( "[0-9]+", meshval[offset:] )
        Nxval = Nxval.group(0)

        Nyval = re.search( "Ny=", meshval )
        offset = Nyval.start() + 3
        Nyval = re.search( "[0-9]+", meshval[offset:] )
        Nyval = Nyval.group(0)

        textFileName = textfoldername + "mixed_" + varName + "values_Nx=" + Nxval + "Ny=" + Nyval + ".txt"

    if re.search( "regular", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        textFileName = textfoldername + "regular_" + varName + "values_N=" + regularNval + ".txt"

    if re.search( "triangle", meshval  ):
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        textFileName = textfoldername + "triangle_" + varName + "values_N=" + regularNval + ".txt"

    return textFileName

def getMinMaxRange( sortedMeshVals1, sortedMeshVals2 ):

    minMaxRangeVals = []
    meshValsLen = len( sortedMeshVals1 )

    for idx in range( meshValsLen ):

        meshval1 = sortedMeshVals1[idx]
        meshval2 = sortedMeshVals2[idx]

        textFile1 = getTextFileName( meshval1, "error" )
        textFile2 = getTextFileName( meshval2, "error" )

        errVals1 = getData( textFile1 )
        errVals2 = getData( textFile2 )

        minVal = np.min( [np.min( errVals1 ), np.min( errVals2 )] )
        maxVal = np.max( [np.max( errVals1 ), np.max( errVals2 )] )

        minMaxRangeVals.append( (minVal, maxVal) )

    return minMaxRangeVals


def getData( filename ):

    uvals = []

    with open(filename) as fval:

        uvals = fval.readlines()

    uvals = uvals[0].split(",")

    uvals[0] = uvals[0][1:]
    uvals[-1] = uvals[-1][:-1]

    for idx, val in enumerate(uvals):

        uvals[idx] = float( uvals[idx] )

    return uvals

# def fmt(x, pos):
#     a, b = '{:.2e}'.format(x).split('e')
#     b = int(b)
#     return r'${} \times 10^{{{}}}$'.format(a, b)

def buildMesh( gmshfilecmd, gmshfileargs ):

    gmshbuildFolder = "/home/gaurav/gmshAutoScripts/build/"
    compilecmd = ["make", "-j", "4", gmshfilecmd]
    subprocess.run( compilecmd, cwd = gmshbuildFolder )

    gmshfileargs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshlstcmd = [ gmshbuildFolder + gmshfilecmd, gmshfileargs ]
    subprocess.run( meshlstcmd, cwd = gmshbuildFolder )

def runJulia( exefilename ):

    juliapath = "/home/gaurav/Downloads/julia-1.9.2-linux-x86_64/julia-1.9.2/bin/julia"
    julialstcmd = [juliapath, exefilename]
    subprocess.run( julialstcmd )

def removeFiles( dirval ):

    meshvals = [join(dirval, f) for f in listdir(dirval) if isfile(join(dirval, f))]

    for meshfilename in meshvals:
        cmdlist = [ "rm", meshfilename ]
        subprocess.run( cmdlist )

def runSim( simPlotFolderName ):

    gmshfileargs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    removeFiles( gmshfileargs )

    gmshfilecmd = "hangingMeshv8"
    buildMesh( gmshfilecmd, gmshfileargs )

    gmshfilecmd = "regularMeshv3"
    buildMesh( gmshfilecmd, gmshfileargs )

    exefilename = "example-mixed-element-2d.jl"
    runJulia( exefilename )

    # showplot( simPlotFolderName )
    # showplotTriangle( simPlotFolderName )
    showParaviewPlot( simPlotFolderName )
    # removeFiles( gmshfileargs )

def showParaviewPlot( simPlotFolderName ):

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshVizPath = "/home/gaurav/Finch/src/examples/Mesh/MeshViz/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    sortedMeshValsHanging = getSortedMeshVals( meshvals, "hanging" )
    sortedMeshValsRegular = getSortedMeshVals( meshvals, "regular" )

    minMaxRangeVals = getMinMaxRange( sortedMeshValsHanging, sortedMeshValsRegular )

    for index, meshval in enumerate( sortedMeshValsHanging + sortedMeshValsRegular ):
        mesh = meshio.read( meshpath + meshval )
        meshio.write( meshVizPath + meshval[:-4] + ".vtu", mesh )

    for index, meshval in enumerate(sortedMeshValsHanging):

        meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
        Nxval = re.search( "Nx=", meshval )
        offset = Nxval.start() + 3
        Nxval = re.search( "[0-9]+", meshval[offset:] )
        Nxval = Nxval.group(0)

        Nyval = re.search( "Ny=", meshval )
        offset = Nyval.start() + 3
        Nyval = re.search( "[0-9]+", meshval[offset:] )
        Nyval = Nyval.group(0)

        hangingfilename = textfoldername + "mixed_" + "errorAnduNx=" + Nxval + "Ny=" + Nyval + ".vtu"

        # view = CreateRenderView( )
        solfile = OpenDataFile( hangingfilename )
        display = Show(solfile)
        ColorBy(display, ('POINTS', 'err'))
        minVal, maxVal = minMaxRangeVals[ index ]
        colorMap = GetColorTransferFunction('err')
        colorMap.RescaleTransferFunction( minVal, maxVal )        

        gmshfile = OpenDataFile( meshvtkName )
        dpGmsh = GetDisplayProperties( gmshfile )
        dpGmsh.Representation = 'Wireframe'
        gmshdisplay = Show(gmshfile)
        
        Render()
        
        myview = GetActiveView()
        myview.ViewSize = [1920, 1080]
        
        dpSol = GetDisplayProperties(solfile, myview)
        # to show the color legend
        dpSol.SetScalarBarVisibility(myview, True)

        curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"   
        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        SaveScreenshot( curPlotFolderName + "paraview_mixed_" + "errorAnduNx=" + Nxval + "Ny=" + Nyval + ".png", myview)

        Hide( solfile )
        Hide( gmshfile )

    for index, meshval in enumerate( sortedMeshValsRegular ):

        meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        regularfilename = textfoldername + "regular_" + "errorAnduN=" + regularNval + ".vtu"

        # view = CreateRenderView(  )
        solfile = OpenDataFile( regularfilename )
        display = Show(solfile)
        ColorBy(display, ('POINTS', 'err'))
        minVal, maxVal = minMaxRangeVals[ index ]
        colorMap = GetColorTransferFunction('err')
        colorMap.RescaleTransferFunction( minVal, maxVal )  

        gmshfile = OpenDataFile( meshvtkName )
        dpGmsh = GetDisplayProperties( gmshfile )
        dpGmsh.Representation = 'Wireframe'
        gmshdisplay = Show(gmshfile)
        
        Render()

        myview = GetActiveView()
        myview.ViewSize = [1920, 1080]

        dpSol = GetDisplayProperties(solfile, myview)
        # to show the color legend
        dpSol.SetScalarBarVisibility(myview, True)

        curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"   
        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        SaveScreenshot( curPlotFolderName + "paraview_regular_" + "errorAnduN=" + regularNval + ".png", myview)

        Hide( solfile )
        Hide( gmshfile )


    # if re.search( "triangle", meshval ):
    #     meshvtkName = meshVizPath + meshval[:-4] + ".vtu"
    #     regularNval = re.search( "[0-9]+"  , meshval  )
    #     regularNval = regularNval.group(0)

    #     regularfilename = textfoldername + "triangle_" + "errorAnduN=" + regularNval + ".vtu"

    #     # view = CreateRenderView(  )
    #     solfile = OpenDataFile( regularfilename )
    #     display = Show(solfile)
    #     ColorBy(display, ('POINTS', 'err'))

    #     gmshfile = OpenDataFile( meshvtkName )
    #     dpGmsh = GetDisplayProperties( gmshfile )
    #     dpGmsh.Representation = 'Wireframe'
    #     gmshdisplay = Show(gmshfile)

    #     Render()

    #     myview = GetActiveView()
    #     myview.ViewSize = [1920, 1080]

    #     dpSol = GetDisplayProperties(solfile, myview)
    #     # to show the color legend
    #     dpSol.SetScalarBarVisibility(myview, True)

    #     SaveScreenshot( plotfoldername + "triangle_" + "errorAnduN=" + regularNval + ".png", myview)    

    #     Hide( solfile )
    #     Hide( gmshfile )
    

def showplot( simPlotFolderName ):

    if not os.path.exists(simPlotFolderName):
        os.mkdir( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfileregular.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]

    sortedMeshValsRegular = getSortedMeshVals( meshvals, "regular" )
    sortedMeshValsHanging = getSortedMeshVals( meshvals, "hanging" )

    minMaxRangeVals = getMinMaxRange( sortedMeshValsRegular, sortedMeshValsHanging )



    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    hangingMaxErrorList = []
    regularMaxErrorList = []
    hangingl2ErrorList = []
    regularl2ErrorList = []

    for index, meshval in enumerate( sortedMeshValsHanging ):

        Nxval = re.search( "Nx=", meshval )
        offset = Nxval.start() + 3
        Nxval = re.search( "[0-9]+", meshval[offset:] )
        Nxval = Nxval.group(0)

        Nyval = re.search( "Ny=", meshval )
        offset = Nyval.start() + 3
        Nyval = re.search( "[0-9]+", meshval[offset:] )
        Nyval = Nyval.group(0)

        xvalsHanging = getData( textfoldername + "mixed_xvalues_Nx=" + Nxval + "Ny=" + Nyval + ".txt" )
        yvalsHanging = getData( textfoldername + "mixed_yvalues_Nx=" + Nxval + "Ny=" + Nyval + ".txt" )
        hangingErrvals = getData( textfoldername + "mixed_errorvalues_Nx=" + Nxval + "Ny=" + Nyval + ".txt" )

        curmaxval = np.max(hangingErrvals)
        curminval = np.min(hangingErrvals)

        hangingMaxErrorList.append( curmaxval )
        hangingl2ErrorList.append( np.sum( np.array(hangingErrvals)**2 ) ) 

        curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"

        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
        # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
        # print( np.max(errvals) )

        minVal, maxVal = minMaxRangeVals[ index ]

        levels = np.linspace( minVal, maxVal, 40 )
        plt.figure()
        # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = levelsmin[-1], vmax = levelsmax[-1] )
        plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = minVal, vmax = maxVal, cmap = cm )
        plt.colorbar()
        # plt.scatter( xvalsHanging, yvalsHanging )
        # plt.show()
        plt.savefig( curPlotFolderName + "hangingErrorcontour_Nx=" + Nxval + "Ny=" + Nyval + ".png" )
        plt.close()

        plt.figure()
        # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
        levels = np.linspace( curminval*0.4+ curmaxval*0.6, curmaxval, 40)
        # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
        # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
        plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors = 'r')
        plt.colorbar()
        plt.savefig( curPlotFolderName + "large_hangingErrorcontour_Nx=" + Nxval + "Ny=" + Nyval + ".png" )
        plt.close()
        # plt.figure()
        # plt.tricontourf( xvals, yvals, uvals, levels = 20 )
        # plt.colorbar()
        # plt.scatter( xvals, yvals )
        # plt.show()
        # plt.savefig( plotfoldername + "ucontour_index=" + str(index) + ".png" )

        # plt.figure()
        # plt.tricontourf( xvals, yvals, uexactvals, levels = 20 )
        # plt.colorbar()
        # plt.scatter( xvals, yvals )
        # plt.show()
        # plt.savefig( plotfoldername + "uexactcontour_index=" + str(index) + ".png" )
    
    for index, meshval in enumerate(sortedMeshValsRegular):

        regularNval = re.search( "[0-9]+"  , meshval  )
        regularNval = regularNval.group(0)

        xvals = getData( textfoldername + "regular_xvalues_N=" + regularNval + ".txt" )
        yvals = getData( textfoldername + "regular_yvalues_N=" + regularNval + ".txt" )
        regularErrvals = getData( textfoldername + "regular_errorvalues_N=" + regularNval + ".txt" )

        curPlotFolderName = simPlotFolderName + "Plot" + str(index) + "/"
        
        if not os.path.exists( curPlotFolderName ):
            os.mkdir( curPlotFolderName )

        # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

        curminval = np.min(regularErrvals)
        curmaxval = np.max(regularErrvals)

        regularMaxErrorList.append( curmaxval )
        regularl2ErrorList.append( np.sum( np.array(regularErrvals)**2 ) ) 

        minVal, maxVal = minMaxRangeVals[ index ]

        levels = np.linspace(minVal, maxVal, 40)
        plt.figure()
        # plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, vmin = levelsmin[indexval], vmax = levelsmax[indexval] )
        plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, vmin = minVal, vmax = maxVal, cmap = cm )
        plt.colorbar()
        # plt.scatter( xvals, yvals )
        # plt.show()
        plt.savefig( curPlotFolderName + "regularErrorcontour_N=" + regularNval + ".png" )
        plt.close()
        
        plt.figure()
        levels = np.linspace(curminval*0.4 + curmaxval*0.6, curmaxval, 40)
        # plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors='r')
        plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors = 'r')
        plt.colorbar()
        plt.savefig( curPlotFolderName + "large_regularErrorcontour_N=" + regularNval + ".png" )
        plt.close()

    plt.figure()
    plt.loglog( dxvals, hangingMaxErrorList, "-o", label = "Boundary Refined $L^{\infty}$ Error" )
    plt.loglog( dxvals, regularMaxErrorList, "-o", label = "Regular $L^{\infty}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( simPlotFolderName + "maxError" + ".png" )

    plt.figure()
    plt.loglog( dxvals, hangingl2ErrorList, "-o", label = "Boundary Refined $L^{2}$ Error" )
    plt.loglog( dxvals, regularl2ErrorList, "-o", label = "Regular $L^{2}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( simPlotFolderName + "l2Error" + ".png" )

    plt.show()
    return 1


def showplotTriangle():

    simPlotFolderName = plotfoldername + "Plot8/"

    os.mkdir( simPlotFolderName )

    dxfilename = "/home/gaurav/gmshAutoScripts/build/outfileregular.txt"

    dxvals = []

    with open( dxfilename ) as dxfilehandle:

        dxvals = dxfilehandle.readlines()

    dxvals = np.array([ float( dxval[ :-1 ] ) for dxval in dxvals ])

    import matplotlib.ticker as ticker

    triangleMaxErrorVals = dict()
    trianglel2ErrorVals = dict()

    regularMaxErrorVals = dict()
    regularl2ErrorVals = dict()

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]
    levelsmintriangle = dict()
    levelsmaxtriangle = dict()

    levelsminregular = dict()
    levelsmaxregular = dict()

    cdict = {
        'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
        'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
        'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }
 
    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
 
    indexval = 0
    for index, meshval in enumerate(meshvals):

        if re.search( "triangle", meshval  ):
            # print(meshval)
            triangleNval = re.search( "[0-9]+"  , meshval  )
            triangleNval = triangleNval.group(0)
            triangleErrvals = getData( textfoldername + "triangle_errorvalues_N=" + triangleNval + ".txt" )
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

            curminval = np.min(triangleErrvals)
            curmaxval = np.max(triangleErrvals)
            N = int(triangleNval)
            levelsmintriangle[N] = curminval, meshval 
            levelsmaxtriangle[N] = curmaxval, meshval 
            indexval += 1


    indexval = 0
    for index, meshval in enumerate(meshvals):

        if re.search( "regular", meshval  ):
            # print(meshval)
            regularNval = re.search( "[0-9]+"  , meshval  )
            regularNval = regularNval.group(0)

            regularErrvals = getData( textfoldername + "regular_errorvalues_N=" + regularNval + ".txt" )
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

            curminval = np.min(regularErrvals)
            curmaxval = np.max(regularErrvals)
            N = int(regularNval)
            levelsminregular[N] = curminval, meshval 
            levelsmaxregular[N] = curmaxval, meshval 
            indexval += 1

    sortedkeystriangle = sorted( levelsmintriangle.keys() )
    sortedkeysregular = sorted( levelsminregular.keys() )    

    levelsmin = dict()
    levelsmax = dict()

    for idx in range( len(sortedkeystriangle) ):

        minval = min( levelsminregular[ sortedkeysregular[idx] ][0], levelsmintriangle[ sortedkeystriangle[idx] ][0] )
        maxval = max( levelsmaxregular[ sortedkeysregular[idx] ][0], levelsmaxtriangle[ sortedkeystriangle[idx] ][0] )

        levelsmin[ levelsminregular[ sortedkeysregular[idx] ][1] ] = minval
        levelsmax[ levelsmaxregular[ sortedkeysregular[idx] ][1] ] = maxval

        levelsmin[ levelsmintriangle[ sortedkeystriangle[idx] ][1] ] = minval
        levelsmax[ levelsmaxtriangle[ sortedkeystriangle[idx] ][1] ] = maxval


    indexval = 0
    for index, meshval in enumerate(meshvals):

        if re.search( "triangle", meshval  ):

            triangleNval = re.search( "[0-9]+"  , meshval  )
            triangleNval = triangleNval.group(0)

            xvals = getData( textfoldername + "triangle_xvalues_N=" + triangleNval + ".txt" )
            yvals = getData( textfoldername + "triangle_yvalues_N=" + triangleNval + ".txt" )
            triangleErrvals = getData( textfoldername + "triangle_errorvalues_N=" + triangleNval + ".txt" )
            keyval = int(triangleNval)

            folderIndex = sortedkeystriangle.index( keyval )
            curPlotFolderName = simPlotFolderName + "Plot" + str(folderIndex) + "/"
            
            if not os.path.exists( curPlotFolderName ):
                os.mkdir( curPlotFolderName )

            triangleMaxErrorVals[ keyval ] = np.max(triangleErrvals) 
            trianglel2ErrorVals[ keyval ] = np.sum( np.array(triangleErrvals)**2 )
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

            curminval = np.min(triangleErrvals)
            curmaxval = np.max(triangleErrvals)
            # print(curminval, curmaxval)
            # print(levelsmin[indexval], levelsmax[indexval])
            levels = np.linspace(levelsmin[meshval], levelsmax[meshval], 40)
            plt.figure()
            # plt.tricontourf( xvals, yvals, triangleErrvals, levels = levels, vmin = levelsmin[indexval], vmax = levelsmax[indexval] )
            plt.tricontourf( xvals, yvals, triangleErrvals, levels = levels, vmin = levelsmin[meshval], vmax = levelsmax[meshval], cmap = cm )
            plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            plt.savefig( curPlotFolderName + "triangleErrorcontour_N=" + triangleNval + ".png" )
            plt.close()
            
            plt.figure()
            levels = np.linspace(curminval*0.4 + curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvals, yvals, triangleErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, triangleErrvals, levels = levels, colors = 'r')
            plt.colorbar()
            plt.savefig( curPlotFolderName + "large_triangleErrorcontour_N=" + triangleNval + ".png" )
            plt.close()
            indexval += 1

    indexval = 0
    for index, meshval in enumerate(meshvals):

        if re.search( "regular", meshval  ):
            regularNval = re.search( "[0-9]+"  , meshval  )
            regularNval = regularNval.group(0)

            xvals = getData( textfoldername + "regular_xvalues_N=" + regularNval + ".txt" )
            yvals = getData( textfoldername + "regular_yvalues_N=" + regularNval + ".txt" )
            regularErrvals = getData( textfoldername + "regular_errorvalues_N=" + regularNval + ".txt" )
            
            keyval = int(regularNval)
            folderIndex = sortedkeysregular.index( keyval )
            curPlotFolderName = simPlotFolderName + "Plot" + str(folderIndex) + "/"
            
            if not os.path.exists( curPlotFolderName ):
                os.mkdir( curPlotFolderName )

            regularMaxErrorVals[ keyval ] = np.max(regularErrvals) 
            regularl2ErrorVals[ keyval ] = np.sum( np.array(regularErrvals)**2 )
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

            curminval = np.min(regularErrvals)
            curmaxval = np.max(regularErrvals)
            # print(curminval, curmaxval)
            # print(levelsmin[indexval], levelsmax[indexval])
            levels = np.linspace(levelsmin[meshval], levelsmax[meshval], 40)
            plt.figure()
            # plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, vmin = levelsmin[indexval], vmax = levelsmax[indexval] )
            plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, vmin = levelsmin[meshval], vmax = levelsmax[meshval], cmap = cm )
            plt.colorbar()
            # plt.scatter( xvals, yvals )
            # plt.show()
            plt.savefig( curPlotFolderName + "regularErrorcontour_N=" + regularNval + ".png" )
            plt.close()
            
            plt.figure()
            levels = np.linspace(curminval*0.4 + curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors = 'r')
            plt.colorbar()
            plt.savefig( curPlotFolderName + "large_regularErrorcontour_N=" + regularNval + ".png" )
            plt.close()
            indexval += 1

    triangleMaxErrorList = []
    regularMaxErrorList = []
    trianglel2ErrorList = []
    regularl2ErrorList = []

    for idx, keyval in enumerate( sortedkeystriangle ):

        triangleMaxErrorList.append( triangleMaxErrorVals[ keyval ] )
        trianglel2ErrorList.append( trianglel2ErrorVals[ keyval ] )

    for idx, keyval in enumerate( sortedkeysregular ):

        regularMaxErrorList.append( regularMaxErrorVals[ keyval ] )
        regularl2ErrorList.append(  regularl2ErrorVals[ keyval ] )

    plt.figure()
    plt.loglog( dxvals, triangleMaxErrorList, "-o", label = "Triangle $L^{\infty}$ Error" )
    plt.loglog( dxvals, regularMaxErrorList, "-o", label = "Regular $L^{\infty}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( simPlotFolderName + "maxError" + ".png" )

    plt.figure()
    plt.loglog( dxvals, trianglel2ErrorList, "-o", label = "Triangle $L^{2}$ Error" )
    plt.loglog( dxvals, regularl2ErrorList, "-o", label = "Regular $L^{2}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( simPlotFolderName + "l2Error" + ".png" )

    plt.show()
    return 1

if __name__ == "__main__":

    simPlotFolderName = gmshImageFolderName + "Plot11_2pi_Hanging/"
    # runSim( simPlotFolderName )
    showParaviewPlot( simPlotFolderName )
    # showplot( simPlotFolderName )
    # showplotTriangle()

