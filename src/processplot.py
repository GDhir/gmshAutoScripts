import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import re
from os import listdir
from os.path import isfile, join
import subprocess

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

    juliapath = "/home/gaurav/julia-1.8.5/bin/julia"
    julialstcmd = [juliapath, exefilename]
    subprocess.run( julialstcmd )

def removeFiles( dirval ):

    meshvals = [join(dirval, f) for f in listdir(dirval) if isfile(join(dirval, f))]

    for meshfilename in meshvals:
        cmdlist = [ "rm", meshfilename ]
        subprocess.run( cmdlist )

def runSim():

    gmshfileargs = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    removeFiles( gmshfileargs )

    gmshfilecmd = "hangingMeshv7"
    buildMesh( gmshfilecmd, gmshfileargs )

    gmshfilecmd = "regularMeshv2"
    buildMesh( gmshfilecmd, gmshfileargs )

    exefilename = "example-mixed-element-2d.jl"
    runJulia( exefilename )

    showplot()
    # removeFiles( gmshfileargs )

def showplot():

    rootfoldername = "/media/gaurav/easystore/Finch/MixedElement/"

    textfoldername = rootfoldername + "TextFiles/"

    plotfoldername = rootfoldername + "PlotFiles/SimPlots/"

    # indexes = [ 1, 2, 3, 4, 5 ]
    # Nyvals = np.array( [ 9, 17, 33, 65, 129 ] )
    Nvals = np.array( [ 11, 17, 21, 27, 31 ] )
    dxvals = 2/( 4*Nvals - 6 )
    # indexes = [1, 2]
    import matplotlib.ticker as ticker

    hangingMaxErrorVals = []
    hangingl2ErrorVals = []

    regularMaxErrorVals = []
    regularl2ErrorVals = []

    meshpath = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/"
    meshvals = [f for f in listdir(meshpath) if isfile(join(meshpath, f))]
    levelsminhanging = dict()
    levelsmaxhanging = dict()

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

        if not re.search( "regular", meshval  ):
            Nxval = re.search( "Nx=", meshval )
            offset = Nxval.start() + 3
            Nxval = re.search( "[0-9]+", meshval[offset:] )
            Nxval = Nxval.group(0)

            Nyval = re.search( "Ny=", meshval )
            offset = Nyval.start() + 3
            Nyval = re.search( "[0-9]+", meshval[offset:] )
            Nyval = Nyval.group(0)

            hangingErrvals = getData( textfoldername + "mixed_errorvalues_Nx=" + Nxval + "Ny=" + Nyval + ".txt" )
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )

            curminval = np.min(hangingErrvals)
            curmaxval = np.max(hangingErrvals)

            Nx = int(Nxval)
            Ny = int(Nyval)
            levelsminhanging[Nx*Ny] = curminval, meshval 
            levelsmaxhanging[Nx*Ny] = curmaxval, meshval 
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

    sortedkeyshanging = sorted( levelsminhanging.keys() )
    sortedkeysregular = sorted( levelsminregular.keys() )    

    levelsmin = dict()
    levelsmax = dict()

    for idx in range( len(sortedkeyshanging) ):

        minval = min( levelsminregular[ sortedkeysregular[idx] ][0], levelsminhanging[ sortedkeyshanging[idx] ][0] )
        maxval = max( levelsmaxregular[ sortedkeysregular[idx] ][0], levelsmaxhanging[ sortedkeyshanging[idx] ][0] )

        levelsmin[ levelsminregular[ sortedkeysregular[idx] ][1] ] = minval
        levelsmax[ levelsmaxregular[ sortedkeysregular[idx] ][1] ] = maxval

        levelsmin[ levelsminhanging[ sortedkeyshanging[idx] ][1] ] = minval
        levelsmax[ levelsmaxhanging[ sortedkeyshanging[idx] ][1] ] = maxval


    indexval = 0
    for index, meshval in enumerate(meshvals):

        if not re.search( "regular", meshval  ):

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
            hangingMaxErrorVals.append( np.max(hangingErrvals) )
            hangingl2ErrorVals.append( np.sum( np.array(hangingErrvals)**2 ) )
            curminval = np.min( hangingErrvals )
            curmaxval = np.max(hangingErrvals)
            # uvals = getData( textfoldername + "uvalues_index=" + str(index) + ".txt" )
            # uexactvals = getData( textfoldername + "uexactvalues_index=" + str(index) + ".txt" )        
            # print( np.max(errvals) )

            levels = np.linspace(levelsmin[meshval], levelsmax[meshval], 40)
            plt.figure()
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = levelsmin[-1], vmax = levelsmax[-1] )
            plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, vmin = levelsmin[meshval], vmax = levelsmax[meshval], cmap = cm )
            plt.colorbar()
            # plt.scatter( xvalsHanging, yvalsHanging )
            # plt.show()
            plt.savefig( plotfoldername + "hangingErrorcontour_Nx=" + Nxval + "Ny=" + Nyval + ".png" )
            plt.close()

            plt.figure()
            # levels = np.linspace( levelsmin[indexval]*0.4 + levelsmax[indexval]*0.6, levelsmax[indexval], 40)
            levels = np.linspace( curminval*0.4 + curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, vmin = levelsmin[-1]*0.2 + levelsmax[-1]*0.8, vmax = levelsmax[-1], colors='r')
            # plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors='r')
            plt.tricontourf( xvalsHanging, yvalsHanging, hangingErrvals, levels = levels, colors = 'r')
            plt.colorbar()
            plt.savefig( plotfoldername + "large_hangingErrorcontour_Nx=" + Nxval + "Ny=" + Nyval + ".png" )
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
    

    indexval = 0
    for index, meshval in enumerate(meshvals):

        if re.search( "regular", meshval  ):
            regularNval = re.search( "[0-9]+"  , meshval  )
            regularNval = regularNval.group(0)

            xvals = getData( textfoldername + "regular_xvalues_N=" + regularNval + ".txt" )
            yvals = getData( textfoldername + "regular_yvalues_N=" + regularNval + ".txt" )
            regularErrvals = getData( textfoldername + "regular_errorvalues_N=" + regularNval + ".txt" )
            regularMaxErrorVals.append( np.max(regularErrvals) )
            regularl2ErrorVals.append( np.sum( np.array(regularErrvals)**2 ) )
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
            plt.savefig( plotfoldername + "regularErrorcontour_N=" + regularNval + ".png" )
            plt.close()
            
            plt.figure()
            levels = np.linspace(curminval*0.4 + curmaxval*0.6, curmaxval, 40)
            # plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors='r')
            plt.tricontourf( xvals, yvals, regularErrvals, levels = levels, colors = 'r')
            plt.colorbar()
            plt.savefig( plotfoldername + "large_regularErrorcontour_N=" + regularNval + ".png" )
            plt.close()
            indexval += 1


    plt.figure()
    plt.loglog( dxvals, hangingMaxErrorVals, "-o", label = "Boundary Refined $L^{\infty}$ Error" )
    plt.loglog( dxvals, regularMaxErrorVals, "-o", label = "Regular $L^{\infty}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( plotfoldername + "maxError" + ".png" )

    plt.figure()
    plt.loglog( dxvals, hangingl2ErrorVals, "-o", label = "Boundary Refined $L^{2}$ Error" )
    plt.loglog( dxvals, regularl2ErrorVals, "-o", label = "Regular $L^{2}$ Error" )
    plt.loglog( dxvals, dxvals**2, "-x", label = "$h^2$" )
    plt.legend()
    plt.savefig( plotfoldername + "l2Error" + ".png" )

    plt.show()
    return 1

if __name__ == "__main__":

    # runSim()
    showplot()
    # a = 3
    # print(a)
