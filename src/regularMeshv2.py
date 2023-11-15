# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 2
#
#  Transformations, extruded geometries, volumes
#
# ------------------------------------------------------------------------------
import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh
import math
import gmshUtils

# If sys.argv is passed to gmsh.initialize(), Gmsh will parse the command line
# in the same way as the standalone Gmsh app:

def singleZoneMesh():

    Ndashvals = [ 9 ]

    foldername = "/home/gaurav/gmshAutoScripts/build/"

    lvl = 0

    for Ndash in Ndashvals:

        Nval = int( ( Ndash + 5 )/2 )
        lc = 1/( Ndash + 2*Nval + 1 )
        print(lc)

        N = int( 1/lc + 1 )
        n3DElements = N

        gmsh.initialize(sys.argv)
        gmsh.model.add("t2")

        pointSet = [N, N, lc]
        xoffset = 0
        yoffset = 0

        zone1 = gmshUtils.Zone( pointSet, xoffset, yoffset )
        gmsh.model.geo.synchronize()
        linesVec = [ zone1.linesx, zone1.linesy ]

        gmshUtils.setTransfiniteCurves( linesVec, 2 )
        gmshUtils.setTransfiniteSurfaces( zone1.surfaces )

        surfacePairs = [ ( 2, val ) for val in zone1.surfaces ]
        ov2 = gmsh.model.geo.extrude( surfacePairs, 0, 0, 1, numElements = [ n3DElements ], recombine = True )

        gmsh.model.geo.synchronize()
        gmshUtils.recombineSurfaces( zone1.surfaces )
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

        gmsh.model.mesh.generate(2)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2)

        regMeshFileName = foldername + "regularMesh_lvl" + str(lvl) + ".msh"
        gmsh.write( regMeshFileName )

        lvl = lvl + 1

        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

        gmsh.finalize()

def doubleZoneMesh():

    Ndash = 9
    foldername = "/home/gaurav/gmshAutoScripts/build/"
    lvl = 0

    Nx1 = 5
    Nx2 = 7
    lc = 1/( Nx1 + 2*Nx2 - 1 )
    Ny1 = int( 1/lc + 1 )
    Ny2 = int( 1/2/lc + 1 )
    print(lc)

    gmsh.initialize(sys.argv)
    gmsh.model.add("t2")

    pointSet1 = [Nx1, Ny1, lc]
    xoffset1 = 0
    yoffset1 = 0
    zone1 = gmshUtils.Zone( pointSet1, xoffset1, yoffset1 )
    gmsh.model.geo.synchronize()

    pointSet2 = [Nx2, Ny2, 2*lc]
    xoffset2 = (Nx1 + 1)*lc
    yoffset2 = 0
    zone2 = gmshUtils.Zone( pointSet2, xoffset2, yoffset2 )
    gmsh.model.geo.synchronize()

    linesVec = [ zone1.linesx, zone1.linesy, zone2.linesx, zone2.linesy ]

    gmshUtils.setTransfiniteCurves( linesVec, 2 )
    gmshUtils.setTransfiniteSurfaces( zone1.surfaces )
    gmshUtils.setTransfiniteSurfaces( zone2.surfaces )

    surfacePairs1 = [ ( 2, val ) for val in zone1.surfaces ]
    ov2 = gmsh.model.geo.extrude( surfacePairs1, 0, 0, 1, numElements = [ Ny1 ], recombine = True )

    surfacePairs2 = [ ( 2, val ) for val in zone2.surfaces ]
    ov2 = gmsh.model.geo.extrude( surfacePairs2, 0, 0, 1, numElements = [ Ny2 ], recombine = True )

    gmsh.model.geo.synchronize()
    gmshUtils.recombineSurfaces( zone1.surfaces )
    gmshUtils.recombineSurfaces( zone2.surfaces )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = foldername + "regularMesh_lvl" + str(lvl) + ".msh"
    gmsh.write( regMeshFileName )

    lvl = lvl + 1

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

if __name__ == "__main__":

    doubleZoneMesh()