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

Ndashvals = [ 9 ]

foldername = "/home/gaurav/gmshAutoScripts/build/"

lvl = 0

for Ndash in Ndashvals:

    Nval = int( ( Ndash + 5 )/2 )
    lc = 1/( Ndash + 2*Nval + 1 )
    print(lc)

    N = int( 1/lc + 1 )
    n3DElements = N

    linesval = []
    ptsval = []

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


