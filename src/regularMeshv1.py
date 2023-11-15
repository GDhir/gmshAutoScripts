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

Ndashvals = [ 9, 17, 33, 65, 129 ]

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

    ptsval.append( gmsh.model.geo.addPoint( 0, 0, 0, lc ) )
    ptsval.append( gmsh.model.geo.addPoint( 1, 0, 0, lc ) )
    ptsval.append( gmsh.model.geo.addPoint( 1, 1, 0, lc ) )
    ptsval.append( gmsh.model.geo.addPoint( 0, 1, 0, lc ) )

    linesval.append( gmsh.model.geo.addLine(ptsval[0], ptsval[1]) )
    linesval.append( gmsh.model.geo.addLine(ptsval[1], ptsval[2]) )
    linesval.append( gmsh.model.geo.addLine(ptsval[2], ptsval[3]) )
    linesval.append( gmsh.model.geo.addLine(ptsval[3], ptsval[0]) )
    
    cl = gmsh.model.geo.addCurveLoop( linesval )
    pl = gmsh.model.geo.addPlaneSurface( [ cl ])

    planeIds = [pl]
    linesVec = [ linesval ]
    gmshUtils.setTransfiniteCurves( linesVec, N )

    gmshUtils.setTransfiniteSurfaces( planeIds )

    gmsh.model.geo.synchronize()

    ov = gmsh.model.geo.copy( [ ( 2, pl ) ] )
    ov2 = gmsh.model.geo.extrude( [ov[0]], 0, 0, 1, numElements = [ n3DElements ], recombine = True )

    gmsh.model.geo.synchronize()
    gmshUtils.recombineSurfaces( planeIds )
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    regMeshFileName = foldername + "regularMesh_lvl" + str(lvl) + ".msh"
    gmsh.write( regMeshFileName )

    lvl = lvl + 1

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()


