import sys

sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh

def createBox( xoffset, yoffset, pts, lines, lc, Nx, Ny ):

    side = 1
    i = 0

    while i < Nx:

        x = xoffset + lc * i
        y = yoffset + 0.0

        pts.append( gmsh.model.geo.addPoint( x, y, 0, lc ) )
        i += 1

    side = 2

    j = 1

    while i < Nx + Ny - 1:

        x = xoffset + (Nx - 1)*lc
        y = yoffset + j * lc

        pts.append( gmsh.model.geo.addPoint( x, y, 0, lc ) )

        j += 1
        i += 1

    side = 3
    j = 1

    while i < (Nx - 1) * 2 + Ny:

        x = xoffset + (Nx - 1) * lc - j * lc
        y = yoffset + (Ny - 1) * lc

        pts.append( gmsh.model.geo.addPoint( x, y, 0, lc ) )

        j = j + 1
        i = i + 1
    
    side = 4
    j = 1

    while i < (Nx - 1) * 2 + (Ny - 1) * 2:
    
        x = xoffset + 0
        y = yoffset + (Ny - 1) * lc - j * lc

        pts.append( gmsh.model.geo.addPoint( x, y, 0, lc ) )
        j = j + 1
        i = i + 1

    for i in range( len(pts) - 1 ):
        lines.append( gmsh.model.geo.addLine(pts[i], pts[i + 1]) )

    lines.append( gmsh.model.geo.addLine( pts[-1], pts[0] ) )
    
def setTransfiniteCurves( linesVec, N ):

    for lines in linesVec:
        [ gmsh.model.geo.mesh.setTransfiniteCurve( line, N ) for line in lines ]
        
def setTransfiniteSurfaces( planeIds ):

    [ gmsh.model.geo.mesh.setTransfiniteSurface(id) \
      for id in planeIds ]


def recombineSurfaces( planeIds ):

    [ gmsh.model.mesh.setRecombine( 2, id ) for id in planeIds ]