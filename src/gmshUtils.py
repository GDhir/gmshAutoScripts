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

class Zone:

    def __init__( self, pointSet, xoffset, yoffset ):

        self.pts = []
        self.Nx = pointSet[0]
        self.Ny = pointSet[1]
        self.h = pointSet[2]
        
        self.addPointSet( xoffset, yoffset )

        self.linesx = []
        self.linesy = []
        self.addLinesInZone( pointSet )

        self.curveloops = []
        self.surfaces = []
        self.addSurfaces( pointSet )

    def addPointSet( self, xoffset, yoffset ):

        yval = yoffset
        xval = xoffset

        for yidx in range(self.Ny):

            yval = yoffset + yidx * self.h
            for xidx in range( self.Nx ):

                xval = xoffset + xidx * self.h
                self.pts.append( gmsh.model.geo.addPoint( xval, yval, 0, self.h ) )

    def addLinesInZone( self, pointSet ):

        for idy in range( self.Ny ):
            for idx in range( self.Nx - 1 ):

                indexval = idy * ( self.Nx ) + idx

                self.linesx.append( gmsh.model.geo.addLine( self.pts[ indexval ], self.pts[ indexval + 1 ] ) )

        for idy in range( self.Ny - 1 ):
            for idx in range( self.Nx ):

                indexval = idy * self.Nx + idx
                self.linesy.append( gmsh.model.geo.addLine( self.pts[ indexval ], self.pts[ indexval + self.Nx ] ) )

    def addSurfaces( self ):

        Nx1 = self.Nx - 1
        Nx2 = self.Nx

        for idy in range( self.Ny - 1 ):
            for idx in range( self.Nx - 1 ):

                index1 = idy * Nx1 + idx
                index2 = idy * Nx2 + idx

                downEdge = self.linesx[ index1 ]
                upEdge = self.linesx[ index1 + Nx1 ]
                leftEdge = self.linesy[ index2 ]
                rightEdge = self.linesy[ index2 + 1 ]

                linesVal = [ downEdge, rightEdge, -upEdge, -leftEdge ]

                cl = gmsh.model.geo.addCurveLoop( linesVal )
                self.curveloops.append( cl )
                
                pl = gmsh.model.geo.addPlaneSurface( [ cl ])
                self.surfaces.append( pl )    

    def getBoundaryPoints( self, boundaryVal = "left" ):

        if boundaryVal == "left":
            
            indexStart = 0
            indexEnd = self.Ny * self.Nx

            bdryPts = self.pts[ indexStart : indexEnd : self.Nx ]        

        elif boundaryVal == "right":

            indexStart = self.Nx - 1
            indexEnd = ( self.Ny + 1 ) * self.Nx

            bdryPts = self.pts[ indexStart : indexEnd : self.Nx ] 
            

        elif boundaryVal == "top":

            indexStart = ( self.Ny - 1 ) * self.Nx
            bdryPts = self.pts[ indexStart : indexStart + self.Nx ]

        elif boundaryVal == "bottom":

            bdryPts = self.pts[ 0 : self.Nx ]   

        return bdryPts