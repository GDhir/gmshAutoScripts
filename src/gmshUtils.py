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

def recombineZones( zones ):

    for zone in zones:

        surfaceVals = zone.surfaces
        recombineSurfaces( surfaceVals )

def recombineLayers( layers ):

    for layer in layers:

        surfacesX = layer.surfacesz_xperp
        surfacesY = layer.surfacesz_yperp
        
        recombineSurfaces( surfacesX )
        recombineSurfaces( surfacesY )

def recombine3DZone( zone ):

    recombineZones( zone.zones2D )
    recombineLayers( zone.layers )

class Layer3D:

    def __init__( self, zone1, zone2, Nx, Ny, transfinite = False ):

        self.linesz = []
        self.Nx = Nx
        self.Ny = Ny
        self.surfacesz_yperp = []
        self.surfacesz_xperp = []
        self.curveloopsz = []
        self.surfaceloops = []
        self.volumes = []

        self.addLinesBetweenLayers( zone1, zone2 )
        self.addSurfacesBetweenLayers( zone1, zone2 )
        self.addVolumesBetweenLayers( zone1, zone2 )

        if transfinite:
            linesVec = [ self.linesz ]
            setTransfiniteCurves( linesVec, 2 )
            
            setTransfiniteSurfaces( self.surfacesz_xperp )
            setTransfiniteSurfaces( self.surfacesz_yperp )

    def addLinesBetweenLayers( self, zone1, zone2 ):

        for idy in range( self.Ny ):
            for idx in range( self.Nx ):

                indexval = idy * ( self.Nx ) + idx
                self.linesz.append( gmsh.model.geo.addLine( zone1.pts[ indexval ], zone2.pts[ indexval ] ) )

    def addSurfacesBetweenLayers( self, zone1, zone2 ):

        for idy in range( self.Ny ):
            for idx in range( self.Nx - 1 ):

                indexval1 = idy * ( self.Nx - 1 ) + idx
                indexval2 = idy * ( self.Nx ) + idx

                downEdge = zone1.linesx[ indexval1 ]
                upEdge = zone2.linesx[ indexval1 ]
                zEdge1 = self.linesz[ indexval2 ]
                zEdge2 = self.linesz[ indexval2 + 1 ]

                linesVal = [ -downEdge, zEdge1, upEdge, -zEdge2 ]

                cl = gmsh.model.geo.addCurveLoop( linesVal )
                self.curveloopsz.append( cl )
                
                pl = gmsh.model.geo.addPlaneSurface( [ cl ])
                self.surfacesz_yperp.append( pl )

        for idy in range( self.Ny - 1 ):
            for idx in range( self.Nx ):

                indexval = idy * ( self.Nx ) + idx

                downEdge = zone1.linesy[ indexval ]
                upEdge = zone2.linesy[ indexval ]
                zEdge1 = self.linesz[ indexval ]
                zEdge2 = self.linesz[ indexval + self.Nx ]

                linesVal = [ -downEdge, zEdge1, upEdge, -zEdge2 ]

                print(linesVal)
                cl = gmsh.model.geo.addCurveLoop( linesVal )
                self.curveloopsz.append( cl )
                
                pl = gmsh.model.geo.addPlaneSurface( [ cl ])
                self.surfacesz_xperp.append( pl )   

    def addVolumesBetweenLayers( self, zone1, zone2 ):

        Nx1 = self.Nx - 1
        Nx2 = self.Nx

        for idy in range( self.Ny - 1 ):
            for idx in range( self.Nx - 1 ):

                indexval1 = idy * Nx1 + idx
                indexval2 = idy * Nx2 + idx

                surface_zperp1 = zone1.surfaces[ indexval1 ]
                surface_zperp2 = zone2.surfaces[ indexval1 ]
                surface_yperp1 = self.surfacesz_yperp[ indexval1 ]
                surface_yperp2 = self.surfacesz_yperp[ indexval1 + Nx1 ]
                surface_xperp1 = self.surfacesz_xperp[ indexval2 ]
                surface_xperp2 = self.surfacesz_xperp[ indexval2 + 1 ]

                sl = gmsh.model.geo.addSurfaceLoop( [ surface_yperp1, surface_yperp2, surface_xperp1, surface_xperp2, \
                                                     surface_zperp1, surface_zperp2 ] )
                
                vl = gmsh.model.geo.addVolume( [sl] )
                self.volumes.append( vl )

class Zone3D:

    def __init__( self, zones, Nx, Ny, transfinite = False ):

        nlayers = len( zones ) - 1

        self.layers = []
        self.zones2D = zones

        for idx in range( nlayers ):
            self.layers.append( Layer3D( zones[idx], zones[idx + 1], Nx, Ny, transfinite = transfinite ) )

class Zone2D:

    def __init__( self, pointSet, xoffset, yoffset, zoffset = 0, transfinite = False ):

        self.pts = []
        self.Nx = pointSet[0]
        self.Ny = pointSet[1]
        self.h = pointSet[2]

        self.xoffset = xoffset
        self.yoffset = yoffset
        self.zoffset = zoffset
        
        self.addPointSet( xoffset, yoffset, zoffset )

        self.linesx = []
        self.linesy = []
        self.addLinesInZone( )

        self.curveloops = []
        self.surfaces = []
        self.addSurfaces( )

        if transfinite:
            linesVec = [ self.linesx, self.linesy ]
            setTransfiniteCurves( linesVec, 2 )
            setTransfiniteSurfaces( self.surfaces )

    def addPointSet( self, xoffset, yoffset, zval = 0 ):

        yval = yoffset
        xval = xoffset

        for yidx in range(self.Ny):

            yval = yoffset + yidx * self.h
            for xidx in range( self.Nx ):

                xval = xoffset + xidx * self.h
                self.pts.append( gmsh.model.geo.addPoint( xval, yval, zval, self.h ) )

    def addLinesInZone( self, startIdx = 0 ):

        for idy in range( self.Ny ):
            for idx in range( self.Nx - 1 ):

                indexval = startIdx + idy * ( self.Nx ) + idx

                self.linesx.append( gmsh.model.geo.addLine( self.pts[ indexval ], self.pts[ indexval + 1 ] ) )

        for idy in range( self.Ny - 1 ):
            for idx in range( self.Nx ):

                indexval = startIdx + idy * self.Nx + idx
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

    def addLayer( self, xoffset, yoffset, zval ):

        idxOffset = len( self.pts )

        self.addPointSet( xoffset, yoffset, zval )
        self.addLinesInZone( idxOffset )

        return