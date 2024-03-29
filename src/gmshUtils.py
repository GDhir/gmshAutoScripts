import sys

sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh

def createCorners( xoffset, yoffset, pts, lines, lc, Nx, Ny ):

    pts.append( gmsh.model.geo.addPoint( xoffset, yoffset, 0, lc ) )
    pts.append( gmsh.model.geo.addPoint( xoffset + (Nx - 1)*lc, yoffset, 0, lc ) )
    pts.append( gmsh.model.geo.addPoint( xoffset + (Nx - 1)*lc, yoffset + (Ny - 1)*lc, 0, lc ) )
    pts.append( gmsh.model.geo.addPoint( xoffset, yoffset + (Ny - 1)*lc, 0, lc ) )

    for i in range( len(pts) - 1 ):
        lines.append( gmsh.model.geo.addLine(pts[i], pts[i + 1]) )

    lines.append( gmsh.model.geo.addLine( pts[-1], pts[0] ) )

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
    
def setTransfiniteCurves( linesVec, N, occ = False ):

    if not occ:
        for lines in linesVec:
            [ gmsh.model.geo.mesh.setTransfiniteCurve( line, N ) for line in lines ]
    else:
        for lines in linesVec:
            [ gmsh.model.mesh.setTransfiniteCurve( line, N ) for line in lines ]
        
def setTransfiniteSurfaces( planeIds, cornerPts = [], occ = False ):

    if not occ:
        if cornerPts:
            [ gmsh.model.geo.mesh.setTransfiniteSurface(id, "Left", cornerPts ) \
            for id in planeIds ]
        else:
            [ gmsh.model.geo.mesh.setTransfiniteSurface(id) \
            for id in planeIds ]
    else:
        if cornerPts:
            [ gmsh.model.mesh.setTransfiniteSurface(id, "Left", cornerPts ) \
            for id in planeIds ]
        else:
            [ gmsh.model.mesh.setTransfiniteSurface(id) \
            for id in planeIds ]

    
def setTransfiniteVolumes( planeIds ):

    [ gmsh.model.geo.mesh.setTransfiniteVolume(id) \
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
            setTransfiniteVolumes( self.volumes )

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

                # print(linesVal)
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

        self.nLayers = len( zones ) - 1
        self.Nx = Nx
        self.Ny = Ny

        self.layers = []
        self.zones2D = zones

        for idx in range( self.nLayers ):
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
    
class ZoneConnection:

    def __init__( self, zoneNeg, zonePos, dxn, transfinite = False ):

        self.surfacesNeg = []
        self.surfacesPos = []
        self.lineConnections = []
        self.volumes = []

        self.getZoneSurfaces( zoneNeg, zonePos, dxn )
        self.connectLines( zoneNeg, zonePos, dxn )

        self.surfaces_Type1 = []
        self.surfaces_Type2 = []

        self.connectSurfaces( zoneNeg, zonePos, dxn, False )
        self.addVolumes( zoneNeg, zonePos, dxn )

        if transfinite:
            
            linesVec = [ self.lineConnections ]
            setTransfiniteCurves( linesVec, 2 )
            # setTransfiniteSurfaces( self.surfaces_Type1 )
            # setTransfiniteSurfaces( self.surfaces_Type2 )

    def getZoneSurfaces( self, zoneNeg, zonePos, dxn = "x" ):

        if dxn == "x":

            for zIdx in range( zoneNeg.nLayers ):

                layerVal = zoneNeg.layers[ zIdx ]

                for yidx in range( zoneNeg.Ny - 1 ):

                    self.surfacesNeg.append( layerVal.surfacesz_xperp[ ( yidx + 1 ) * zoneNeg.Nx - 1 ] )                    

            for zIdx in range( zonePos.nLayers ):

                layerVal = zonePos.layers[ zIdx ]

                for yidx in range( zonePos.Ny - 1 ):

                    self.surfacesPos.append( layerVal.surfacesz_xperp[ ( yidx ) * zonePos.Nx ] )

        elif dxn == "y":
            
            for zIdx in range( zoneNeg.nLayers ):

                layerVal = zoneNeg.layers[ zIdx ]

                for xidx in range( zoneNeg.Nx - 1 ):

                    self.surfacesNeg.append( layerVal.surfacesz_yperp[ ( zoneNeg.Ny - 1 ) * ( zoneNeg.Nx - 1 ) + xidx ] )

            for zIdx in range( zonePos.nLayers ):

                layerVal = zonePos.layers[ zIdx ]

                for xidx in range( zonePos.Nx - 1 ):

                    self.surfacesPos.append( layerVal.surfacesz_yperp[ xidx ] )

        elif dxn == "z":

            zone2DNeg = zoneNeg.zones2D[ -1 ]
            zone2DPos = zonePos.zones2D[ 0 ]

            for yIdx in range( zoneNeg.Ny - 1 ):
                for xIdx in range( zoneNeg.Nx - 1 ):

                    indexval = yidx * ( zoneNeg.Nx - 1 ) + xidx

                    self.surfacesNeg.append( zone2DNeg.surfaces[ indexval ] )

            for yIdx in range( zonePos.Ny - 1 ):
                for xIdx in range( zonePos.Nx - 1 ):

                    indexval = yidx * ( zonePos.Nx - 1 ) + xidx
                    
                    self.surfacesPos.append( zone2DPos.surfaces[ indexval ] )


    def connectLines( self, zoneNeg, zonePos, dxn ):

        nZonesNeg = len( zoneNeg.zones2D )
        nZonesPos = len( zonePos.zones2D )

        nZonesToJoin = min( nZonesNeg, nZonesPos )
        
        if dxn == "x":

            nPtsToJoin = min( zoneNeg.Ny, zonePos.Ny )

            for zoneIdx in range( nZonesToJoin ):
                for yIdx in range( nPtsToJoin ):

                    if nZonesNeg < nZonesPos:

                        zone2DNeg = zoneNeg.zones2D[ zoneIdx ]
                        zone2DPos = zonePos.zones2D[ 2 * zoneIdx ]

                        negIdx = ( yIdx + 1 ) * zone2DNeg.Nx - 1
                        posIdx = yIdx * zone2DPos.Nx * 2

                    else:

                        zone2DNeg = zoneNeg.zones2D[ 2 * zoneIdx ]
                        zone2DPos = zonePos.zones2D[ zoneIdx ]

                        negIdx = ( 2 * yIdx + 1 ) * zone2DNeg.Nx - 1
                        posIdx = yIdx * zone2DPos.Nx

                    self.lineConnections.append( gmsh.model.geo.addLine( zone2DNeg.pts[ negIdx ], \
                                                                        zone2DPos.pts[ posIdx ] ) )

        elif dxn == "y":

            nPtsToJoin = min( zoneNeg.Nx, zonePos.Nx )

            for zoneIdx in range( nZonesToJoin ):
                for xIdx in range( nPtsToJoin ):

                    if nZonesNeg < nZonesPos:

                        zone2DNeg = zoneNeg.zones2D[ zoneIdx ]
                        zone2DPos = zonePos.zones2D[ 2 * zoneIdx ]

                        negIdx = ( zone2DNeg.Ny - 1 ) * zone2DNeg.Nx + xIdx
                        posIdx = 2 * xIdx

                    else:

                        zone2DNeg = zoneNeg.zones2D[ 2 * zoneIdx ]
                        zone2DPos = zonePos.zones2D[ zoneIdx ]

                        negIdx = ( zone2DNeg.Ny - 1 ) * zone2DNeg.Nx + 2 * xIdx
                        posIdx = xIdx

                    self.linesConnections.append( gmsh.model.geo.addLine( zone2DNeg.pts[ negIdx ], \
                                                                        zone2DPos.pts[ posIdx ] ) )

        elif dxn == "z":

            zone2DNeg = zoneNeg.zones2D[ -1 ]
            zone2DPos = zonePos.zones2D[ 0 ]

            Nxpts = min( zone2DNeg.Nx, zone2DPos.Nx )
            Nypts = min( zone2DNeg.Ny, zone2DPos.Ny )

            for yIdx in range( Nypts ):
                for xIdx in range( Nxpts ):

                    if zone2DNeg.Nx < zone2DPos.Nx:

                        negIdx = yIdx * Nxpts + xIdx
                        posIdx = 2 * ( yIdx * Nxpts + xIdx )
                        
                    else:

                        negIdx = 2 * ( yIdx * Nxpts + xIdx )
                        posIdx = yIdx * Nxpts + xIdx

                    self.linesConnections.append( gmsh.model.geo.addLine( zone2DNeg.pts[ negIdx ], \
                                                                        zone2DPos.pts[ posIdx ] ) )


    def connectSurfaces( self, zoneNeg, zonePos, dxn, transfinite = False ):

        nLayersNeg = zoneNeg.nLayers
        nLayersPos = zonePos.nLayers

        nZonesNeg = nLayersNeg + 1
        nZonesPos = nLayersPos + 1

        nLayersToJoin = min( nLayersNeg, nLayersPos )
        nZonesToJoin = nLayersToJoin + 1

        if dxn == "x":
            nSurfacesToJoin = min( zoneNeg.Ny, zonePos.Ny )

            for zIdx in range( nLayersToJoin ):
                for yIdx in range( nSurfacesToJoin ):

                    linesVal = []
                    downEdge = -self.lineConnections[ zIdx * nSurfacesToJoin + yIdx ]
                    linesVal.append( downEdge )

                    if nLayersNeg < nLayersPos:
                        linesVal.append( zoneNeg.layers[ zIdx ].linesz[ zoneNeg.Nx * ( yIdx + 1 ) - 1 ] )

                        cornerPts = [ zoneNeg.zones2D[ zIdx ].pts[ zoneNeg.Nx * ( yIdx + 1 ) - 1 ], \
                                    zonePos.zones2D[ 2 * ( zIdx ) ].pts[ zonePos.Nx * ( yIdx * 2 ) ], \
                                    zonePos.zones2D[ 2 * ( zIdx + 1 ) ].pts[ zonePos.Nx * ( yIdx * 2 ) ], \
                                    zoneNeg.zones2D[ zIdx + 1 ].pts[ zoneNeg.Nx * ( yIdx + 1 ) - 1 ] ]
                    else:
                        linesVal.append( zoneNeg.layers[ 2 * zIdx ].linesz[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ] )
                        linesVal.append( zoneNeg.layers[ 2 * zIdx + 1 ].linesz[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ] )

                        cornerPts = [ zoneNeg.zones2D[ 2 * zIdx ].pts[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ], \
                                    zonePos.zones2D[ zIdx ].pts[ zonePos.Nx * ( yIdx ) ], \
                                    zonePos.zones2D[ zIdx + 1 ].pts[ zonePos.Nx * ( yIdx ) ], \
                                    zoneNeg.zones2D[ 2 * ( zIdx + 1 ) ].pts[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ] ]

                    upEdge = self.lineConnections[ ( zIdx + 1 ) * nSurfacesToJoin + yIdx ]
                    linesVal.append( upEdge )

                    if nLayersPos < nLayersNeg:
                        linesVal.append( -zonePos.layers[ zIdx ].linesz[ zonePos.Nx * yIdx ] )
                    else:
                        linesVal.append( -zonePos.layers[ 2 * zIdx + 1 ].linesz[ zoneNeg.Nx * yIdx * 2 ] )
                        linesVal.append( -zonePos.layers[ 2 * zIdx ].linesz[ zoneNeg.Nx * yIdx * 2 ] )

                    cl = gmsh.model.geo.addCurveLoop( linesVal )
                    pl = gmsh.model.geo.addPlaneSurface( [ cl ] )
                    self.surfaces_Type1.append( pl ) 

                    if transfinite:
                        # print(pl)
                        # print( cornerPts )
                        setTransfiniteSurfaces( [pl], cornerPts )

            for zIdx in range( nZonesToJoin ):
                for yIdx in range( nSurfacesToJoin - 1 ):

                    linesVal = []
                    downEdge = self.lineConnections[ zIdx * nSurfacesToJoin + yIdx ]
                    linesVal.append( -downEdge )

                    if nZonesNeg < nZonesPos:
                        linesVal.append( zoneNeg.zones2D[ zIdx ].linesy[ zoneNeg.Nx * ( yIdx + 1 ) - 1 ] )

                        cornerPts = [ zoneNeg.zones2D[ zIdx ].pts[ zoneNeg.Nx * ( yIdx + 1 ) - 1 ], \
                                    zonePos.zones2D[ 2 * ( zIdx ) ].pts[ zonePos.Nx * ( yIdx * 2 ) ], \
                                    zonePos.zones2D[ 2 * ( zIdx ) ].pts[ zonePos.Nx * ( ( yIdx + 1 ) * 2 ) ], \
                                    zoneNeg.zones2D[ zIdx ].pts[ zoneNeg.Nx * ( yIdx + 2 ) - 1 ] ]
                 
                    else:
                        linesVal.append( zoneNeg.zones2D[ 2 * zIdx ].linesy[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ] )
                        linesVal.append( zoneNeg.zones2D[ 2 * zIdx ].linesy[ zoneNeg.Nx * ( yIdx * 2 + 2 ) - 1 ] )

                        cornerPts = [ zoneNeg.zones2D[ 2 * zIdx ].pts[ zoneNeg.Nx * ( yIdx * 2 + 1 ) - 1 ], \
                                    zonePos.zones2D[ zIdx ].pts[ zonePos.Nx * ( yIdx ) ], \
                                    zonePos.zones2D[ zIdx ].pts[ zonePos.Nx * ( yIdx + 1 ) ], \
                                    zoneNeg.zones2D[ 2 * ( zIdx ) ].pts[ zoneNeg.Nx * ( yIdx * 2 + 2 ) - 1 ] ]

                    upEdge = self.lineConnections[ zIdx * nSurfacesToJoin + yIdx + 1 ]
                    linesVal.append( upEdge )

                    if nZonesPos < nZonesNeg:
                        linesVal.append( -zonePos.zones2D[ zIdx ].linesy[ zonePos.Nx * yIdx ] )
                    else:
                        linesVal.append( -zonePos.zones2D[ 2 * zIdx ].linesy[ zoneNeg.Nx * ( yIdx * 2 + 1 ) ] )
                        linesVal.append( -zonePos.zones2D[ 2 * zIdx ].linesy[ zoneNeg.Nx * ( yIdx * 2 ) ] )

                    cl = gmsh.model.geo.addCurveLoop( linesVal )
                    pl = gmsh.model.geo.addPlaneSurface( [ cl ] )
                    self.surfaces_Type2.append( pl )

                    if transfinite:
                        # print(pl)
                        # print( cornerPts )
                        setTransfiniteSurfaces( [pl], cornerPts )

        elif dxn == "y":
            pass

        elif dxn == "z":
            pass


    def addVolumes( self, zoneNeg, zonePos, dxn ):

        if dxn == "x":

            Nz = min( len( zoneNeg.zones2D ), len( zonePos.zones2D ) )
            Ny = min( zoneNeg.Ny, zonePos.Ny )

            for zIdx in range( Nz - 1 ):
                for yIdx in range( Ny - 1 ):

                    surfaceVals = []
                    surfaceVals.append( self.surfaces_Type1[ zIdx * Ny + yIdx ] )
                    surfaceVals.append( self.surfaces_Type1[ zIdx * Ny + yIdx + 1 ] )

                    if zoneNeg.nLayers < zonePos.nLayers:

                        surfaceVals.append( self.surfacesNeg[ zIdx * ( Ny - 1 ) + yIdx ] )
                        surfaceVals.append( self.surfacesPos[ 2 * zIdx * 2 *( Ny - 1 ) + 2 * yIdx ] )
                        surfaceVals.append( self.surfacesPos[ 2 * zIdx * 2 *( Ny - 1 ) + 2 * yIdx + 1 ] )
                        surfaceVals.append( self.surfacesPos[ ( 2 * zIdx + 1 ) *  2 * ( Ny - 1 ) + 2 * yIdx ] )
                        surfaceVals.append( self.surfacesPos[ ( 2 * zIdx + 1 ) * 2 * ( Ny - 1 ) + 2 * yIdx + 1 ] )

                    else:

                        surfaceVals.append( self.surfacesNeg[ 2 * zIdx * 2 * ( Ny - 1 ) + 2 * yIdx ] )
                        surfaceVals.append( self.surfacesNeg[ 2 * zIdx * 2 * ( Ny - 1 ) + 2 * yIdx + 1 ] )
                        surfaceVals.append( self.surfacesNeg[ ( 2 * zIdx + 1 ) * 2 * ( Ny - 1 ) + 2 * yIdx ] )
                        surfaceVals.append( self.surfacesNeg[ ( 2 * zIdx + 1 ) * 2 * ( Ny - 1 ) + 2 * yIdx + 1 ] )
                        surfaceVals.append( self.surfacesPos[ zIdx * ( Ny - 1 ) + yIdx ] )

                    surfaceVals.append( self.surfaces_Type2[ zIdx * ( Ny - 1 ) + yIdx ] )
                    surfaceVals.append( self.surfaces_Type2[ ( zIdx + 1 ) * ( Ny - 1 ) + yIdx ] )

                    sl = gmsh.model.geo.addSurfaceLoop( surfaceVals )
                 
                    vl = gmsh.model.geo.addVolume( [sl] )
                    self.volumes.append( vl )

        else:

            pass

def checkPtOnBoundary( pointCoords ):

    boundaryCoordVals = [0, 1]

    for pointCoord in pointCoords:

        for bcval in boundaryCoordVals:
            if abs( pointCoord - bcval ) < 1e-6:
                return True
        
    return False

def checkElementOnBoundary( vertexIndices, vertexCoords ):

    for vertexIndex in vertexIndices:

        if not checkPtOnBoundary( vertexCoords[ vertexIndex ] ):
            return False
        
    return True

def parseGMSHFile( gmshFileName ):

    etype2nv = dict( [("1", "2"), ("2", "3"), ("3", "4"), ("4", "4"), ("5", "8"), ("6", "6"),\
                       ("7", "5"), ("15", "1") ] )
    
    etype2dim = dict( [("1", 1), ("2", 2), ("3", 2), ("4", 3), ("5", 3), ("6", 3),\
                       ("7", 3), ("15", 0) ] )
    lineIndicesToDelete = set()

    vertexCoords = dict()    

    with open( gmshFileName ) as gmshFileHandle:

        allLines = gmshFileHandle.readlines()

    for idx, lineval in enumerate( allLines ):
        if lineval == "$Nodes\n":

            nodeStart = idx + 2

    for idx, lineval in enumerate( allLines[ nodeStart: ] ):

        if lineval == "$EndNodes\n":
            break

        splitvals = lineval.split()
        vertexCoords[ splitvals[0] ] = [ float( val ) for val in splitvals[ 1: ] ]

    for idx, lineval in enumerate( allLines ):
        if lineval == "$Elements\n":

            elementStart = idx + 2
            nElements = int( allLines[ idx + 1 ] )

    for idx, lineval in enumerate( allLines[ elementStart: ] ):

        if lineval == "$EndElements\n":
            break

        splitvals = lineval.split()
        ntags = int( splitvals[ 2 ] )

        etype = splitvals[1] 

        if etype2dim[ etype ] < 3:
            vertexIndices = splitvals[ ntags + 3: ]

            if not checkElementOnBoundary( vertexIndices, vertexCoords ):

                lineIndicesToDelete.add( elementStart + idx )

    nElements = nElements - len( lineIndicesToDelete )
    allLines[ elementStart - 1 ] = str( nElements ) + "\n"
    newlines = []

    for idx, lineval in enumerate( allLines ):

        if idx not in lineIndicesToDelete:

            newlines.append( lineval )

    with open( gmshFileName, "w" ) as gmshFileHandle:

        gmshFileHandle.writelines( newlines )

def addLine( nodeBitVals, isPointPresent, curIdx, linesIdx, isPresentBitVals, allLines, allPts, ptEdgeMap, dxn ):

    offset = 3**dxn

    if nodeBitVals[ dxn ] != 2:
        if isPointPresent:    

            nextPoint = curIdx + offset
            
            if isPresentBitVals[ nextPoint ]:
                
                allLines[ linesIdx ] = gmsh.model.occ.addLine( allPts[ curIdx ], allPts[ nextPoint ] )
                
                ptEdgeMap[ (curIdx, dxn, 0) ] = allLines[ linesIdx ]
                ptEdgeMap[ (nextPoint, dxn, 1) ] = -allLines[ linesIdx ]

            elif nodeBitVals[ dxn ] == 0 and isPresentBitVals[ nextPoint + offset ] and 1 not in nodeBitVals:
                allLines[ linesIdx ] = gmsh.model.occ.addLine( allPts[ curIdx ], allPts[ nextPoint + offset ] )

                ptEdgeMap[ (curIdx, dxn, 0) ] = allLines[ linesIdx ]
                ptEdgeMap[ (nextPoint + offset, dxn, 1) ] = -allLines[ linesIdx ]

        return 1
    
    else:
        return 0

if __name__ == "__main__":

    gmshFileName = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/regularMesh3D_lvl0.msh"

    parseGMSHFile( gmshFileName )





    