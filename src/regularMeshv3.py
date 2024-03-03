# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 17
#
#  Anisotropic background mesh
#
# ------------------------------------------------------------------------------

# As seen in `t7.py', mesh sizes can be specified very accurately by providing a
# background mesh, i.e., a post-processing view that contains the target mesh
# sizes.

# Here, the background mesh is represented as a metric tensor field defined on a
# square. One should use bamg as 2d mesh generator to enable anisotropic meshes
# in 2D.

import sys
sys.path.append( "/home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/" )
import gmsh
import math
import os
import sys

def occHexBoxMesh( folderName ):

    lvl = 0
    Ndashvals = [10, 20, 30, 40]

    for Ndash in Ndashvals:

        lc = 1/( Ndash )
        print(lc)

        N = int( 1/lc + 1 )

        gmsh.initialize()
        gmsh.model.add("t17")

        # Create a square
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

        gmsh.model.mesh.setTransfiniteAutomatic()

        # Merge a post-processing view containing the target anisotropic mesh sizes
        # path = os.path.dirname(os.path.abspath(__file__))
        # gmsh.merge(os.path.join(path, os.pardir, 't17_bgmesh.pos'))

        # # Apply the view as the current background mesh
        # bg_field = gmsh.model.mesh.field.add("PostView")
        # gmsh.model.mesh.field.setNumber(bg_field, "ViewIndex", 0)
        # gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

        # Use bamg
        # gmsh.option.setNumber("Mesh.SmoothRatio", 30)
        # gmsh.option.setNumber("Mesh.AnisoMax", 10)
        gmsh.option.setNumber("Mesh.Algorithm", 5)

        gmsh.model.mesh.generate(3)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2)

        regMeshFileName = folderName + "regularMesh3D_lvl" + str(lvl) + ".msh"
        gmsh.write( regMeshFileName )

        # Launch the GUI to see the results:
        # if '-nopopup' not in sys.argv:
        #     gmsh.fltk.run()

        gmsh.finalize()

        lvl = lvl + 1

def occTetBoxMesh( folderName ):

    Ndashvals = [8, 12, 16, 20]

    for lvl, Ndash in enumerate( Ndashvals ):

        lc = 1/( Ndash )
        print(lc)

        N = int( 1/lc + 1 )

        gmsh.initialize(sys.argv)
        gmsh.model.add("t2")

        xoffset = 0
        yoffset = 0
        pts = []
        lines = []
        Nx = N
        Ny = N

        # Create a square
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

        # Merge a post-processing view containing the target anisotropic mesh sizes
        # path = os.path.dirname(os.path.abspath(__file__))
        # gmsh.merge(os.path.join(path, os.pardir, 't17_bgmesh.pos'))

        # # Apply the view as the current background mesh
        # bg_field = gmsh.model.mesh.field.add("PostView")
        # gmsh.model.mesh.field.setNumber(bg_field, "ViewIndex", 0)
        # gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

        # Use bamg
        gmsh.option.setNumber("Mesh.SmoothRatio", 30)
        gmsh.option.setNumber("Mesh.AnisoMax", 10)
        gmsh.option.setNumber("Mesh.Algorithm", 4)

        gmsh.model.mesh.generate(3)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2)

        regMeshFileName = folderName + "tetMesh3D_lvl" + str(lvl) + ".msh"
        gmsh.write( regMeshFileName )

        # Launch the GUI to see the results:
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()

        gmsh.finalize()

if __name__ == "__main__":

    folderName = "/home/gaurav/gmshAutoScripts/build/"
    occHexBoxMesh( folderName )