// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 11
//
//  Unstructured quadrangular meshes
//
// -----------------------------------------------------------------------------

#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"

int main(int argc, char **argv)
{
    gmsh::initialize();

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    int N = 21;
    double dx = 1.0/( N - 1 );
    double lc = dx;
    int Nright = N/2 + 1;
    
    std::vector<int> linesleft;
    std::vector<int> ptsleft;

    double xoffset = 0;
    double yoffset = 0;

    createBox( xoffset, yoffset, ptsleft, linesleft, lc, N, N );

    std::vector<int> linesright;
    std::vector<int> ptsright;

    xoffset = (N - 1)*dx + 2*dx;
    yoffset = 0;

    createBox( xoffset, yoffset, ptsright, linesright, lc*2, N, N/2 + 1 );

    std::vector<int> linesmid;
    std::vector< int > midcurves;
    std::vector< int > midplanes;

    createMidObjects(linesmid, midcurves, midplanes, ptsleft, ptsright, linesleft, linesright, Nright, N, std::make_pair(N - 1, 0), std::make_pair(0, 0), false, true);

    int cl = gmsh::model::geo::addCurveLoop( linesleft );
    int pl = gmsh::model::geo::addPlaneSurface({cl});

    int cr = gmsh::model::geo::addCurveLoop( linesright );
    int pr = gmsh::model::geo::addPlaneSurface({cr});

    std::vector< std::vector<int> > linesVec{ linesleft, linesright, linesmid };
    setTransfiniteCurves( linesVec );

    std::vector<int> planeIds{ pl, pr };
    std::vector< std::vector<int> > ptsVec{ ptsleft, ptsright };
    std::vector<int> Nxvals{ N, N };  
    std::vector<int> Nyvals{ N, Nright };

    setTransfiniteSurfaces( planeIds, ptsVec, Nxvals, Nyvals );

    gmsh::model::geo::synchronize();

    // To generate quadrangles instead of triangles, we can simply add
    recombineSurfaces( planeIds );

    // If we'd had several surfaces, we could have used the global option
    // "Mesh.RecombineAll":
    //
    // gmsh::option::setNumber("Mesh.RecombineAll", 1);

    // The default recombination algorithm is called "Blossom": it uses a minimum
    // cost perfect matching algorithm to generate fully quadrilateral meshes from
    // triangulations. More details about the algorithm can be found in the
    // following paper: J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise,
    // A. Johnen and C. Geuzaine, "Blossom-Quad: a non-uniform quadrilateral mesh
    // generator using a minimum cost perfect matching algorithm", International
    // Journal for Numerical Methods in Engineering 89, pp. 1102-1119, 2012.

    // For even better 2D (planar) quadrilateral meshes, you can try the
    // experimental "Frontal-Delaunay for quads" meshing algorithm, which is a
    // triangulation algorithm that enables to create right triangles almost
    // everywhere: J.-F. Remacle, F. Henrotte, T. Carrier-Baudouin, E. Bechet,
    // E. Marchandise, C. Geuzaine and T. Mouton. A frontal Delaunay quad mesh
    // generator using the L^inf norm. International Journal for Numerical Methods
    // in Engineering, 94, pp. 494-512, 2013. Uncomment the following line to try
    // the Frontal-Delaunay algorithms for quads:
    //
    gmsh::option::setNumber("Mesh.Algorithm", 8);

    // The default recombination algorithm might leave some triangles in the mesh,
    // if recombining all the triangles leads to badly shaped quads. In such
    // cases, to generate full-quad meshes, you can either subdivide the resulting
    // hybrid mesh (with `Mesh.SubdivisionAlgorithm' set to 1), or use the
    // full-quad recombination algorithm, which will automatically perform a
    // coarser mesh followed by recombination, smoothing and
    // subdivision. Uncomment the following line to try the full-quad algorithm:
    //
    gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // or 3

    // You can also set the subdivision step alone, with
    //
    // gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);

    gmsh::model::mesh::generate(2);

    // Note that you could also apply the recombination algorithm and/or the
    // subdivision step explicitly after meshing, as follows:
    //
    // gmsh::model::mesh::generate(2);
    // gmsh::model::mesh::recombine();
    // gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
    // gmsh::model::mesh::refine();

    // Launch the GUI to see the results:
    gmsh::write("hangingMeshv1.msh");
    // gmsh::view::write(t1, "TextFiles/hangingMeshv1.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
