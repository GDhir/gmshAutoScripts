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

    gmsh::model::add("t11");

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

    // gmsh::model::geo::addPoint(1.25, 0.25, 0, lc/2);
    // gmsh::model::geo::addPoint(1.25, 0.75, 0, lc/2);

    std::vector<int> linesmid;

    // linesmid.push_back( gmsh::model::geo::addLine( ptsleft[ N - 1 ], ptsright[ 0 ] ) );

    // for( int i = 1; i < Nright; i++ ) {
    //     linesmid.push_back( gmsh::model::geo::addLine( ptsleft[ (N - 1) + i*2 ], ptsright[ ptsright.size() - i ] ) );
    //     // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
    // }

    std::vector< int > midcurves;
    std::vector< int > midplanes;

    createMidObjects(linesmid, midcurves, midplanes, ptsleft, ptsright, linesleft, linesright, Nright, N);

    // for( int i = 0; i < Nright - 1; i++ ) {

    //     midcurves.push_back( gmsh::model::geo::addCurveLoop( { linesmid[i], -*( linesright.end() - i - 1 ), -linesmid[i + 1], -linesleft[ N - 1 + i*2 + 1 ], -linesleft[ N - 1 + i*2 ] } ) );
    //     midplanes.push_back( gmsh::model::geo::addPlaneSurface({midcurves.back()}) );

    // }

    int cl = gmsh::model::geo::addCurveLoop( linesleft );
    int pl = gmsh::model::geo::addPlaneSurface({cl});

    int cr = gmsh::model::geo::addCurveLoop( linesright );
    int pr = gmsh::model::geo::addPlaneSurface({cr});

    // int cmid0 = gmsh::model::geo::addCurveLoop( linesmid0 );
    // int pmid0 = gmsh::model::geo::addPlaneSurface({cmid0});

    // int cmid1 = gmsh::model::geo::addCurveLoop( linesmid1 );
    // int pmid1 = gmsh::model::geo::addPlaneSurface({cmid1});

    // for( int i = 0; i < linesleft.size(); i++ ) {
    //     gmsh::model::geo::mesh::setTransfiniteCurve( linesleft[i], 2 );
    // }
    // for( int i = 0; i < linesright.size(); i++ ) {
    //     gmsh::model::geo::mesh::setTransfiniteCurve( linesright[i], 2 );
    // }
    // for( int i = 0; i < linesmid.size(); i++ ) {
    //     gmsh::model::geo::mesh::setTransfiniteCurve( linesmid[i], 2 );
    // }

    // gmsh::model::geo::mesh::setTransfiniteSurface(pl, "Left", {ptsleft[0], ptsleft[ (N - 1) ], ptsleft[ (N - 1)*2 ], ptsleft[ (N - 1)*3 ]});
    // gmsh::model::geo::mesh::setTransfiniteSurface(pr, "Left", {ptsright[0], ptsright[ (Nright - 1) ], ptsright[ (Nright - 1)*2 ], ptsright[ (Nright - 1)*3 ]});

    std::vector< std::vector<int> > linesVec{ linesleft, linesright, linesmid };
    setTransfiniteCurves( linesVec );

    std::vector<int> planeIds{ pl, pr };
    std::vector< std::vector<int> > ptsVec{ ptsleft, ptsright };
    std::vector<int> Nxvals{ N, N };  
    std::vector<int> Nyvals{ N, Nright };

    setTransfiniteSurfaces( planeIds, ptsVec, Nxvals, Nyvals );

    gmsh::model::geo::synchronize();

    // To generate quadrangles instead of triangles, we can simply add
    // gmsh::model::mesh::setRecombine(2, pr);
    // gmsh::model::mesh::setRecombine(2, pl);
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
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
