#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("hangingMeshv2");


    // gmsh::model::add("t11");

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    int Nx1 = 3;
    int Nx2 = 7;
    int Nx3 = 11;
    int Ny = 21;
    double dx = 1.0 / (Ny - 1);
    double lc = dx;
    int Nright = Ny / 2 + 1;

    std::vector<int> lines1;
    std::vector<int> pts1;

    double xoffset = 0;
    double yoffset = 0;

    createBox(xoffset, yoffset, pts1, lines1, lc, Nx1, Ny);

    std::vector<int> lines2;
    std::vector<int> pts2;

    xoffset = (Nx1 - 1) * dx + 2 * dx;
    yoffset = 0;

    createBox(xoffset, yoffset, pts2, lines2, lc * 2, Nx2, Ny / 2 + 1);

    std::vector<int> lines3;
    std::vector<int> pts3;

    xoffset = (Nx1 - 1) * dx + 2 * dx + (Nx2 - 1) * 2 * dx + 2 * dx;
    yoffset = 0;

    createBox(xoffset, yoffset, pts3, lines3, lc, Nx3, Ny );

    // gmsh::model::geo::addPoint(1.25, 0.25, 0, lc/2);
    // gmsh::model::geo::addPoint(1.25, 0.75, 0, lc/2);

    std::vector<int> linesmid1;

    std::vector<int> midcurves1;
    std::vector<int> midplanes1;

    createMidObjects(linesmid1, midcurves1, midplanes1, pts1, pts2, lines1, lines2, Nright, Nx1, std::make_pair(Nx1 - 1, 0), std::make_pair(0, 0), false, true);

    std::vector<int> linesmid2;

    std::vector<int> midcurves2;
    std::vector<int> midplanes2;

    createMidObjects(linesmid2, midcurves2, midplanes2, pts2, pts3, lines2, lines3, Nright, Nx2, std::make_pair(Nx2 - 1, 0), std::make_pair(0, 0), true, true);

    int c1 = gmsh::model::geo::addCurveLoop(lines1);
    int p1 = gmsh::model::geo::addPlaneSurface({c1});

    int c2 = gmsh::model::geo::addCurveLoop(lines2);
    int p2 = gmsh::model::geo::addPlaneSurface({c2});

    int c3 = gmsh::model::geo::addCurveLoop(lines3);
    int p3 = gmsh::model::geo::addPlaneSurface({c3});

    std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, linesmid1, linesmid2 };
    // std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, linesmid2 };
    setTransfiniteCurves( linesVec );

    std::vector<int> planeIds{ p1, p2, p3 };
    std::vector< std::vector<int> > ptsVec{ pts1, pts2, pts3 };
    std::vector<int> Nxvals{ Nx1, Nx2, Nx3 };  
    std::vector<int> Nyvals{ Ny, Nright, Ny };

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
    gmsh::write("hangingMeshv3.msh");

    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
