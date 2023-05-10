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

void createBox(double xoffset, double yoffset, std::vector<int> &pts, std::vector<int> &lines, double lc, int Nx, int Ny)
{

    int side = 1;
    int i = 0;
    double x, y;

    while (i < Nx)
    {
        x = xoffset + lc * i;
        y = yoffset + 0.0;
        // std::cout << x << "\t" << y << "\n";
        // pts.push_back( gmsh::model::geo::addPoint(x, y, 0, lc, pttagoffset + (N - 1)*side + i) );
        pts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));
        i++;
    }

    side = 2;

    int j = 1;

    while (i < Nx + Ny - 1)
    {
        x = xoffset + (Nx - 1) * lc;
        y = yoffset + j * lc;
        // std::cout << x << "\t" << y << "\n";
        // pts.push_back( gmsh::model::geo::addPoint( x, y, 0, lc, pttagoffset + (N - 1)*side + i) );
        pts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));

        j++;
        i++;
    }

    side = 3;
    j = 1;

    while (i < (Nx - 1) * 2 + Ny)
    {
        x = xoffset + (Nx - 1) * lc - j * lc;
        y = yoffset + (Ny - 1) * lc;
        // std::cout << x << "\t" << y << "\n";
        // pts.push_back( gmsh::model::geo::addPoint( x, y, 0, lc, pttagoffset + (N - 1)*side + i) );
        pts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));

        j++;
        i++;
    }

    side = 4;
    j = 1;

    while (i < (Nx - 1) * 2 + (Ny - 1) * 2)
    {
        x = xoffset + 0;
        y = yoffset + (Ny - 1) * lc - j * lc;
        // std::cout << x << "\t" << y << "\n";
        // pts.push_back( gmsh::model::geo::addPoint( 0, (N - 1)*lc - j*lc, 0, lc, pttagoffset + (N - 1)*side + i) );
        pts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));
        j++;
        i++;
    }

    for (int i = 0; i < pts.size() - 1; i++)
    {

        lines.push_back(gmsh::model::geo::addLine(pts[i], pts[i + 1]));
    }

    lines.push_back(gmsh::model::geo::addLine(pts.back(), pts[0]));
}

void createMidObjects(std::vector<int> &linesmid, std::vector<int> &midcurves, std::vector<int> &midplanes,
                      std::vector<int> &ptsleft, std::vector<int> &ptsright, std::vector<int> &linesleft, std::vector<int> &linesright, int Nright, int Nx)
{

    linesmid.push_back(gmsh::model::geo::addLine(ptsleft[Nx - 1], ptsright[0]));

    for (int i = 1; i < Nright; i++)
    {
        linesmid.push_back(gmsh::model::geo::addLine(ptsleft[(Nx - 1) + i * 2], ptsright[ptsright.size() - i]));
        // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
    }

    for (int i = 0; i < Nright - 1; i++)
    {

        midcurves.push_back(gmsh::model::geo::addCurveLoop({linesmid[i], -*(linesright.end() - i - 1), -linesmid[i + 1], -linesleft[Nx - 1 + i * 2 + 1], -linesleft[Nx - 1 + i * 2]}));
        midplanes.push_back(gmsh::model::geo::addPlaneSurface({midcurves.back()}));
    }
}

void setTransfiniteCurves( std::vector< std::vector<int> >& linesVec ) {

    for( auto& lines: linesVec ) {
        for (int i = 0; i < lines.size(); i++)
        {
            gmsh::model::geo::mesh::setTransfiniteCurve(lines[i], 2);
        }
    }

}

void setTransfiniteSurfaces( std::vector<int>& planeIds, std::vector< std::vector<int> >& ptsVec,
     std::vector<int>& Nxvals, std::vector<int>& Nyvals ) {

    int id{0};
    std::vector<int> pts;
    int Nx{0};
    int Ny{0};

    for( int i = 0; i < planeIds.size(); i++ ) {    

        id = planeIds[i];
        pts = ptsVec[i];

        Nx = Nxvals[i];
        Ny = Nyvals[i];

        gmsh::model::geo::mesh::setTransfiniteSurface(id, "Left", {pts[0], pts[(Nx - 1)], pts[(Nx - 1) + (Ny - 1)], pts[(Nx - 1) * 2 + Ny - 1]});

    }

}

void recombineSurfaces( std::vector<int>& planeIds ) {

    for( auto& id: planeIds ) {
        gmsh::model::mesh::setRecombine(2, id);
    }
}

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("t11");

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    int Nx = 5;
    int Ny = 21;
    double dx = 1.0 / (Ny - 1);
    double lc = dx;
    int Nright = Ny / 2 + 1;

    std::vector<int> linesleft;
    std::vector<int> ptsleft;

    double xoffset = 0;
    double yoffset = 0;

    createBox(xoffset, yoffset, ptsleft, linesleft, lc, Nx, Ny);

    std::vector<int> linesright;
    std::vector<int> ptsright;

    xoffset = (Nx - 1) * dx + 2 * dx;
    yoffset = 0;

    createBox(xoffset, yoffset, ptsright, linesright, lc * 2, Nx, Ny / 2 + 1);

    // gmsh::model::geo::addPoint(1.25, 0.25, 0, lc/2);
    // gmsh::model::geo::addPoint(1.25, 0.75, 0, lc/2);

    std::vector<int> linesmid;

    std::vector<int> midcurves;
    std::vector<int> midplanes;

    createMidObjects(linesmid, midcurves, midplanes, ptsleft, ptsright, linesleft, linesright, Nright, Nx);

    int cl = gmsh::model::geo::addCurveLoop(linesleft);
    int pl = gmsh::model::geo::addPlaneSurface({cl});

    int cr = gmsh::model::geo::addCurveLoop(linesright);
    int pr = gmsh::model::geo::addPlaneSurface({cr});

    std::vector< std::vector<int> > linesVec{ linesleft, linesright, linesmid };
    setTransfiniteCurves( linesVec );

    std::vector<int> planeIds{ pl, pr };
    std::vector< std::vector<int> > ptsVec{ ptsleft, ptsright };
    std::vector<int> Nxvals{ Nx, Nx };  
    std::vector<int> Nyvals{ Ny, Nright };

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
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
