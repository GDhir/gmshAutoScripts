#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"
#include <fstream>

int main(int argc, char **argv)
{
    // gmsh::model::add("t11");

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    // int Nx1{0};
    std::vector<int> Nvals = { 11, 17 };
    // std::vector<int> Nx2vals = {9, 17, 33, 65, 129};
    // std::vector<int> Nx2vals = {3};
    
    int Nx2{0}, Nx4{0};
    std::ofstream outhandle{ "outfile.txt" };

    std::string foldername{ argv[1] };

    for( auto& N: Nvals ) {
        gmsh::initialize();
        gmsh::model::add("hangingMeshv4");


        int Nx3 = N;
        int Ny4 = N - 2;
        int Nx1 = N;
        int Ny = 3*N - 4;

        Nx2 = (Ny + Ny4)/2 - (Nx3 - 1)/2 - (Nx1 - 1)/2 - 1;
        Nx4 = Nx1 + 2*Nx2 + 1 + Nx3;
        
        // std::string filenamesuffix = 

        double dx = 1.0 / (Ny + Ny4);
        outhandle << 2*dx << "\n";
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

        std::vector<int> lines4;
        std::vector<int> pts4;

        xoffset = 0;
        yoffset = (Ny - 1)*lc + 2*lc;

        createBox(xoffset, yoffset, pts4, lines4, lc, Nx4, Ny4 );

        xoffset = 0;
        yoffset = (Ny - 1)*lc;
        // std::vector<int> connectPtsLeft;
        std::vector< std::vector<int> > connectLinesLeft( 2, std::vector<int>() );

        std::vector<int> topPtsLeft( Nx1, 0 );
        std::vector<int> botPtsLeft( Nx1, 0 );

        std::vector<int> topLinesLeft;
        std::vector<int> botLinesLeft;

        for( int i = 0; i < Nx1; i++ ) {
            botPtsLeft[i] = pts1[ Nx1 + Ny - 1 + Nx1 - 1 - 1 - i ];
            topPtsLeft[i] = pts4[ i ];
        }

        for( int i = 0; i < Nx1 - 1; i++ ) {
            botLinesLeft.push_back( -lines1[ Nx1 - 1 + Ny - 1 + Nx1 - 1 - 1 - i ] );
        }

        for( int i = 0; i < Nx1 - 1; i++ ) {
            topLinesLeft.push_back( -lines4[ Nx1 - 1 - 1 - i ] );
        }

        std::vector<int> pts5;
        std::vector<int> lines5;

        connectRefinedRegion( Nx1, Ny, xoffset, yoffset, lc, 
            botPtsLeft, topPtsLeft, botLinesLeft, topLinesLeft, pts5, lines5, connectLinesLeft );    

        // connect right refined region

        xoffset = (Nx1 - 1)*lc + 2*lc + (Nx2 - 1)*2*lc + 2*lc;
        yoffset = (Ny - 1)*lc;
        // std::vector<int> connectPtsLeft;
        std::vector< std::vector<int> > connectLinesRight( 2, std::vector<int>() );

        std::vector<int> topPtsRight( Nx3, 0 );
        std::vector<int> botPtsRight( Nx3, 0 );

        std::vector<int> topLinesRight;
        std::vector<int> botLinesRight;

        for( int i = 0; i < Nx3; i++ ) {
            botPtsRight[i] = pts3[ Nx3 + Ny - 1 + Nx3 - 1 - 1 - i ];
            topPtsRight[i] = pts4[ Nx1 + 2*Nx2 + 1 + i ];
        }

        for( int i = 0; i < Nx3 - 1; i++ ) {
            botLinesRight.push_back( -lines3[ Nx3 - 1 + Ny - 1 + Nx3 - 1 - 1 - i ] );
        }

        for( int i = 0; i < Nx3 - 1; i++ ) {
            topLinesRight.push_back( -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 2 + 2 + Nx3 - 1 - 1 - i ] );
        }

        std::vector<int> pts6;
        std::vector<int> lines6;

        connectRefinedRegion( Nx3, Ny, xoffset, yoffset, lc, 
            botPtsRight, topPtsRight, botLinesRight, topLinesRight, pts6, lines6, connectLinesRight );    

        std::vector<int> linesmid1;

        std::vector<int> midcurves1;
        std::vector<int> midplanes1;

        std::vector<int> ptsleft;
        std::vector<int> ptsright;

        std::vector<int> linesleft;
        std::vector<int> linesright;

        for( int i = 0; i < Ny - 1; i++ ) {

            linesleft.push_back( lines1[ Nx1 - 1 + i ] );

        }

        for( int i = 0; i < Nright - 1; i++ ) {

            linesright.push_back( -lines2[  lines2.size() - 1 - i ] );

        }

        for( int i = 0; i < Ny; i += 2 ) {

            ptsleft.push_back( pts1[ Nx1 + i - 1 ] );

        }

        ptsright.push_back( pts2[ 0 ] );
        for( int i = 0; i < Nright - 1; i += 1 ) {

            ptsright.push_back( pts2[ pts2.size() - 1 - i ] );

        }

        createMidObjects(linesmid1, midcurves1, midplanes1, ptsleft, ptsright, linesleft, linesright, Nright, Nx1, false, true);

        std::vector<int> linesmid2;

        std::vector<int> midcurves2;
        std::vector<int> midplanes2;

        linesleft.clear();
        linesright.clear();
        ptsleft.clear();
        ptsright.clear();

        for( int i = 0; i < Nright - 1; i++ ) {

            linesleft.push_back( lines2[ Nx2 - 1 + i ] );
            // std::cout << linesleft[i] << "\n";

        }

        // std::cout << "lines3 \n";

        for( int i = 0; i < lines3.size(); i++ ) {
            // std::cout << lines3[i] << "\n";
        }

        // std::cout << "linesright \n";

        for( int i = 0; i < Ny - 1; i++ ) {

            linesright.push_back( -lines3[  lines3.size() - 1 - i ] );
            // std::cout << linesright[i] << "\n";
        }

        for( int i = 0; i < Nright; i += 1 ) {

            ptsleft.push_back( pts2[ Nx2 + i - 1 ] );

        }

        ptsright.push_back( pts3[ 0 ] );
        for( int i = 0; i < Ny - 1; i += 2 ) {

            ptsright.push_back( pts3[ pts3.size() - 2 - i ] );

        }

        // createMidObjects(linesmid2, midcurves2, midplanes2, pts2, pts3, lines2, lines3, Nright, Nx2, std::make_pair(Nx2 - 1, 0), std::make_pair(0, 0), true, true);
        createMidObjects(linesmid2, midcurves2, midplanes2, ptsleft, ptsright, linesleft, linesright, Nright, Nx2, true, true);

        // new

        std::vector<int> linesmid3;

        std::vector<int> midcurves3;
        std::vector<int> midplanes3;

        linesleft.clear();
        linesright.clear();
        ptsleft.clear();
        ptsright.clear();

        for( int i = 0; i < 2*Nx2 - 2; i++ ) {

            linesleft.push_back( lines4[ Nx1 - 1 + 2 + i ] );
            // std::cout << linesleft[i] << "\n";

        }

        // std::cout << "lines3 \n";

        // for( int i = 0; i < lines3.size(); i++ ) {
        //     std::cout << lines3[i] << "\n";
        // }

        // std::cout << "linesright \n";

        for( int i = 0; i < Nx2 - 1; i++ ) {

            linesright.push_back( -lines2[  lines2.size() - ( Nright - 1 ) - 1 - i ] );
            // std::cout << linesright[i] << "\n";
        }

        for( int i = 0; i < 2*Nx2 - 1; i += 2 ) {

            ptsleft.push_back( pts4[ Nx1 + 2 + i - 1 ] );

        }

        // ptsright.push_back( pts3[ 0 ] );
        for( int i = 0; i < Nx2; i += 1 ) {

            ptsright.push_back( pts2[ pts2.size() - ( Nright - 1) - i ] );

        }

        // createMidObjects(linesmid2, midcurves2, midplanes2, pts2, pts3, lines2, lines3, Nright, Nx2, std::make_pair(Nx2 - 1, 0), std::make_pair(0, 0), true, true);
        createMidObjects(linesmid3, midcurves3, midplanes3, ptsleft, ptsright, linesleft, linesright, Nx2, Nx2, false, true);

        int c1 = gmsh::model::geo::addCurveLoop(lines1);
        int p1 = gmsh::model::geo::addPlaneSurface({c1});

        int c2 = gmsh::model::geo::addCurveLoop(lines2);
        int p2 = gmsh::model::geo::addPlaneSurface({c2});

        int c3 = gmsh::model::geo::addCurveLoop(lines3);
        int p3 = gmsh::model::geo::addPlaneSurface({c3});

        int c4 = gmsh::model::geo::addCurveLoop(lines4);
        int p4 = gmsh::model::geo::addPlaneSurface({c4});

        int c5 = gmsh::model::geo::addCurveLoop(lines5);
        int p5 = gmsh::model::geo::addPlaneSurface({c5});

        int c6 = gmsh::model::geo::addCurveLoop(lines6);
        int p6 = gmsh::model::geo::addPlaneSurface({c6});

        int c7 = gmsh::model::geo::addCurveLoop( { linesmid1.back(), -linesmid3[0], -lines4[ Nx1 ], -lines4[ Nx1 - 1 ],
        -connectLinesLeft[1][1], -connectLinesLeft[1][0] } );
        int p7 = gmsh::model::geo::addPlaneSurface({c7});

        int c8 = gmsh::model::geo::addCurveLoop( { linesmid2.back(), connectLinesRight[0][0], connectLinesRight[0][1],
        -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 1 ], -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 2 ], linesmid3.back() } );
        int p8 = gmsh::model::geo::addPlaneSurface({c8});

        std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, lines4, linesmid1, linesmid2, linesmid3,
        connectLinesLeft[0], connectLinesLeft[1],connectLinesRight[0], connectLinesRight[1] };
        // std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, linesmid2 };
        setTransfiniteCurves( linesVec );

        std::vector<int> planeIds{ p1, p2, p3, p4, p5, p6};
        std::vector< std::vector<int> > ptsVec{ pts1, pts2, pts3, pts4, pts5, pts6};
        std::vector<int> Nxvals{ Nx1, Nx2, Nx3, Nx4, Nx1, Nx3 };  
        std::vector<int> Nyvals{ Ny, Nright, Ny, Ny4, 3, 3 };

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
        gmsh::model::mesh::setOrder( 1 );

        gmsh::model::mesh::generate(2);

        // Note that you could also apply the recombination algorithm and/or the
        // subdivision step explicitly after meshing, as follows:
        //
        // gmsh::model::mesh::generate(2);
        // gmsh::model::mesh::recombine();
        // gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
        // gmsh::model::mesh::refine();
        gmsh::option::setNumber("Mesh.MshFileVersion", 2);
        std::string meshfilename = foldername + "hangingMeshNx=" + std::to_string( Nx2 ) + "Ny=" + std::to_string( Ny ) + ".msh";
        gmsh::write(meshfilename);

        // Launch the GUI to see the results:
        std::set<std::string> args(argv, argv + argc);
        if (!args.count("-nopopup"))
            gmsh::fltk::run();

        gmsh::finalize();
    }
    
    return 0;
}
