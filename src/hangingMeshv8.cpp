#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"
#include <fstream>

// void connectRefined( std::vector< std::vector<int> >& connectLines, std::vector<int>& pts, std::vector<int>& lines,
//      std::vector<int>& ptsleft, std::vector<int>& ptsright, std::vector<int>& linesleft, std::vector<int>& linesright,
//     int xoffset, int yoffset, int Nxleft, int Nyleft, int Nxright, int Nyright, double lc  ) {

//         std::vector<int> ptsleftStored( Nxleft, 0 );
//         std::vector<int> ptsrightStored( Nxright, 0 );

//         std::vector<int> linesleftStored;
//         std::vector<int> linesrightStored;

//         for( int i = 0; i < Nxright; i++ ) {
//             ptsrightStored[i] = ptsright[ Nxright + Nyright - 1 + Nxright - 1 - 1 - i ];
//             ptsleftStored[i] = ptsleft[ i ];
//         }

//         for( int i = 0; i < Nxright - 1; i++ ) {
//             linesrightStored.push_back( -linesright[ Nxright - 1 + Nyright - 1 + Nxright - 1 - 1 - i ] );
//         }

//         for( int i = 0; i < Nxright - 1; i++ ) {
//             linesleftStored.push_back( -linesleft[ Nxright - 1 - 1 - i ] );
//         }

//         connectRefinedRegion( Nxright, Nyright, xoffset, yoffset, lc, 
//             ptsrightStored, ptsleftStored, linesrightStored, linesleftStored, pts, lines, connectLines );    

// }

int main(int argc, char **argv)
{
    // gmsh::model::add("t11");

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    // std::vector<int> Nyvals = {33, 65, 129, 257, 513};
    // std::vector<int> Nyvals = {13, 15, 17};
    std::vector<int> Nyvals = { 9, 17, 33, 65, 129 };
    // std::vector<int> Nyvals = {17, 33};
    // std::vector<int> Nx2vals = {3};
    
    int N, Nx1, Nx3, Ny4, Ny5, Nx2{0}, Nx4{0}, Nx5{0};
    std::ofstream outhandle;
    std::string foldername{ argv[1] };

    for( auto& Ny: Nyvals ) {
        gmsh::initialize();
        gmsh::model::add("hangingMeshv8");
        
        N = (Ny + 5)/2;
        Nx1 = N;
        Nx3 = N;
        Ny4 = N;
        Ny5 = N;

        // Nx2 = (Ny + Ny4)/2 - (Nx3 - 1)/2 - (Nx1 - 1)/2 - 1;
        // Nx2 = (Ny + Ny4 + Ny5)/2 - (Nx3)/2 - (Nx1)/2;
        Nx2 = Ny/2 + 1;
        Nx4 = Nx1 + 2*Nx2 + 1 + Nx3;
        Nx5 = Nx1 + 2*Nx2 + 1 + Nx3;
        
        // std::string filenamesuffix = 
        outhandle.open( "outfile.txt", std::ios::app );
        double dx = 1.0 / (Ny + Ny4 + Ny5 + 1);
        outhandle << 2*dx << "\n";
        outhandle << Nx1 << "\t" << Nx2 << "\t" << Nx3 << "\t" << Nx4 << "\t" << Nx5 << "\n";
        outhandle << Ny << "\t" << Ny4 << "\t" << Ny5 << "\n";
        outhandle << "end \n";
        outhandle.close();

        double lc = dx;
        int Nright = Ny / 2 + 1;

        std::vector<int> lines1;
        std::vector<int> pts1;

        double xoffset = 0;
        double yoffset = (Ny5 + 1)*lc;

        createBox(xoffset, yoffset, pts1, lines1, lc, Nx1, Ny);

        std::vector<int> lines2;
        std::vector<int> pts2;

        xoffset = (Nx1 - 1) * lc + 2 * lc;
        yoffset = (Ny5 + 1)*lc;

        createBox(xoffset, yoffset, pts2, lines2, lc * 2, Nx2, Ny / 2 + 1);

        std::vector<int> lines3;
        std::vector<int> pts3;

        xoffset = (Nx1 - 1) * lc + 2 * lc + (Nx2 - 1) * 2 * lc + 2 * lc;
        yoffset = (Ny5 + 1)*lc;

        createBox(xoffset, yoffset, pts3, lines3, lc, Nx3, Ny );

        // gmsh::model::geo::addPoint(1.25, 0.25, 0, lc/2);
        // gmsh::model::geo::addPoint(1.25, 0.75, 0, lc/2);

        std::vector<int> lines4;
        std::vector<int> pts4;

        xoffset = 0;
        yoffset = (Ny - 1)*lc + 2*lc + (Ny5 - 1)*lc + 2*lc;

        createBox(xoffset, yoffset, pts4, lines4, lc, Nx4, Ny4 );

        std::vector<int> linesbot;
        std::vector<int> ptsbot;

        xoffset = 0;
        yoffset = 0;

        createBox(xoffset, yoffset, ptsbot, linesbot, lc, Nx5, Ny5 );

        xoffset = 0;
        yoffset = (Ny - 1)*lc + (Ny5 - 1)*lc + 2*lc;
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
        yoffset = (Ny - 1)*lc + (Ny5 - 1)*lc + 2*lc;
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

        // Connect refined region bottom left

        xoffset = 0;
        yoffset = (Ny5 - 1)*lc;
        // std::vector<int> connectPtsLeft;
        std::vector< std::vector<int> > connectLinesBotLeft( 2, std::vector<int>() );

        // topPtsLeft.clear();
        // botPtsLeft.clear();

        topLinesLeft.clear();
        botLinesLeft.clear();

        for( int i = 0; i < Nx1; i++ ) {
            botPtsLeft[i] = ptsbot[ Nx5 + Ny5 - 1 + Nx5 - 1 - 1 - i ];
            topPtsLeft[i] = pts1[ i ];
        }

        for( int i = 0; i < Nx1 - 1; i++ ) {
            botLinesLeft.push_back( -linesbot[ Nx5 - 1 + Ny5 - 1 + Nx5 - 1 - 1 - i ] );
        }

        for( int i = 0; i < Nx1 - 1; i++ ) {
            topLinesLeft.push_back( -lines1[ Nx1 - 1 - 1 - i ] );
        }

        std::vector<int> ptsConnectBotLeft;
        std::vector<int> linesConnectBotLeft;

        connectRefinedRegion( Nx1, Ny, xoffset, yoffset, lc, 
            botPtsLeft, topPtsLeft, botLinesLeft, topLinesLeft, ptsConnectBotLeft, linesConnectBotLeft, connectLinesBotLeft );    

        // for( auto& val: ptsConnectBotLeft ) {
        //     std::cout << val << "\n";
        // }

        // connect bottom right refined region

        xoffset = (Nx1 - 1)*lc + 2*lc + (Nx2 - 1)*2*lc + 2*lc;
        yoffset = (Ny5 - 1)*lc;
        // std::vector<int> connectPtsLeft;
        std::vector< std::vector<int> > connectLinesBotRight( 2, std::vector<int>() );

        // topPtsRight.clear();
        // botPtsRight.clear();

        topLinesRight.clear();
        botLinesRight.clear();

        for( int i = 0; i < Nx3; i++ ) {
            botPtsRight[i] = ptsbot[ Nx5 - 1 + Ny5 - 1 + Nx3 - 1 - i ];
            topPtsRight[i] = pts3[ i ];
        }

        for( int i = 0; i < Nx3 - 1; i++ ) {
            botLinesRight.push_back( -linesbot[ Nx5 - 1 + Ny5 - 1 + Nx3 - 1 - 1 - i ] );
        }

        for( int i = 0; i < Nx3 - 1; i++ ) {
            topLinesRight.push_back( -lines3[ Nx3 - 1 - 1 - i ] );
        }

        std::vector<int> ptsConnectBotRight;
        std::vector<int> linesConnectBotRight;

        connectRefinedRegion( Nx3, Ny5, xoffset, yoffset, lc, 
            botPtsRight, topPtsRight, botLinesRight, topLinesRight, ptsConnectBotRight, linesConnectBotRight, connectLinesBotRight );    

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

        // new bottom mid objects

        std::vector<int> linesmid4;

        std::vector<int> midcurves4;
        std::vector<int> midplanes4;

        linesleft.clear();
        linesright.clear();
        ptsleft.clear();
        ptsright.clear();

        for( int i = 0; i < Nx2 - 1; i++ ) {

            linesleft.push_back( lines2[ i ] );
            // std::cout << linesleft[i] << "\n";

        }

        // std::cout << "lines3 \n";

        // for( int i = 0; i < lines3.size(); i++ ) {
        //     std::cout << lines3[i] << "\n";
        // }

        // std::cout << "linesright \n";

        for( int i = 0; i < 2*Nx2 - 2; i++ ) {

            linesright.push_back( -linesbot[  linesbot.size() - ( Ny5 - 1 ) - (Nx1 - 1) - 2 - 1 - i ] );
            // std::cout << linesright[i] << "\n";
        }
        // std::cout << "end \n";
        // for( int i = 0; i < 2*Nx2 - 1; i += 2 ) {

        //     ptsleft.push_back( pts4[ Nx1 + 2 + i - 1 ] );

        // }

        for( int i = 0; i < Nx2; i += 1 ) {

            ptsleft.push_back( pts2[ i ] );

        }

        // ptsright.push_back( pts3[ 0 ] );
        for( int i = 0; i < 2*Nx2 - 1; i += 2 ) {

            ptsright.push_back( ptsbot[ ptsbot.size() - ( Ny5 - 1) - Nx1 - 1 - i ] );
            // std::cout << Nx2 << "\n";
        }

        // createMidObjects(linesmid2, midcurves2, midplanes2, pts2, pts3, lines2, lines3, Nright, Nx2, std::make_pair(Nx2 - 1, 0), std::make_pair(0, 0), true, true);
        createMidObjects(linesmid4, midcurves4, midplanes4, ptsleft, ptsright, linesleft, linesright, Nx2, Nx2, true, true);

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

        int cbot = gmsh::model::geo::addCurveLoop(linesbot);
        int pbot = gmsh::model::geo::addPlaneSurface({cbot});

        int cbotleft = gmsh::model::geo::addCurveLoop(linesConnectBotLeft);
        int pbotleft = gmsh::model::geo::addPlaneSurface({cbotleft});

        int cbotright = gmsh::model::geo::addCurveLoop(linesConnectBotRight);
        int pbotright = gmsh::model::geo::addPlaneSurface({cbotright});

        int c7 = gmsh::model::geo::addCurveLoop( { linesmid1.back(), -linesmid3[0], -lines4[ Nx1 ], -lines4[ Nx1 - 1 ],
        -connectLinesLeft[1][1], -connectLinesLeft[1][0] } );
        int p7 = gmsh::model::geo::addPlaneSurface({c7});

        int c8 = gmsh::model::geo::addCurveLoop( { linesmid2.back(), connectLinesRight[0][0], connectLinesRight[0][1],
        -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 1 ], -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 2 ], linesmid3.back() } );
        int p8 = gmsh::model::geo::addPlaneSurface({c8});

        std::cout << -linesbot[linesbot.size() - (Ny5 - 1) - (Nx1 - 1) - 1 ] << "\t" <<
         -linesbot[linesbot.size() - (Ny5 - 1) - (Nx1 - 1) - 2 ] << "\t" << 
         -linesmid4[0] << "\t" << -connectLinesBotLeft[1][1] << "\t" << -connectLinesBotLeft[1][0] << "\n";

        int c9 = gmsh::model::geo::addCurveLoop( { -linesbot[linesbot.size() - (Ny5 - 1) - (Nx1 - 1) - 1 ],
         -linesbot[linesbot.size() - (Ny5 - 1) - (Nx1 - 1) - 2 ],
         -linesmid4[0], -linesmid1[0], -connectLinesBotLeft[1][1], -connectLinesBotLeft[1][0] });
        int p9 = gmsh::model::geo::addPlaneSurface({c9});

        std::cout << -linesbot[(Nx5 - 1) + (Ny5 - 1) + (Nx3 - 1) + 1 ] << "\t" <<
         -linesbot[(Nx5 - 1) + (Ny5 - 1) + (Nx3 - 1) ] << "\t" <<
         connectLinesBotRight[0][0] << "\t" << connectLinesBotRight[0][1] << "\t" <<
         -linesmid2[0] << "\t" << linesmid4.back() << "\n";

        int c10 = gmsh::model::geo::addCurveLoop( { -linesbot[(Nx5 - 1) + (Ny5 - 1) + (Nx3 - 1) + 1 ],
         -linesbot[(Nx5 - 1) + (Ny5 - 1) + (Nx3 - 1) ],
         connectLinesBotRight[0][0], connectLinesBotRight[0][1],
         -linesmid2[0], linesmid4.back() });
        int p10 = gmsh::model::geo::addPlaneSurface({c10});

        // int c10 = gmsh::model::geo::addCurveLoop( { linesmid2.back(), connectLinesRight[0][0], connectLinesRight[0][1],
        // -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 1 ], -lines4[ Nx1 - 1 + 2 + 2*Nx2 - 2 ], linesmid3.back() } );
        // int p10 = gmsh::model::geo::addPlaneSurface({c10});

        std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, lines4, linesbot, linesmid1, linesmid2, linesmid3, linesmid4,
        connectLinesLeft[0], connectLinesLeft[1],connectLinesRight[0], connectLinesRight[1],
        connectLinesBotLeft[0], connectLinesBotLeft[1], connectLinesBotRight[0], connectLinesBotRight[1]};
        // std::vector< std::vector<int> > linesVec{ lines1, lines2, lines3, linesmid2 };
        setTransfiniteCurves( linesVec );

        std::vector<int> planeIds{ p1, p2, p3, p4, p5, p6, pbot, pbotleft, pbotright};
        std::vector< std::vector<int> > ptsVec{ pts1, pts2, pts3, pts4, pts5, pts6, ptsbot, ptsConnectBotLeft, ptsConnectBotRight};
        std::vector<int> Nxvals{ Nx1, Nx2, Nx3, Nx4, Nx1, Nx3, Nx5, Nx1, Nx3};  
        std::vector<int> Nyvals{ Ny, Nright, Ny, Ny4, 3, 3, Ny5, 3, 3};

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
        gmsh::option::setNumber("Mesh.MshFileVersion", 2);
        std::string meshfilename = foldername + "hangingMeshv8Nx=" + std::to_string( Nx2 ) + "Ny=" + std::to_string( Ny ) + ".msh";
        gmsh::write(meshfilename);

        // Launch the GUI to see the results:
        // std::set<std::string> args(argv, argv + argc);
        // if (!args.count("-nopopup"))
        //     gmsh::fltk::run();

        gmsh::finalize();
    }
    
    return 0;
}
