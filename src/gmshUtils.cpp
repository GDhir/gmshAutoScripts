#include "gmshUtils.hpp"

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