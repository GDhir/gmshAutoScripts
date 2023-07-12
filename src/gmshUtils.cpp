#include "gmshUtils.hpp"
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

void connectRefinedRegion(int Nx, int Ny, double xoffset, double yoffset, double lc,
                          std::vector<int> &botPts, std::vector<int> &topPts, std::vector<int> &botLines, std::vector<int> &topLines,
                          std::vector<int> &pts5, std::vector<int> &lines5, std::vector<std::vector<int>> &connectLines)
{

    double x = xoffset;
    double y = yoffset + lc;
    std::vector<int> connectPts;

    connectPts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));

    x = xoffset + (Nx - 1) * lc;
    connectPts.push_back(gmsh::model::geo::addPoint(x, y, 0, lc));

    // int botIdx = Nx + Ny - 1 + Nx - 1 - 1;
    // int topIdx = 0;

    connectLines[0].push_back(gmsh::model::geo::addLine(botPts[0], connectPts[0]));
    connectLines[0].push_back(gmsh::model::geo::addLine(connectPts[0], topPts[0]));

    // // botIdx = Nx + Ny - 1 - 1;
    // // topIdx = Nx - 1;

    connectLines[1].push_back(gmsh::model::geo::addLine(botPts[Nx - 1], connectPts[1]));
    connectLines[1].push_back(gmsh::model::geo::addLine(connectPts[1], topPts[Nx - 1]));

    pts5.insert(pts5.end(), botPts.begin(), botPts.end());
    pts5.push_back(connectPts[1]);
    pts5.insert(pts5.end(), topPts.rbegin(), topPts.rend());
    pts5.push_back(connectPts[0]);

    lines5.insert(lines5.end(), botLines.begin(), botLines.end());
    lines5.insert(lines5.end(), connectLines[1].begin(), connectLines[1].end());
    lines5.insert(lines5.end(), topLines.begin(), topLines.end());

    for (int i = 0; i < 2; i++)
    {
        lines5.push_back(-connectLines[0][1 - i]);
    }

    // std::cout << "begin \n";
    // for( auto& val: pts5 ) {
    //     std::cout << val << "\n";
    // }
    // std::cout << "end \n";
}

void createMidObjects(std::vector<int> &linesmid, std::vector<int> &midcurves, std::vector<int> &midplanes,
                      std::vector<int> &ptsleft, std::vector<int> &ptsright, std::vector<int> &linesleft, std::vector<int> &linesright,
                      int Ncoarse, int Nx, std::pair<int, int> leftOffset, std::pair<int, int> rightOffset, bool leftCoarse, bool dxn)
{

    int xleft{leftOffset.first};
    int xright{rightOffset.first};
    int yleft{leftOffset.second};
    int yright{rightOffset.second};

    if (dxn)
    {
        if (leftCoarse)
        {

            linesmid.push_back(gmsh::model::geo::addLine(ptsleft[xleft], ptsright[xright]));

            for (int i = 1; i < Ncoarse; i++)
            {
                linesmid.push_back(gmsh::model::geo::addLine(ptsleft[xleft + i], ptsright[ptsright.size() - i * 2]));
                // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
            }

            for (int i = 0; i < Ncoarse - 1; i++)
            {

                std::cout << linesmid[i] << "\t" << -*(linesright.end() - i * 2 - 1) << "\t" << -*(linesright.end() - i * 2 - 2) << "\t" << -linesmid[i + 1] << "\t" << -linesleft[xleft + i] << "\n";

                midcurves.push_back(gmsh::model::geo::addCurveLoop({linesmid[i], -*(linesright.end() - i * 2 - 1), -*(linesright.end() - i * 2 - 2), -linesmid[i + 1], -linesleft[xleft + i]}));
                midplanes.push_back(gmsh::model::geo::addPlaneSurface({midcurves.back()}));

                // gmsh::model::geo::mesh::setTransfiniteSurface(midplanes[i], "Left", {ptsleft[xleft + i], ptsright[ptsright.size() - i*2 ], ptsright[ptsright.size() - i*2 - 2], ptsleft[xleft + i + 1]});
            }
            // gmsh::model::geo::mesh::setTransfiniteSurface(midplanes[0], "Left", {ptsleft[xleft], ptsright[0], ptsright[ptsright.size() - 2], ptsleft[xleft + 1]});
        }
        else
        {
            linesmid.push_back(gmsh::model::geo::addLine(ptsleft[xleft], ptsright[xright]));

            for (int i = 1; i < Ncoarse; i++)
            {
                linesmid.push_back(gmsh::model::geo::addLine(ptsleft[xleft + i * 2], ptsright[ptsright.size() - i]));
                // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
            }

            for (int i = 0; i < Ncoarse - 1; i++)
            {

                midcurves.push_back(gmsh::model::geo::addCurveLoop({linesmid[i], -*(linesright.end() - i - 1), -linesmid[i + 1], -linesleft[xleft + i * 2 + 1], -linesleft[xleft + i * 2]}));
                midplanes.push_back(gmsh::model::geo::addPlaneSurface({midcurves.back()}));
            }
        }
    }
}

void createMidObjects(std::vector<int> &linesmid, std::vector<int> &midcurves, std::vector<int> &midplanes,
                      std::vector<int> &ptsleft, std::vector<int> &ptsright, std::vector<int> &linesleft, std::vector<int> &linesright,
                      int Ncoarse, int Nx, bool leftCoarse, bool dxn)
{

    if(dxn) {

        if( leftCoarse ) {

            // std::cout << "othermid \n";

            for (int i = 0; i < Ncoarse; i++)
            {
                linesmid.push_back(gmsh::model::geo::addLine( ptsleft[i], ptsright[i] ) );
                // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
            }

            for (int i = 0; i < Ncoarse - 1; i++)
            {

                // std::cout << linesmid[i] << "\t" << linesright[2*i] << "\t" << linesright[2*i + 1] << "\t" << -linesmid[i + 1] << "\t" << -linesleft[i] << "\n";

                midcurves.push_back(gmsh::model::geo::addCurveLoop({linesright[ 2*i ], linesright[ 2*i + 1 ], -linesmid[i + 1], -linesleft[i], linesmid[i]}));
                midplanes.push_back(gmsh::model::geo::addPlaneSurface({midcurves.back()}));
                
                // gmsh::model::geo::mesh::setTransfiniteSurface(midplanes[i], "Left", {ptsleft[xleft + i], ptsright[ptsright.size() - i*2 ], ptsright[ptsright.size() - i*2 - 2], ptsleft[xleft + i + 1]});
            }
        }
        else {

            for (int i = 0; i < Ncoarse; i++)
            {
                linesmid.push_back(gmsh::model::geo::addLine( ptsleft[i], ptsright[i] ) );
                // std::cout << ptsleft[ (N - 1) + i*2 ] << "\t" << ptsright[ ptsright.size() - i - 1 ] << "\n";
            }

            for (int i = 0; i < Ncoarse - 1; i++)
            {

                midcurves.push_back(gmsh::model::geo::addCurveLoop({linesmid[i], linesright[ i ], -linesmid[i + 1], -linesleft[2*i + 1], -linesleft[2*i] } ));
                midplanes.push_back(gmsh::model::geo::addPlaneSurface({midcurves.back()}));

                // gmsh::model::geo::mesh::setTransfiniteSurface(midplanes[i], "Left", {ptsleft[xleft + i], ptsright[ptsright.size() - i*2 ], ptsright[ptsright.size() - i*2 - 2], ptsleft[xleft + i + 1]});
            }

        }

    }

}

void setTransfiniteCurves(std::vector<std::vector<int>> &linesVec)
{

    for (auto &lines : linesVec)
    {
        for (int i = 0; i < lines.size(); i++)
        {
            gmsh::model::geo::mesh::setTransfiniteCurve(lines[i], 2);
        }
    }
}

void setTransfiniteSurfaces(std::vector<int> &planeIds, std::vector<std::vector<int>> &ptsVec,
                            std::vector<int> &Nxvals, std::vector<int> &Nyvals)
{

    int id{0};
    std::vector<int> pts;
    int Nx{0};
    int Ny{0};

    for (int i = 0; i < planeIds.size(); i++)
    {

        id = planeIds[i];
        pts = ptsVec[i];

        Nx = Nxvals[i];
        Ny = Nyvals[i];

        gmsh::model::geo::mesh::setTransfiniteSurface(id, "Left", {pts[0], pts[(Nx - 1)], pts[(Nx - 1) + (Ny - 1)], pts[(Nx - 1) * 2 + Ny - 1]});
    }
}

void recombineSurfaces(std::vector<int> &planeIds)
{

    for (auto &id : planeIds)
    {
        gmsh::model::mesh::setRecombine(2, id);
    }
}