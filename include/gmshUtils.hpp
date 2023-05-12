#ifndef utils
#define utils "utils"

#include <vector>
#include <gmsh.h>

void createBox(double xoffset, double yoffset, std::vector<int> &pts, std::vector<int> &lines, double lc, int Nx, int Ny);

void createMidObjects(std::vector<int> &linesmid, std::vector<int> &midcurves, std::vector<int> &midplanes,
                      std::vector<int> &ptsleft, std::vector<int> &ptsright, std::vector<int> &linesleft, std::vector<int> &linesright,
                      int Ncoarse, int Nx, std::pair<int, int> leftOffset, std::pair<int, int> rightOffset, bool leftCoarse, bool dxn);

void setTransfiniteCurves( std::vector< std::vector<int> >& linesVec );

void setTransfiniteSurfaces( std::vector<int>& planeIds, std::vector< std::vector<int> >& ptsVec,
     std::vector<int>& Nxvals, std::vector<int>& Nyvals );

void recombineSurfaces( std::vector<int>& planeIds );

void connectRefinedRegion( int Nx, int Ny, double xoffset, double yoffset, double lc, 
    std::vector<int>& botPts, std::vector<int>& topPts, std::vector<int>& botLinesLeft, std::vector<int>& topLinesLeft, 
    std::vector<int>& pts5, std::vector<int>& lines5, std::vector< std::vector<int> >& connectLines );

#endif