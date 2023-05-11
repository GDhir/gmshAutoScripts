#ifndef utils
#define utils "utils"

#include <vector>
#include <gmsh.h>

void createBox(double xoffset, double yoffset, std::vector<int> &pts, std::vector<int> &lines, double lc, int Nx, int Ny);

void createMidObjects(std::vector<int> &linesmid, std::vector<int> &midcurves, std::vector<int> &midplanes,
                      std::vector<int> &ptsleft, std::vector<int> &ptsright, std::vector<int> &linesleft,
                       std::vector<int> &linesright, int Nright, int Nx);

void setTransfiniteCurves( std::vector< std::vector<int> >& linesVec );

void setTransfiniteSurfaces( std::vector<int>& planeIds, std::vector< std::vector<int> >& ptsVec,
     std::vector<int>& Nxvals, std::vector<int>& Nyvals );

void recombineSurfaces( std::vector<int>& planeIds );

#endif