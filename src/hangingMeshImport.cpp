#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"

int main(int argc, char **argv)
{

    gmsh::initialize();

    std::string foldername = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/";
    std::string filename = foldername + "regularMesh3D_lvl0.msh";
    // std::string filename = "importMesh.msh";
    // You can run this tutorial on any file that Gmsh can read, e.g. a mesh file
    // in the MSH format: `t1.exe file.msh'
    gmsh::open( filename );

    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();


    return 0;
}