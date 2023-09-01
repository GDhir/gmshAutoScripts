#include <set>
#include <gmsh.h>
#include <iostream>
#include "gmshUtils.hpp"

int main(int argc, char **argv)
{

    gmsh::initialize();

    std::string filename = "hangingMeshv7Nx=8Ny=29.msh";
    // std::string filename = "importMesh.msh";
    // You can run this tutorial on any file that Gmsh can read, e.g. a mesh file
    // in the MSH format: `t1.exe file.msh'
    gmsh::open( filename );

    gmsh::merge( filename );

    int order = 2;
    gmsh::model::mesh::setOrder(order);

    // gmsh::model::mesh::generate(2);

    filename = "hangingMeshv7Order=" + std::to_string(order) + "Nx=8Ny=29.msh";
    gmsh::write( filename );
    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}