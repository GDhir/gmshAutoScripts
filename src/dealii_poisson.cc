/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */

// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <cmath>
// This one is new. We want to read a triangulation from disk, and the class
// which does this is declared in the following file:
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
// We will use a circular domain, and the object describing the boundary of it
// comes from this file:
#include <deal.II/grid/manifold_lib.h>

// This is C++ ...
#include <fstream>
#include <iostream>
#include <vector>
#include <deal.II/base/hdf5.h>
#include <filesystem>
#include <regex>

// Finally, this has been discussed in previous tutorial programs before:
using namespace dealii;
using numbers::PI;
namespace fs = std::filesystem;

template <int dim>
class RightHandSide : public Function<dim>
{
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;

  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;

    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
    {
      std::cout << pair.first << '(' << pair.second << " times) ";
    }
    std::cout << std::endl;
  }

  std::ofstream out(filename);
  GridOut grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << filename << std::endl
            << std::endl;
}

// @sect3{The <code>Step5</code> class template}

// The main class is mostly as in the previous example. The most visible
// change is that the function <code>make_grid</code> has been
// removed, since creating the grid is now done in the <code>run</code>
// function and the rest of its functionality is now in
// <code>setup_system</code>. Apart from this, everything is as before.
template <int dim>
class Step5
{
public:
  Step5(FiniteElement<dim> &finiteElement, const Quadrature<dim> &quadObj);
  void run(std::string filename);

private:
  void setup_system();
  void assemble_system();
  void solve();
  void output_results(std::string meshFileName) const;

  Triangulation<dim> triangulation;
  // FE_Q<dim>          fe;
  const MappingFE<2> mapping;
  FiniteElement<2> &fe;
  const Quadrature<2> &quadrature_formula;
  DoFHandler<dim> dof_handler;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};

template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
  double return_value = 0.0;

  if (dim == 2)
  {
    return_value = -8 * std::pow(PI, 2) * sin(2 * PI * p[0]) * sin(2 * PI * p[1]);
  }

  return return_value;
}

template <int dim>
double BoundaryValues<dim>::value(const Point<dim> &p,
                                  const unsigned int /*component*/) const
{

  double return_value = 0.0;

  if (dim == 2)
  {
    return_value = sin(2 * PI * p[0]) * sin(2 * PI * p[1]);
  }

  return return_value;
}

// @sect3{The <code>Step5</code> class implementation}

// @sect4{Step5::Step5}

// This function is as before.
// template <int dim>
// Step5<dim>::Step5()
//   : fe(1)
//   , dof_handler(triangulation)
// {}
template <int dim>
Step5<dim>::Step5(FiniteElement<dim> &finiteElement, const Quadrature<dim> &quadObj)
    : mapping(finiteElement), fe(finiteElement), quadrature_formula(quadObj), dof_handler(triangulation)
{
}

// @sect4{Step5::setup_system}

// This is the function <code>make_grid</code> from the previous
// example, minus the generation of the grid. Everything else is unchanged:
template <int dim>
void Step5<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

// @sect4{Step5::assemble_system}

// As in the previous examples, this function is not changed much with regard
// to its functionality, but there are still some optimizations which we will
// show. For this, it is important to note that if efficient solvers are used
// (such as the preconditioned CG method), assembling the matrix and right hand
// side can take a comparable time, and you should think about using one or
// two optimizations at some places.
//
// The first parts of the function are completely unchanged from before:
template <int dim>
void Step5<dim>::assemble_system()
{
  // QGauss<dim> quadrature_formula(fe.degree + 1);
  RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Next is the typical loop over all cells to compute local contributions
  // and then to transfer them into the global matrix and vector. The only
  // change in this part, compared to step-4, is that we will use the
  // <code>coefficient()</code> function defined above to compute the
  // coefficient value at each quadrature point.
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    cell_matrix = 0.;
    cell_rhs = 0.;

    fe_values.reinit(cell);

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
    {

      for (const unsigned int i : fe_values.dof_indices())
      {
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
              (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
               fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
               fe_values.JxW(q_index));           // dx

        const auto &x_q = fe_values.quadrature_point(q_index);
        cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                        right_hand_side.value(x_q) *        // f(x_q)
                        fe_values.JxW(q_index));            // dx
      }
    }

    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices())
    {
      for (const unsigned int j : fe_values.dof_indices())
        system_matrix.add(local_dof_indices[i],
                          local_dof_indices[j],
                          cell_matrix(i, j));

      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  // With the matrix so built, we use zero boundary values again:
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(mapping,
                                           dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           boundary_values);

  MatrixTools::apply_boundary_values(
      boundary_values,
      system_matrix,
      solution,
      system_rhs);
}

// @sect4{Step5::solve}

// The solution process again looks mostly like in the previous
// examples. However, we will now use a preconditioned conjugate gradient
// algorithm. It is not very difficult to make this change. In fact, the only
// thing we have to alter is that we need an object which will act as a
// preconditioner. We will use SSOR (symmetric successive overrelaxation),
// with a relaxation factor of 1.2. For this purpose, the
// <code>SparseMatrix</code> class has a function which does one SSOR step,
// and we need to package the address of this function together with the
// matrix on which it should act (which is the matrix to be inverted) and the
// relaxation factor into one object. The <code>PreconditionSSOR</code> class
// does this for us. (<code>PreconditionSSOR</code> class takes a template
// argument denoting the matrix type it is supposed to work on. The default
// value is <code>SparseMatrix@<double@></code>, which is exactly what we need
// here, so we simply stick with the default and do not specify anything in
// the angle brackets.)
//
// Note that for the present case, SSOR doesn't really perform much better
// than most other preconditioners (though better than no preconditioning at
// all). A brief comparison of different preconditioners is presented in the
// Results section of the next tutorial program, step-6.
//
// With this, the rest of the function is trivial: instead of the
// <code>PreconditionIdentity</code> object we have created before, we now use
// the preconditioner we have declared, and the CG solver will do the rest for
// us:
template <int dim>
void Step5<dim>::solve()
{
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}

// @sect4{Step5::output_results and setting output flags}

// Writing output to a file is mostly the same as for the previous tutorial.
// The only difference is that we now need to construct a different filename
// for each refinement cycle.
//
// The function writes the output in VTU format, a variation of the VTK format
// that requires less disk space because it compresses the data. Of course,
// there are many other formats supported by the DataOut class if you
// desire to use a program for visualization that doesn't understand
// VTK or VTU.
template <int dim>
void Step5<dim>::output_results(std::string meshFileName) const
{

  std::string textfoldername{"/media/gaurav/easystore/dealii/MixedElement/TextFiles/"};

  std::regex regexVal("triangleMeshStruct|triangleMeshUnstruct|triangle|regular|mixed",
                      std::regex_constants::ECMAScript | std::regex_constants::icase);

  std::smatch meshTypeMatch;
  std::regex_search(meshFileName, meshTypeMatch, regexVal);

  std::string prefixVal = meshTypeMatch.str();
  std::string solutionFileName;

  std::vector<std::string> regularMeshTypes{"triangleMeshStruct", "triangleMeshUnstruct", "triangle", "regular"};

  if (std::find(regularMeshTypes.begin(), regularMeshTypes.end(), prefixVal) != regularMeshTypes.end())
  {

    std::regex re("N=");
    std::smatch m;
    std::regex_search(meshFileName, m, re);
    std::string suffixVal = m.suffix().str();

    re = "[0-9]+";
    std::regex_search(suffixVal, m, re);
    auto Nval = m.str();

    prefixVal = prefixVal + "_";

    solutionFileName = textfoldername + prefixVal + "solutionvaluesGaussModified_N=" + Nval;
    std::cout << solutionFileName << std::endl;
  }
  else if (prefixVal == "mixed")
  {

    /* To be implemented later */
  }

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(mapping);

  std::ofstream output(solutionFileName + ".vtu");
  data_out.write_vtu(output);

  const std::string filename_h5 = solutionFileName + ".h5";
  DataOutBase::DataOutFilterFlags flags(true, true);
  DataOutBase::DataOutFilter data_filter(flags);
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter, filename_h5, MPI_COMM_WORLD);
}

// @sect4{Step5::run}

// The second to last thing in this program is the definition of the
// <code>run()</code> function. In contrast to the previous programs, we will
// compute on a sequence of meshes that after each iteration is globally
// refined. The function therefore consists of a loop over 6 cycles. In each
// cycle, we first print the cycle number, and then have to decide what to do
// with the mesh. If this is not the first cycle, we simply refine the
// existing mesh once globally. Before running through these cycles, however,
// we have to generate a mesh:

// In previous examples, we have already used some of the functions from the
// <code>GridGenerator</code> class. Here we would like to read a grid from a
// file where the cells are stored and which may originate from someone else,
// or may be the product of a mesh generator tool.
//
// In order to read a grid from a file, we generate an object of data type
// GridIn and associate the triangulation to it (i.e. we tell it to fill our
// triangulation object when we ask it to read the file). Then we open the
// respective file and initialize the triangulation with the data in the file:
template <int dim>
void Step5<dim>::run(std::string filename)
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);

  std::ifstream input_file(filename);
  // We would now like to read the file. However, the input file is only for a
  // two-dimensional triangulation, while this function is a template for
  // arbitrary dimension. Since this is only a demonstration program, we will
  // not use different input files for the different dimensions, but rather
  // quickly kill the whole program if we are not in 2d. Of course, since the
  // main function below assumes that we are working in two dimensions we
  // could skip this check, in this version of the program, without any ill
  // effects.
  //
  // It turns out that more than 90 per cent of programming errors are invalid
  // function parameters such as invalid array sizes, etc, so we use
  // assertions heavily throughout deal.II to catch such mistakes. For this,
  // the <code>Assert</code> macro is a good choice, since it makes sure that
  // the condition which is given as first argument is valid, and if not
  // throws an exception (its second argument) which will usually terminate
  // the program giving information where the error occurred and what the
  // reason was. (A longer discussion of what exactly the @p Assert macro
  // does can be found in the @ref Exceptions "exception documentation module".)
  // This generally reduces the time to find programming errors
  // dramatically and we have found assertions an invaluable means to program
  // fast.
  //
  // On the other hand, all these checks (there are over 10,000 of them in the
  // library at present) should not slow down the program too much if you want
  // to do large computations. To this end, the <code>Assert</code> macro is
  // only used in debug mode and expands to nothing if in optimized
  // mode. Therefore, while you test your program on small problems and debug
  // it, the assertions will tell you where the problems are. Once your
  // program is stable, you can switch off debugging and the program will run
  // your real computations without the assertions and at maximum speed. More
  // precisely: turning off all the checks in the library (which prevent you
  // from calling functions with wrong arguments, walking off of arrays, etc.)
  // by compiling your program in optimized mode usually makes things run
  // about four times faster. Even though optimized programs are more
  // performant, we still recommend developing in debug mode since it allows
  // the library to find lots of common programming errors automatically. For
  // those who want to try: The way to switch from debug mode to optimized
  // mode is to recompile your program with the command <code>make
  // release</code>. The output of the <code>make</code> program should now
  // indicate to you that the program is now compiled in optimized mode, and
  // it will later also be linked to libraries that have been compiled for
  // optimized mode. In order to switch back to debug mode, simply recompile
  // with the command <code>make debug</code>.
  Assert(dim == 2, ExcInternalError());
  // ExcInternalError is a globally defined exception, which may be thrown
  // whenever something is terribly wrong. Usually, one would like to use more
  // specific exceptions, and particular in this case one would of course try
  // to do something else if <code>dim</code> is not equal to two, e.g. create
  // a grid using library functions. Aborting a program is usually not a good
  // idea and assertions should really only be used for exceptional cases
  // which should not occur, but might due to stupidity of the programmer,
  // user, or someone else. The situation above is not a very clever use of
  // Assert, but again: this is a tutorial and it might be worth to show what
  // not to do, after all.

  // So if we got past the assertion, we know that dim==2, and we can now
  // actually read the grid. It is in UCD (unstructured cell data) format
  // (though the convention is to use the suffix <code>inp</code> for UCD
  // files):
  grid_in.read_msh(input_file);
  // grid_in.read_ucd(input_file);
  // If you like to use another input format, you have to use one of the other
  // <code>grid_in.read_xxx</code> function. (See the documentation of the
  // <code>GridIn</code> class to find out what input formats are presently
  // supported.)

  // The grid in the file describes a circle. Therefore we have to use a
  // manifold object which tells the triangulation where to put new points on
  // the boundary when the grid is refined. Unlike step-1, since GridIn does
  // not know that the domain has a circular boundary (unlike
  // GridGenerator::hyper_shell) we have to explicitly attach a manifold to
  // the boundary after creating the triangulation to get the correct result
  // when we refine the mesh.
  // const SphericalManifold<dim> boundary;
  // triangulation.set_all_manifold_ids_on_boundary(0);
  // triangulation.set_manifold(0, boundary);

  print_mesh_info(triangulation, "grid.vtu");
  // std::cout << "   Number of active cells: "  //
  //           << triangulation.n_active_cells() //
  //           << std::endl                      //
  //           << "   Total number of cells: "   //
  //           << triangulation.n_cells()        //
  //           << std::endl;

  setup_system();
  assemble_system();
  solve();
  output_results(filename);
}

// template<typename dim>
void runFEM(FiniteElement<2> &fem, const Quadrature<2> &quadObj, std::string filename)
{

  Step5<2> laplace_problem_2d(fem, quadObj);
  laplace_problem_2d.run(filename);
}

// @sect3{The <code>main</code> function}

// The main function looks mostly like the one in the previous example, so we
// won't comment on it further:
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::string path = "/home/gaurav/Finch/src/examples/Mesh/MeshRun";
  std::string filename = "/home/gaurav/gmshAutoScripts/build/triangleMeshStructN=25.msh";

  std::regex tri_regex("triangle",
                       std::regex_constants::ECMAScript | std::regex_constants::icase);

  std::regex quad_regex("regular",
                        std::regex_constants::ECMAScript | std::regex_constants::icase);

  int idx = 0

      for (const auto &filenameObj : fs::directory_iterator(path))
  {

    filename = filenameObj.path();

    if (std::regex_search(filename, tri_regex))
    {
      FE_SimplexP<2> fem(1);
      const QGaussSimplex<2> quadrature_formula(fem.degree + 1);

      runFEM(fem, quadrature_formula, filename);

      const std::vector<Point<dim>> quadraturePoints = quadrature_formula.get_points();
      const std::vector<double> quadratureWeights = quadrature_formula.get_weights();

      std::cout << "idx = " << idx << " filename = " + filename << "\n";

      std::ofstream fileHandle("pointValuesTriangle_idx=" + std::to_string(idx) + ".txt", std::ios::out);

      fileHandle << filename << "\n";
      fileHandle << "Points Start Here \n";

      for (auto &pt : quadraturePoints)
      {

        fileHandle << pt[0] << "\t" << pt[1] << "\n";
      }

      fileHandle << "Points End Here \n";

      fileHandle << "Weights Start Here \n";

      for (auto &wt : quadratureWeights)
      {

        fileHandle << wt << "\n";
      }

      fileHandle << "Weights End Here \n";
    }
    else if (std::regex_search(filename, quad_regex))
    {
      FE_Q<2> fem(1);
      const QGauss<2> quadrature_formula(fem.degree + 1);
      runFEM(fem, quadrature_formula, filename);

      const std::vector<Point<dim>> quadraturePoints = quadrature_formula.get_points();
      const std::vector<double> quadratureWeights = quadrature_formula.get_weights();

      std::cout << "idx = " << idx << " filename = " + filename << "\n";

      std::ofstream fileHandle("pointValuesQuad_idx=" + std::to_string(idx) + ".txt", std::ios::out);

      fileHandle << filename << "\n";
      fileHandle << "Points Start Here \n";

      for (auto &pt : quadraturePoints)
      {

        fileHandle << pt[0] << "\t" << pt[1] << "\n";
      }

      fileHandle << "Points End Here \n";

      fileHandle << "Weights Start Here \n";

      for (auto &wt : quadratureWeights)
      {

        fileHandle << wt << "\n";
      }

      fileHandle << "Weights End Here \n";
    }

    idx += 1;

    // std::cout << filename << std::endl;
    // laplace_problem_2d.run( filename );
  }
  return 0;
}