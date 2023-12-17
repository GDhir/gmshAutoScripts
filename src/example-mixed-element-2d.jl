#=
# 2D Poisson, Dirichlet bc
# Using a mixed-element mesh (quads and triangles)
=#
### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################


functionSpace(order=2) # basis function polynomial order

function getDerivatives( num_elements, mesh, geometric_factors, nodes_per_element, refel, solution )

    qnodes_per_element = refel.Nqp

    RQ1::Matrix{Float64} = zeros(Float64, qnodes_per_element, nodes_per_element)
    RQ2::Matrix{Float64} = zeros(Float64, qnodes_per_element, nodes_per_element)

    deriv_vals = zeros( Float64, num_elements, 2 )
    centroid_vals = zeros( Float64, num_elements, 2 )

    for ei = 1:num_elements
        eid = mesh.elemental_order[ei]
        index_offset = 0
        build_derivative_matrix(refel, geometric_factors, 1, eid, 0, RQ1)
        build_derivative_matrix(refel, geometric_factors, 2, eid, 0, RQ2)
        #= Prepare derivative matrices. =#
        #= Evaluate coefficients. =#
        for ni = 1:nodes_per_element
            nodeID = mesh.loc2glb[ni, eid]
            x = mesh.allnodes[1, nodeID]
            y = mesh.allnodes[2, nodeID]
            
            centroid_vals[ ei, 1 ] += x
            centroid_vals[ ei, 2 ] += y

            deriv_vals[ ei, 1 ] += RQ1[ 1, ni ]*solution[ ni ]
            deriv_vals[ ei, 2 ] += RQ2[ 1, ni ]*solution[ ni ]

        end

        centroid_vals = centroid_vals/nodes_per_element
    end

    return deriv_vals, centroid_vals

end

function regularElement( k, c )

    meshvals = readdir( "/home/gaurav/Finch/src/examples/Mesh/MeshRun", join = true )

    for (index, meshval) in enumerate(meshvals)

        val = match( r"regular", meshval )

        if !isnothing( val )

            initFinch("mixed2d");

            useLog("mixed2dlog", level=3)

            domain(2)

            mesh(meshval)

            u = variable("u")
            err = variable("err")
            err = variable("err")
            testSymbol("v")
            exactval = variable("exact")
            # exact solution is sin(2*pi*x)*sin(2*pi*y)
            maxerr = 0;

            coeff = parse( Int64, k )
            exact(x,y) = sin(pi*x*coeff)*sin(pi*y*coeff);

            linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

            boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # boundaryStr = "sin(pi*x)*sin(pi*y)"
            boundary(u, 1, DIRICHLET, boundaryStr)

            strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # strval = "-pi*pi*sin(pi*x)*sin(pi*y)"
            coefficient("f", strval)

            weakForm(u, "dot(grad(u), grad(v)) - f*v")

            exportCode("mixed2dcode");
            # importCode("mixed2dcodein");

            solve(u);

            finalizeFinch()
            
            xy = Finch.finch_state.grid_data.allnodes;
            global maxerr;
            global L2Error = 0;

            for i=1:size(xy,2)
                x = xy[1,i];
                y = xy[2,i];
                err.values[i] = abs(u.values[i] - exact(x,y));
                exactval.values[i] = exact(x, y)
                global maxerr;
                maxerr = max(err.values[i], maxerr);
                global L2Error;
                L2Error += err.values[i]^2
            end
            println("max error = "*string(maxerr));
            println("L2 error = "*string(L2Error));

            # Uncomment below to plot or output result

            # using Plots
            # using PyPlot
            # # using Triplot
            # pyplot();
            # using DelimitedFiles

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # figure(2)
            # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
            # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
            # display(h)

            # a = u.values[:]

            Nval =  match( r"[0-9]+", meshval )
            Nval = Nval.match
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "regular_"
            file = open( textfoldername * prefix * "uvalues_lvl=" * Nval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_lvl=" * Nval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_lvl=" * Nval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            meshval = Finch.finch_state.grid_data
            num_elements = meshval.nel_owned
            refel = Finch.finch_state.refel
            nodes_per_element = refel.Np

            deriv_vals, centroid_vals = getDerivatives( num_elements, meshval, Finch.finch_state.geo_factors, nodes_per_element, refel, u.values[:] )

            file = open( textfoldername * prefix * "centroid_xvalues_lvl=" * Nval * ".txt", "w")
            print( file, centroid_vals[:, 1])

            close(file)

            file = open( textfoldername * prefix * "centroid_yvalues_lvl=" * Nval * ".txt", "w")
            print( file, centroid_vals[:, 2])

            close(file)

            file = open( textfoldername * prefix * "deriv_xvalues_lvl=" * Nval * ".txt", "w")
            print( file, deriv_vals[:, 1])

            close(file)

            file = open( textfoldername * prefix * "deriv_yvalues_lvl=" * Nval * ".txt", "w")
            print( file, deriv_vals[:, 2])

            close(file)

            outputValues( [u, err], textfoldername * prefix * "errorAndu_lvl=" * Nval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end

    end
end

function mixedElement( k, c )

    meshvals = readdir( "/home/gaurav/Finch/src/examples/Mesh/MeshRun", join = true )

    for (index, meshval) in enumerate(meshvals)

        val = match( r"hanging", meshval )

        if !isnothing( val )
            Nxstr =  match( r"Nx=", meshval )
            Nxoffset = Nxstr.offset + 3
            Nxval = match( r"[0-9]+", meshval[ Nxoffset:end ] )
            Nxval = Nxval.match

            Nystr =  match( r"Ny=", meshval )
            Nyoffset = Nystr.offset + 3
            Nyval = match( r"[0-9]+", meshval[ Nyoffset:end ] )
            Nyval = Nyval.match

            initFinch("mixed2d");

            useLog("mixed2dlog", level=3)

            domain(2)

            mesh(meshval)

            u = variable("u")
            err = variable("err")
            err = variable("err")
            testSymbol("v")
            exactval = variable("exact")
            # exact solution is sin(2*pi*x)*sin(2*pi*y)
            maxerr = 0;
            exact(x,y) = sin(pi*x*3)*sin(pi*y*3);
            linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )
            boundary(u, 1, DIRICHLET, "sin(pi*x*3)*sin(pi*y*3)")

            coefficient("f", "-18*pi*pi*sin(3*pi*x)*sin(3*pi*y)")

            weakForm(u, "dot(grad(u), grad(v)) - f*v")

            # exportCode("mixed2dcode");
            importCode("mixed2dcodein");

            solve(u);

            finalizeFinch()
            
            xy = Finch.finch_state.grid_data.allnodes;
            global maxerr;
            global L2Error = 0;

            for i=1:size(xy,2)
                x = xy[1,i];
                y = xy[2,i];
                err.values[i] = abs(u.values[i] - exact(x,y));
                exactval.values[i] = exact(x, y)
                global maxerr;
                maxerr = max(err.values[i], maxerr);
                global L2Error;
                L2Error += err.values[i]^2
            end
            println("max error = "*string(maxerr));
            println("L2 error = "*string(L2Error));

            # Uncomment below to plot or output result

            # using Plots
            # using PyPlot
            # # using Triplot
            # pyplot();
            # using DelimitedFiles

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # figure(2)
            # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
            # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
            # display(h)

            # a = u.values[:]
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "mixed_"
            file = open( textfoldername * prefix * "uvalues_Nx=" * Nxval * "Ny=" * Nyval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_Nx=" * Nxval * "Ny=" * Nyval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_Nx=" * Nxval * "Ny=" * Nyval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_Nx=" * Nxval * "Ny=" * Nyval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_Nx=" * Nxval * "Ny=" * Nyval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            outputValues( [u, err], textfoldername * prefix * "errorAnduNx=" * Nxval * "Ny=" * Nyval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end
    end
end

function mixedElementWithLevel( k, c )

    meshvals = readdir( "/home/gaurav/Finch/src/examples/Mesh/MeshRun/mix_mesh/", join = true )

    for (index, meshval) in enumerate(meshvals)

        val = match( r"lvl", meshval )

        if !isnothing( val )
            lvlstr =  match( r"lvl", meshval )
            lvloffset = lvlstr.offset + 3
            lvlval = match( r"[0-9]+", meshval[ lvloffset:end ] )
            lvlval = lvlval.match

            initFinch("mixed2d");

            useLog("mixed2dlog", level=3)

            domain(2)

            mesh(meshval)

            u = variable("u")
            err = variable("err")
            err = variable("err")
            testSymbol("v")
            exactval = variable("exact")
            # exact solution is sin(2*pi*x)*sin(2*pi*y)
            maxerr = 0;
            coeff = parse( Int64, k )
            coeff = parse( Int64, k )
            exact(x,y) = sin(pi*x*coeff)*sin(pi*y*coeff);

            linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

            boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # boundaryStr = "sin(pi*x)*sin(pi*y)"
            boundary(u, 1, DIRICHLET, boundaryStr)

            strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # strval = "pi*pi*sin(pi*x)*sin(pi*y)"
            coefficient("f", strval)

            weakForm(u, "dot(grad(u), grad(v)) - f*v")

            # exportCode("mixed2dcode");
            importCode("mixed2dcodein");

            solve(u);

            finalizeFinch()
            
            xy = Finch.finch_state.grid_data.allnodes;
            global maxerr;
            global L2Error = 0;

            for i=1:size(xy,2)
                x = xy[1,i];
                y = xy[2,i];
                err.values[i] = abs(u.values[i] - exact(x,y));
                exactval.values[i] = exact(x, y)
                global maxerr;
                maxerr = max(err.values[i], maxerr);
                global L2Error;
                L2Error += err.values[i]^2
            end
            println("max error = "*string(maxerr));
            println("L2 error = "*string(L2Error));

            # Uncomment below to plot or output result

            # using Plots
            # using PyPlot
            # # using Triplot
            # pyplot();
            # using DelimitedFiles

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # figure(2)
            # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
            # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
            # display(h)

            # a = u.values[:]
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "mesh_"
            file = open( textfoldername * prefix * "uvalues_lvl=" * lvlval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_lvl=" * lvlval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_lvl=" * lvlval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_lvl=" * lvlval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_lvl=" * lvlval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            outputValues( [u, err], textfoldername * prefix * "errorAndu_lvl=" * lvlval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end
    end
end

function triangleElement( k, c )

    meshvals = readdir( "/home/gaurav/Finch/src/examples/Mesh/MeshRun", join = true )

    for (index, meshval) in enumerate(meshvals)

        val = match( r"triangleMeshStruct", meshval )

        if !isnothing( val )

            initFinch("mixed2d");

            useLog("mixed2dlog", level=3)

            domain(2)

            mesh(meshval)

            u = variable("u")
            err = variable("err")
            err = variable("err")
            testSymbol("v")
            exactval = variable("exact")
            # exact solution is sin(3*pi*x)*sin(3*pi*y)
            maxerr = 0;
            coeff = parse( Int64, k )
            exact(x,y) = sin(pi*x*coeff)*sin(pi*y*coeff);

            linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

            boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # boundaryStr = "sin(pi*x)*sin(pi*y)"
            boundary(u, 1, DIRICHLET, boundaryStr)

            strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # strval = "pi*pi*sin(pi*x)*sin(pi*y)"
            coefficient("f", strval)

            weakForm(u, "dot(grad(u), grad(v)) - f*v")

            exportCode("mixed2dcode");
            # importCode("mixed2dcode");

            solve(u);

            finalizeFinch()
            
            xy = Finch.finch_state.grid_data.allnodes;
            global maxerr;
            global L2Error = 0;

            for i=1:size(xy,2)
                x = xy[1,i];
                y = xy[2,i];
                err.values[i] = abs(u.values[i] - exact(x,y));
                exactval.values[i] = exact(x, y)
                global maxerr;
                maxerr = max(err.values[i], maxerr);
                global L2Error;
                L2Error += err.values[i]^2
            end
            println("max error = "*string(maxerr));
            println("L2 error = "*string(L2Error));

            # Uncomment below to plot or output result

            # using Plots
            # using PyPlot
            # # using Triplot
            # pyplot();
            # using DelimitedFiles

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # figure(2)
            # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
            # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
            # display(h)

            # a = u.values[:]

            Nval =  match( r"[0-9]+", meshval )
            Nval = Nval.match
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "triangleMeshStruct_"
            file = open( textfoldername * prefix * "uvalues_lvl=" * Nval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_lvl=" * Nval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_lvl=" * Nval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            outputValues( [u, err], textfoldername * prefix * "errorAndu_lvl=" * Nval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end
    end

    meshvals = readdir( "/home/gaurav/Finch/src/examples/Mesh/MeshRun", join = true )

    for (index, meshval) in enumerate(meshvals)

        val = match( r"triangleMeshUnstruct", meshval )

        if !isnothing( val )

            initFinch("mixed2d");

            useLog("mixed2dlog", level=3)

            domain(2)

            mesh(meshval)

            u = variable("u")
            err = variable("err")
            err = variable("err")
            testSymbol("v")
            exactval = variable("exact")
            # exact solution is sin(2*pi*x)*sin(2*pi*y)
            maxerr = 0;
            coeff = parse( Int64, k )
            exact(x,y) = sin(pi*x*coeff)*sin(pi*y*coeff);

            linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

            boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # boundaryStr = "sin(pi*x)*sin(pi*y)"
            boundary(u, 1, DIRICHLET, boundaryStr)

            strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")"
            # strval = "pi*pi*sin(pi*x)*sin(pi*y)"
            coefficient("f", strval)

            weakForm(u, "dot(grad(u), grad(v)) - f*v")

            exportCode("mixed2dcode");
            # importCode("mixed2dcode");

            solve(u);

            finalizeFinch()
            
            xy = Finch.finch_state.grid_data.allnodes;
            global maxerr;
            global L2Error = 0;

            for i=1:size(xy,2)
                x = xy[1,i];
                y = xy[2,i];
                err.values[i] = abs(u.values[i] - exact(x,y));
                exactval.values[i] = exact(x, y)
                global maxerr;
                maxerr = max(err.values[i], maxerr);
                global L2Error;
                L2Error += err.values[i]^2
            end
            println("max error = "*string(maxerr));
            println("L2 error = "*string(L2Error));

            # Uncomment below to plot or output result

            # using Plots
            # using PyPlot
            # # using Triplot
            # pyplot();
            # using DelimitedFiles

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # figure(2)
            # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

            # figure(1)
            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
            # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
            # display(h)

            # a = u.values[:]

            Nval =  match( r"[0-9]+", meshval )
            Nval = Nval.match
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "triangleMeshUnstruct_"
            file = open( textfoldername * prefix * "uvalues_lvl=" * Nval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_lvl=" * Nval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_lvl=" * Nval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_lvl=" * Nval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            outputValues( [u, err], textfoldername * prefix * "errorAndu_lvl=" * Nval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end
    end
end

function runFunc( k, c, meshval, pythonVarName )

    initFinch("mixed2d");

    useLog("mixed2dlog", level=3)

    domain(2)

    mesh(meshval)

    u = variable("u")
    err = variable("err")
    err = variable("err")
    testSymbol("v")
    exactval = variable("exact")
    # exact solution is sin(2*pi*x)*sin(2*pi*y)
    maxerr = 0;
    coeff = parse( Int64, k )
    coeff = parse( Int64, k )
    exact(x,y) = sin(pi*x*coeff)*sin(pi*y*coeff);

    linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

    boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")"
    # boundaryStr = "sin(pi*x)*sin(pi*y)"
    boundary(u, 1, DIRICHLET, boundaryStr)

    strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")"
    # strval = "pi*pi*sin(pi*x)*sin(pi*y)"
    coefficient("f", strval)

    weakForm(u, "dot(grad(u), grad(v)) - f*v")

    exportCode("mixed2dcode");
    importCode("mixed2dcodein");

    solve(u);

    finalizeFinch()
    
    xy = Finch.finch_state.grid_data.allnodes;
    global maxerr;
    global L2Error = 0;

    for i=1:size(xy,2)
        x = xy[1,i];
        y = xy[2,i];
        err.values[i] = abs(u.values[i] - exact(x,y));
        exactval.values[i] = exact(x, y)
        global maxerr;
        maxerr = max(err.values[i], maxerr);
        global L2Error;
        L2Error += err.values[i]^2
    end
    println("max error = "*string(maxerr));
    println("L2 error = "*string(L2Error));

    # Uncomment below to plot or output result

    # using Plots
    # using PyPlot
    # # using Triplot
    # pyplot();
    # using DelimitedFiles

    # figure(1)
    # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

    # figure(2)
    # display(plot(xy[1,:], xy[2,:], exactval.values[:], st=:contourf))

    # figure(1)
    # display(plot(xy[1,:], xy[2,:], u.values[:], st=:contour))
    # h = tricontourf(xy[1,:], xy[2,:], u.values[:]) 
    # display(h)

    # a = u.values[:]
    textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

    file = open( textfoldername * pythonVarName * "uvalues.txt", "w")
    print( file, u.values[:])
    close(file)

    textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

    file = open( textfoldername * pythonVarName * "uexactvalues.txt", "w")
    print( file, exactval.values[:])
    close(file)

    # print( u.values[:] )

    file = open( textfoldername * pythonVarName * "xvalues.txt", "w")
    print( file, xy[1, :])

    close(file)

    file = open( textfoldername * pythonVarName * "yvalues.txt", "w")
    print( file, xy[2, :])

    close(file)

    file = open( textfoldername * pythonVarName * "errorvalues.txt", "w")
    print( file, err.values[:])

    close(file)
    # for i = 1:size( u.values, 2 )
    #     println( file, u.values[i] )
    # end

    outputValues( [u, err], textfoldername * pythonVarName * "errorAndu", format="vtk" )

    # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
    # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

    # outputValues([u,err], "p2dmixed", format="vtk");
end

# mixedElement()


# for x in ARGS
#     println(x);
# end

k = ARGS[1]
c = ARGS[2]
pythonVarName = ARGS[4]
meshval = ARGS[3]

# regularElement(k, c)
# triangleElement(k, c)
# mixedElementWithLevel(k, c)

runFunc( k, c, meshval, pythonVarName )