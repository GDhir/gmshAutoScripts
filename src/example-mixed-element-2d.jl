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

# meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=3Ny=9.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=7Ny=17.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=15Ny=33.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=31Ny=65.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=63Ny=129.msh" ]

# meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=8Ny=29.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=10Ny=35.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=12Ny=41.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=14Ny=47.msh",
#     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=16Ny=53.msh" ]



# meshvals = [ "/home/gaurav/Finch/src/examples/importMesh.msh", "/home/gaurav/Finch/src/examples/newMeshModifiedAllCorners.msh" ]
# meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv5Nx=3Ny=11.msh"]
# meshvals = [ "/home/gaurav/Finch/src/examples/mixedmesh.mesh" ]
# meshval = "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv433.msh"

# function regularElement()

    # meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=7.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=11.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=19.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=35.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=67.msh" ]

    # meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=41.msh",
    # "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=49.msh",
    # "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=57.msh",
    # "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=65.msh",
    # "/home/gaurav/Finch/src/examples/Mesh/regularMeshN=73.msh" ]

    meshvals = readdir( pwd()*"/Mesh/MeshRun", join = true )

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
            exact(x,y) = sin(pi*x*2)*sin(pi*y*2);

            boundary(u, 1, DIRICHLET, "sin(pi*x*2)*sin(pi*y*2)")

            coefficient("f", "-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)")

            weakForm(u, "dot(grad(u), grad(v)) + f*v")

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

            Nval =  match( r"[0-9]+", meshval )
            Nval = Nval.match
            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            prefix = "regular_"
            file = open( textfoldername * prefix * "uvalues_N=" * Nval * ".txt", "w")
            print( file, u.values[:])
            close(file)

            textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

            file = open( textfoldername * prefix * "uexactvalues_N=" * Nval * ".txt", "w")
            print( file, exactval.values[:])
            close(file)

            # print( u.values[:] )

            file = open( textfoldername * prefix * "xvalues_N=" * Nval * ".txt", "w")
            print( file, xy[1, :])

            close(file)

            file = open( textfoldername * prefix * "yvalues_N=" * Nval * ".txt", "w")
            print( file, xy[2, :])

            close(file)

            file = open( textfoldername * prefix * "errorvalues_N=" * Nval * ".txt", "w")
            print( file, err.values[:])

            close(file)
            # for i = 1:size( u.values, 2 )
            #     println( file, u.values[i] )
            # end

            outputValues( [u, err], textfoldername * prefix * "errorAnduN=" * Nval, format="vtk" )

            # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
            # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

            # outputValues([u,err], "p2dmixed", format="vtk");
        end

    end
# end

# function mixedElement()

    # meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=3Ny=9.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=7Ny=17.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=15Ny=33.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=31Ny=65.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv4Nx=63Ny=129.msh" ]

    # meshvals = [ "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=8Ny=29.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=10Ny=35.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=12Ny=41.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=14Ny=47.msh",
    #     "/home/gaurav/Finch/src/examples/Mesh/hangingMeshv7Nx=16Ny=53.msh" ]

    meshvals = readdir( pwd()*"/Mesh/MeshRun", join = true )

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
            exact(x,y) = sin(pi*x*2)*sin(pi*y*2);

            boundary(u, 1, DIRICHLET, "sin(pi*x*2)*sin(pi*y*2)")

            coefficient("f", "-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)")

            weakForm(u, "dot(grad(u), grad(v)) + f*v")

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
# end


meshvals = readdir( pwd()*"/Mesh/MeshRun", join = true )

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
        # exact solution is sin(2*pi*x)*sin(2*pi*y)
        maxerr = 0;
        exact(x,y) = sin(pi*x*2)*sin(pi*y*2);

        boundary(u, 1, DIRICHLET, "sin(pi*x*2)*sin(pi*y*2)")

        coefficient("f", "-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)")

        weakForm(u, "dot(grad(u), grad(v)) + f*v")

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

        Nval =  match( r"[0-9]+", meshval )
        Nval = Nval.match
        textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

        prefix = "triangleMeshStruct_"
        file = open( textfoldername * prefix * "uvalues_N=" * Nval * ".txt", "w")
        print( file, u.values[:])
        close(file)

        textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

        file = open( textfoldername * prefix * "uexactvalues_N=" * Nval * ".txt", "w")
        print( file, exactval.values[:])
        close(file)

        # print( u.values[:] )

        file = open( textfoldername * prefix * "xvalues_N=" * Nval * ".txt", "w")
        print( file, xy[1, :])

        close(file)

        file = open( textfoldername * prefix * "yvalues_N=" * Nval * ".txt", "w")
        print( file, xy[2, :])

        close(file)

        file = open( textfoldername * prefix * "errorvalues_N=" * Nval * ".txt", "w")
        print( file, err.values[:])

        close(file)
        # for i = 1:size( u.values, 2 )
        #     println( file, u.values[i] )
        # end

        outputValues( [u, err], textfoldername * prefix * "errorAnduN=" * Nval, format="vtk" )

        # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
        # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

        # outputValues([u,err], "p2dmixed", format="vtk");
    end
end

meshvals = readdir( pwd()*"/Mesh/MeshRun", join = true )

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
        exact(x,y) = sin(pi*x*2)*sin(pi*y*2);

        boundary(u, 1, DIRICHLET, "sin(pi*x*2)*sin(pi*y*2)")

        coefficient("f", "-8*pi*pi*sin(2*pi*x)*sin(2*pi*y)")

        weakForm(u, "dot(grad(u), grad(v)) + f*v")

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

        Nval =  match( r"[0-9]+", meshval )
        Nval = Nval.match
        textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

        prefix = "triangleMeshUnstruct_"
        file = open( textfoldername * prefix * "uvalues_N=" * Nval * ".txt", "w")
        print( file, u.values[:])
        close(file)

        textfoldername = "/media/gaurav/easystore/Finch/MixedElement/TextFiles/"

        file = open( textfoldername * prefix * "uexactvalues_N=" * Nval * ".txt", "w")
        print( file, exactval.values[:])
        close(file)

        # print( u.values[:] )

        file = open( textfoldername * prefix * "xvalues_N=" * Nval * ".txt", "w")
        print( file, xy[1, :])

        close(file)

        file = open( textfoldername * prefix * "yvalues_N=" * Nval * ".txt", "w")
        print( file, xy[2, :])

        close(file)

        file = open( textfoldername * prefix * "errorvalues_N=" * Nval * ".txt", "w")
        print( file, err.values[:])

        close(file)
        # for i = 1:size( u.values, 2 )
        #     println( file, u.values[i] )
        # end

        outputValues( [u, err], textfoldername * prefix * "errorAnduN=" * Nval, format="vtk" )

        # display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
        # display(plot(xy[1,:], xy[2,:], err.values[:], st=:surface))

        # outputValues([u,err], "p2dmixed", format="vtk");
    end
end
# mixedElement()
# regularElement()