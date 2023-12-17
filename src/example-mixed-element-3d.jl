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


function runFunc( k, c, meshval, pythonVarName )

    initFinch("mixed2d");

    useLog("mixed2dlog", level=3)

    domain(3)

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
    exact(x,y,z) = sin(pi*x*coeff)*sin(pi*y*coeff)*sin(pi*z*coeff);

    linAlgOptions( matrixFree = false, iterative = false, reltol = 1e-20 )

    boundaryStr = "sin(pi*x*"*k*")*sin(pi*y*"*k*")*sin(pi*z*"*k*")"
    # boundaryStr = "sin(pi*x)*sin(pi*y)"
    boundary(u, 1, DIRICHLET, boundaryStr)

    strval = c*"*pi*pi*sin(pi*x*"*k*")*sin(pi*y*"*k*")*sin(pi*z*"*k*")"
    # strval = "pi*pi*sin(pi*x)*sin(pi*y)"
    coefficient("f", strval)

    weakForm(u, "dot(grad(u), grad(v)) - f*v")

    exportCode("mixed2dcode");
    # importCode("mixed2dcodein");

    solve(u);

    finalizeFinch()
    
    xyz = Finch.finch_state.grid_data.allnodes;
    global maxerr;
    global L2Error = 0;

    for i=1:size(xyz, 2)
        x = xyz[1,i];
        y = xyz[2,i];
        z = xyz[3,i];
        err.values[i] = abs(u.values[i] - exact(x,y,z));
        exactval.values[i] = exact(x, y, z)
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
    print( file, xyz[1, :])

    close(file)

    file = open( textfoldername * pythonVarName * "yvalues.txt", "w")
    print( file, xyz[2, :])

    close(file)

    file = open( textfoldername * pythonVarName * "zvalues.txt", "w")
    print( file, xyz[3, :])

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
meshval = ARGS[3]
pythonVarName = ARGS[4]

# regularElement(k, c)
# triangleElement(k, c)
# mixedElementWithLevel(k, c)

runFunc( k, c, meshval, pythonVarName )