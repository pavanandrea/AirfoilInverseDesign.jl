#======================================================================
    Example script that performs a global optimization of the airfoil
    using the EMFID method and a 10-params flow-feature parametrization

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 08/09/2023
======================================================================#
using AirfoilInverseDesign;
using BlackBoxOptim;
using Plots;
using Xfoil;


#define the objective function
"""
OBJECTIVEFUN implements the objective function to minimize (e.g. drag on cruise)
    f = objectivefun(params)
INPUT:
    params: Vector of size 10 containing the cpgen10h parameters
OUTPUT:
    f: fitness value
"""
function objectivefun(params)
    #define operating conditions
    CLcruise = 0.40;
    CLclimb = 0.80;
    Re = 400_000;
    Ncrit = 9.0;

    #generate airfoil geometry using the MGM inverse design method
    (airfoil0,_) = generatenaca4airfoil("0009", 50);
    cptarget = cpgen10h(atanh.(params), airfoil0[:,1]);
    (airfoil,status) = mgm(cptarget, airfoil0);
    if !startswith(status, "COMPLETED")
        return NaN;
    end

    #analyze airfoil with Xfoil
    Xfoil.set_coordinates(airfoil[:,1],airfoil[:,2]);
    α = -2:0.5:15;              #range of angle of attacks in deg
    CL = ones(length(α));       #preallocating polar
    CD = ones(length(α));
    CM = ones(length(α));
    hasconverged = zeros(Bool,length(α));
    for i=1:lastindex(α)
        (CL[i],CD[i],_,CM[i],hasconverged[i]) = Xfoil.solve_alpha(α[i],Re,iter=100,ncrit=Ncrit,reinit=true);
    end
    if sum(hasconverged)<length(α)/2 || minimum(CD)<=0.0025
        #invalid values from Xfoil
        return NaN;
    end
    CL = CL[findall(hasconverged.==true)];
    CD = CD[findall(hasconverged.==true)];
    CM = CM[findall(hasconverged.==true)];
    CL = CL[findall(CM.!=1)];
    CD = CD[findall(CM.!=1)];
    CM = CM[findall(CM.!=1)];

    #calculate drag on operating conditions
    CDcruise = 0;
    CDclimb = 0;
    if minimum(CL)<=CLcruise<=maximum(CL) && minimum(CL)<=CLclimb<=maximum(CL)
        #interpolate CD values
        for i=1:lastindex(CL)-1
            if CL[i]<=CLcruise<=CL[i+1]
                CDcruise = CD[i] + (CLcruise-CL[i])*(CD[i+1]-CD[i])/(CL[i+1]-CL[i]);
            end
            if CL[i]<=CLclimb<=CL[i+1]
                CDclimb = CD[i] + (CLclimb-CL[i])*(CD[i+1]-CD[i])/(CL[i+1]-CL[i]);
            end
        end
    else
        #unable to interpolate CD values
        return NaN;
    end

    #calculate penalty on the CLmax
    CLmaxtarget = 1.30;         #minimum allowable value for maximum(CL)
    penalty1 = 20*max(0, CLmaxtarget-maximum(CL))^2;

    #calculate penalty on the thickness
    thicknesstarget = 0.12;     #minimum allowable value for the airfoil thickness
    penalty2 = 200*max(0, thicknesstarget-maximumthickness(airfoil))^2;

    #calculate penalty on the trailing-edge angle
    TEangletarget = 4.0;
    TEangle = abs((airfoil[end-2,2]-airfoil[end-4,2])/(airfoil[end-2,1]-airfoil[end-4,1]) - (airfoil[2,2]-airfoil[4,2])/(airfoil[2,1]-airfoil[4,1]))*180/pi;
    penalty3 = 2*max(0, TEangletarget-TEangle)^2;

    #calculate objective function
    return 0.55*CDcruise + 0.45*CDclimb + penalty1 + penalty2 + penalty3;
end


#check if multithreading is enabled
if Threads.nthreads() < 2
    println("WARNING: Julia is running on a single thread, the execution will be very slow. Please use the following command:");
    println("  julia airfoil_blackboxoptim.jl --threads 8");
end


#define a callback function, called by the optimizer after each
#optimization step in order to record and display the progress
history_nevals = Vector{Float64}(undef,0);
history_objectivefun = Vector{Float64}(undef,0);
function progresscallback(currentresult)
    push!(history_nevals, BlackBoxOptim.num_func_evals(currentresult));
    push!(history_objectivefun, best_fitness(currentresult));
    if length(history_nevals)%10 == 0
        println("  Number of evaluations: ",BlackBoxOptim.num_func_evals(currentresult),";\tobjective function: ",best_fitness(currentresult));
    end
end


#run black-box optimization
#see: https://en.wikipedia.org/wiki/Differential_evolution
println("Setting up optimizer");
optimizer = bbsetup(objectivefun;
    #params[1] = x-coordinate of the minimum cp (upper surface)
    #params[2] = minimum cp value (upper surface)
    #params[3] = cp(x) function curvature before the minimum cp (upper surface)
    #params[4] = cp(x) function curvature after the minimum cp (upper surface)
    #params[5] = cp value at x/c=0.9 (upper surface)
    #params[6] = x-coordinate of the minimum cp (lower surface)
    #params[7] = minimum cp value (lower surface)
    #params[8] = cp(x) function curvature before the minimum cp (lower surface)
    #params[9] = cp(x) function curvature after the minimum cp (lower surface)
    #params[10] = cp value at x/c=0.9 (lower surface)
    SearchRange = [                         #upper and lower bounds
        (tanh(0.22), tanh(0.35)),
        (tanh(-1.0), tanh(-0.6)),
        (tanh(0.20), tanh(0.60)),
        (tanh(0.00), tanh(4.00)),
        (tanh(-0.2), tanh(+0.2)),
        (tanh(0.02), tanh(0.20)),
        (tanh(-0.4), tanh(+0.1)),
        (tanh(1.50), tanh(5.00)),
        (tanh(3.00), tanh(7.00)),
        (tanh(-0.1), tanh(+0.3))
    ],
    Method = :adaptive_de_rand_1_bin_radiuslimited,     #choose optimizer from the following list: https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/examples/benchmarking/latest_toplist.csv
    MaxSteps = 500,                                     #maximum number of optimization steps
    TraceMode = :silent,                                #do not show output
    CallbackFunction = progresscallback,
    CallbackInterval = 0.0,
    NThreads = Threads.nthreads()-1
);
startingparams = tanh.([0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, 0.175]);      #initial candidate (NACA-4412)
println("Running black-box optimization...");
result = bboptimize(optimizer, startingparams);
optimparams = atanh.(best_candidate(result));
println("Optimization process completed");
println("  Best candidate: params=",optimparams);
println("  Objective function: ",best_fitness(result));


#plot optimization progress history
plt1 = plot(history_nevals, history_objectivefun,
    title = "Airfoil optimization progress",
    xlabel = "Number of function evaluations",
    ylabel = "Best objective function",
    yaxis = :log
);
display(plt1);


#plot resulting airfoil
(airfoil0,_) = generatenaca4airfoil("0009", 50);
cptarget = cpgen10h(optimparams, airfoil0[:,1]);
(airfoil,_) = mgm(cptarget, airfoil0);
plt2 = plot(airfoil[:,1], airfoil[:,2], label="Optimized",
    title = "Optimized airfoil",
    xlabel = "x/c",
    ylabel = "y/c",
    aspect_ratio = :equal
);
display(plt2);
