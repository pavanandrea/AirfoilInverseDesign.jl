#======================================================================
    Example script that performs a local optimization of the airfoil
    using the EMFID method and a 10-params flow-feature parametrization

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 08/09/2023
======================================================================#
using AirfoilInverseDesign;
using Optim;
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


#=
#Particle Swarm optimization
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
optimizer = ParticleSwarm(;
    lower = tanh.([0.22, -1.0, 0.20, 0.00, -0.2, 0.02, -0.4, 1.50, 3.00, -0.1]),        #params lower bound
    upper = tanh.([0.35, -0.6, 0.60, 4.00, +0.2, 0.20, +0.1, 5.00, 7.00, +0.3]),        #params upper bound
    n_particles = 200
);
options = Optim.Options(
    iterations = 100,
    store_trace = true,
    show_trace = true,
    show_every = 1
);
startingparams = tanh.([0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, 0.175]);   #initial candidate (NACA-4412)
println("Running optimization...");
result = optimize(objectivefun, startingparams, optimizer, options);
println("Optimization process completed");
println("  Best candidate: params=",atanh.(Optim.minimizer(result)));
println("  Converged: ",Optim.converged(result));
println("  Objective function: ",Optim.minimum(result));
=#


#Nelder-Mead local optimization
options = Optim.Options(
    #iterations = 30,
    store_trace = true,
    show_trace = true,
    show_every = 10
);
startingparams = tanh.([0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, 0.175]);   #initial candidate (NACA-4412)
println("Running optimization...");
lb = tanh.([0.22, -1.0, 0.20, 0.00, -0.2, 0.02, -0.4, 1.50, 3.00, -0.1]);        #params lower bound
ub = tanh.([0.35, -0.6, 0.60, 4.00, +0.2, 0.20, +0.1, 5.00, 7.00, +0.3]);        #params upper bound
result = optimize(objectivefun, lb, ub, startingparams, Fminbox(NelderMead()), options);
optimparams = atanh.(Optim.minimizer(result));
println("Optimization process completed");
println("  Best candidate: params=",optimparams);
println("  Converged: ",Optim.converged(result));
println("  Objective function: ",Optim.minimum(result));


#plot optimization progress history
plt1 = plot(1:length(Optim.f_trace(result)), Optim.f_trace(result),
    title = "Airfoil optimization progress",
    legend = false,
    xlabel = "Number of function evaluations",
    ylabel = "Best objective function",
    ylims = [minimum(Optim.f_trace(result))-0.0005, 0.01]
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


#=
#airfoil geometry comparison (NACA-4412 and optimized)
(airfoil1,airfoil1header) = generatenaca4airfoil("4412",50);
plt3 = plot(airfoil1[:,1], airfoil1[:,2], label="NACA-4412",
    title = "Airfoil comparison",
    xlabel = "x/c",
    ylabel = "y/c",
    aspect_ratio = :equal
);
plot!(plt3, airfoil[:,1], airfoil[:,2], label="Optimized");
display(plt3);

#airfoil polar comparison (NACA-4412 and optimized)
Xfoil.set_coordinates(airfoil1[:,1],airfoil1[:,2]);
α = -2:0.5:15;
CL1 = ones(length(α));
CD1 = ones(length(α));
hasconverged = zeros(Bool,length(α));
for i=1:lastindex(α)
    (CL1[i],CD1[i],_,_,hasconverged[i]) = Xfoil.solve_alpha(α[i],400_000,iter=100,ncrit=9,reinit=true);
end
CL1 = CL1[findall(hasconverged.==true)];
CD1 = CD1[findall(hasconverged.==true)];
Xfoil.set_coordinates(airfoil[:,1],airfoil[:,2]);
CL = ones(length(α));
CD = ones(length(α));
hasconverged = zeros(Bool,length(α));
for i=1:lastindex(α)
    (CL[i],CD[i],_,_,hasconverged[i]) = Xfoil.solve_alpha(α[i],400_000,iter=100,ncrit=9,reinit=true);
end
CL = CL[findall(hasconverged.==true)];
CD = CD[findall(hasconverged.==true)];
#plot CL-CD
plt4 = plot(CD1, CL1, label="NACA-4412",
    title = "Airfoils polars comparison",
    xlabel = "Drag coefficient CD",
    ylabel = "Lift coefficient CL"
);
plot!(plt4, CD, CL, label="Optimized");
display(plt4);
=#
