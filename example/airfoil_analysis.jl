#======================================================================
    Example script that generates a NACA-4412 airfoil and calculates
    the polar curves for various Reynolds numbers using Xfoil.

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 07/09/2023
======================================================================#
using AirfoilInverseDesign;
using Plots;
using Xfoil;


#generate NACA airfoil
(airfoil1,airfoil1header) = generatenaca4airfoil("4412",100);
Xfoil.set_coordinates(airfoil1[:,1],airfoil1[:,2]);


#generate airfoil polars
α = -5:0.5:15;                                  #range of angle of attacks in deg
Re = [100_000 200_000 500_000 1_000_000];       #list of Reynolds numbers
CL = ones(length(α),length(Re),);               #preallocating results
CD = ones(length(α),length(Re));
CM = ones(length(α),length(Re));
hasconverged = zeros(Bool,(length(α),length(Re)));
println("Analyzing ",airfoil1header,"...");
@time for i=1:lastindex(Re)
    for j=1:lastindex(α)
        (CL[j,i],CD[j,i],_,CM[j,i],hasconverged[j,i]) = Xfoil.solve_alpha(α[j],Re[i],iter=100,reinit=true);
        if !hasconverged[j,i]
            println("Xfoil is unable to converge at Re=",Re[i],", α=",α[j],"°");
        end
    end
end


#plot CL-α
plt1 = plot(α, CL,
    title = airfoil1header*" lift curves",
    label = "Re = ".*string.(Re),
    xlabel = "Angle of attack α [deg]",
    ylabel = "Lift coefficient CL"
);
display(plt1);

#plot CM-α
plt2 = plot(α, CM,
    title = airfoil1header*" moment curves",
    label = "Re = ".*string.(Re),
    xlabel = "Angle of attack α [deg]",
    ylabel = "Moment coefficient CM"
);
display(plt2);

#plot CL-CD
plt3 = plot(CD, CL,
    title = airfoil1header*" polars",
    label = "Re = ".*string.(Re),
    xlabel = "Drag coefficient CD",
    ylabel = "Lift coefficient CL"
);
display(plt3);

#plot CL-E
plt4 = plot(CL, CL./CD,
    title = airfoil1header*" efficiency",
    label = "Re = ".*string.(Re),
    xlabel = "Lift coefficient CL",
    ylabel = "Aerodynamic efficiency CL/CD"
);
display(plt4);
