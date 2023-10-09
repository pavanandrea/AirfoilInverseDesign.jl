#======================================================================
    Testing script for mgm.jl

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 22/08/2023
======================================================================#
using Plots;
include("../src/airfoilutils.jl");
include("../src/mgm.jl");


#define the target airfoil and analyze it with panel1
(airfoiltarget,airfoiltargetheader) = generatenaca4airfoil("4412",100);
(cptarget,_,_,_,_) = panel1(airfoiltarget,0);

#define the starting airfoil and call the MGM method
#using the pressure distribution from the target airfoil
#hopefully, the MGM method returns an airfoil similar to the target
println("Airfoil inverse design...");
(airfoil0,airfoil0header) = generatenaca4airfoil("0009",100);
(airfoil,status) = @time mgm(cptarget, airfoil0);
println("Completed with status '",status,"'");


#compare generated airfoil with target
plt1 = plot(airfoil[:,1], airfoil[:,2], label="MGM method",
    title = "Inverse design airfoil - geometry comparison",
    xlabel = "x/c",
    ylabel = "y/c",
    #xlims = [-0.1,1.1],
    #ylims = [-0.6,0.6]
);
plot!(plt1, airfoil0[:,1], airfoil0[:,2], label=airfoil0header);
plot!(plt1, airfoiltarget[:,1], airfoiltarget[:,2], label=airfoiltargetheader);
display(plt1);


#compare pressure distributions
(cp0,_,_,_,_) = panel1(airfoil0,0);
(cp,_,_,_,_) = panel1(airfoil,0);
plt2 = plot(cp[:,1], cp[:,2:end], label="MGM method",
    title = "Inverse design airfoil - pressure comparison",
    xlabel = "x/c",
    ylabel = "cp",
    yflip = true
);
plot!(plt2, cp0[:,1], cp0[:,2:end], label=airfoil0header);
plot!(plt2, cptarget[:,1], cptarget[:,2:end], label=airfoiltargetheader);
display(plt2);
