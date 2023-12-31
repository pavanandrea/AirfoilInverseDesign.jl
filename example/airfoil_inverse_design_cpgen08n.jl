#======================================================================
    Example script that given a target pressure distribution as input
    generates the corresponding airfoil geometry

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 30/12/2023
======================================================================#
using AirfoilInverseDesign;
using Plots;

#define the starting airfoil (a NACA-0009 with 20 nodes)
(airfoil0,airfoil0header) = generatenaca4airfoil("0009", 20);

#define the target pressure distribution from a set of 8 parameters
params = [0.037, 0.0, 0.0, 0.0215, -0.226, -0.02, -0.015, 0.0812];
cptarget = cpgen08n(params);

#airfoil inverse design
(airfoil,status) = mgm(cptarget, airfoil0);

#analyze the generated airfoil with an inviscid panel method at α=0°
(cp,CL,CM,CD,_) = panel1(airfoil, 0);

#compare pressure distributions
plt1 = plot(cptarget[:,1], cptarget[:,2], label="Target",
    title = "Pressure distribution comparison",
    xlabel = "x/c",
    ylabel = "cp",
    yflip = true
);
plot!(plt1, cp[:,1], cp[:,2], label="Generated");
display(plt1);

#plot the new airfoil geometry
plt2 = plot(airfoil0[:,1], airfoil0[:,2], label="Starting",
    title = "Airfoil comparison",
    xlabel = "x/c",
    ylabel = "y/c",
    aspect_ratio = :equal
);
plot!(plt2, airfoil[:,1], airfoil[:,2], label="Generated");
display(plt2);
