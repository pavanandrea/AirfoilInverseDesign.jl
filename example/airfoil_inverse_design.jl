#======================================================================
    Example script that generates the airfoil geometry that generates
    a given target pressure distribution

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 08/09/2023
======================================================================#
using AirfoilInverseDesign;
using Plots;

#define the starting airfoil (a NACA-0009 with 100 nodes)
(airfoil0,airfoil0header) = generatenaca4airfoil("0009", 100);

#define the target pressure distribution from a set of 10 parameters
#and evaluate it at the airfoil0 nodes
params = [0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, 0.175];
cptarget = cpgen10h(params, airfoil0[:,1]);

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
