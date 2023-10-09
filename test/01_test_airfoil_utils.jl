#======================================================================
    Testing script for airfoilutils.jl

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 21/08/2023
======================================================================#
using Plots;
include("../src/airfoilutils.jl");


#generate NACA airfoil
println("Generating NACA 4 digits airfoil...");
(airfoil1,airfoil1header) = generatenaca4airfoil("4412",100);
println("Completed");
println("Airfoil '",airfoil1header,"':");
display(airfoil1);


#save airfoil as DAT file
println("Exporting airfoil to DAT file...");
#exportairfoildat(airfoil1, airfoil1header, "01_naca_4412_generated.dat");
exportairfoildat(airfoil1, airfoil1header, joinpath(@__DIR__,"01_naca_4412_generated.dat"));
println("Completed");


#import airfoil from DAT file
println("Importing DAT file...")
#(airfoil2,airfoil2header) = importairfoilfromfile("01_naca_4412_generated.dat");
(airfoil2,airfoil2header) = importairfoilfromfile(joinpath(@__DIR__,"01_naca_4412_generated.dat"));
println("Completed");
plt1 = plot(airfoil2[:,1], airfoil2[:,2],
    title = airfoil2header,
    legend = false,
    xlabel = "x/c",
    ylabel = "y/c",
    xlims = [-0.1,1.1],
    ylims = [-0.6,0.6]
);
display(plt1);


#compare with a NACA-4412 from airfoiltools.com
#(airfoiltest,airfoiltestheader) = importairfoilfromfile("01_test_airfoil_naca4412.dat");
(airfoilreference,airfoilreferenceheader) = importairfoilfromfile(joinpath(@__DIR__,"01_naca_4412_reference.dat"));
plt2 = plot(airfoil2[:,1], airfoil2[:,2], label="Generated",
    title = airfoil2header,
    xlabel = "x/c",
    ylabel = "y/c"
);
scatter!(plt2, airfoilreference[:,1], airfoilreference[:,2], label="Reference");
display(plt2);
