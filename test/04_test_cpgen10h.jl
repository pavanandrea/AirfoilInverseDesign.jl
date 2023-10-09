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


#define a pressure distribution using the cp10param function
#and compare to the one provided by the panel method
#to assess the representativeness of the parametrization
(airfoil,airfoilheader) = generatenaca4airfoil("4412",100);
(cpexact,_,_,_,_) = panel1(airfoil,0);
#params = [0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 2.0, 5.5, +0.175];
params = [0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, +0.175];
println("Generating pressure distribution...");
cpgenerated = @time cpgen10h(params,airfoil[:,1]);
println("Completed");


#compare the pressure distributions
plt1 = plot(cpgenerated[:,1], cpgenerated[:,2:end], label="cp10param",
    title = "Pressure distribution comparison",
    xlabel = "x/c",
    ylabel = "cp",
    yflip = true
);
plot!(plt1, cpexact[:,1], cpexact[:,2:end], label=airfoilheader);
display(plt1);


#=
#trailing edge coefficients regression
#measured from: NACA-0012, NACA-2412, NACA-4412
At = [-1.16726954857861, -1.03370184099848, -1.09166373706037];
Ct = [0.263328716673221, 0.206097805328380, 0.215895587034236];
Ab = [-0.875061336169738, -0.991943600783843, -0.905244368324774];
Cb = [0.145754828082820, 0.188575360101011, 0.137856575768736];
X = vcat(At,Ab);
Y = vcat(Ct,Cb);
A = [6 sum(X); sum(X) sum(X.*X)];
b = [sum(Y), sum(X.*Y)];
coeffs = A\b;
println("coeffs = ");
display(coeffs);
plt2 = scatter(At, Ct,
    legend = false,
    xlabel = "A",
    ylabel = "C"
);
scatter!(plt2, Ab, Cb);
plot!(plt2, LinRange(-1.2,-0.8,100), coeffs[1].+coeffs[2].*LinRange(-1.2,-0.8,100));
display(plt2);
=#
