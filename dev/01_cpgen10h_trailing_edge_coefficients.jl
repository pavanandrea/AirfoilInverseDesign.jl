#======================================================================
    Development script for mgm.jl
    
    While developing the cpgen10h parametrization, it can be observed
    for many airfoils and pressure distributions that the reverse
    Hicks-Henne curves on the trailing are showing a linear relation
    between their curvatures and their minimum values.
    This script makes a simple linear regression to find this relation.

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 22/08/2023
======================================================================#
using Plots;


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
println("Linear regression: C = ",coeffs[1]," + A * ",coeffs[2]);
plt2 = scatter(At, Ct,
    legend = false,
    xlabel = "A",
    ylabel = "C"
);
scatter!(plt2, Ab, Cb);
plot!(plt2, LinRange(-1.2,-0.8,100), coeffs[1].+coeffs[2].*LinRange(-1.2,-0.8,100));
display(plt2);
