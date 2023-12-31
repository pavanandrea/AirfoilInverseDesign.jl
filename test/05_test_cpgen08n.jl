#======================================================================
    Testing script for mgm.jl

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 30/12/2023
======================================================================#
using Plots;
include("../src/airfoilutils.jl");
include("../src/mgm.jl");


#define a pressure distribution using the cpgen08n function
#and compare to the one provided by the panel method
#to assess the representativeness of the parametrization
(airfoil,airfoilheader) = generatenaca4airfoil("4412",20);
(cpexact,_,_,_,_) = panel1(airfoil,0);
#params = 2*(rand(8).-0.5);
params = [0.037, 0.0, 0.0, 0.0215, -0.226, -0.02, -0.015, 0.0812];
println("Generating pressure distribution...");
cpgenerated = @time cpgen08n(params,airfoil[:,1]);
println("Completed");


#compare the pressure distributions
plt1 = plot(cpgenerated[:,1], cpgenerated[:,2:end], label="cpgen08n",
    title = "Pressure distribution comparison",
    xlabel = "x/c",
    ylabel = "cp",
    yflip = true
);
plot!(plt1, cpexact[:,1], cpexact[:,2:end], label=airfoilheader);
display(plt1);


#=
#manually find the params that best approximate the 4412 pressure distribution
function loss(p)
    return sum((cpgen08n(p,airfoil[:,1])-cpexact).^2);
end
prange = -0.02:0.01:+0.02;
minL = 10.0;
minp = [0, 0, 0, 0, 0, 0, 0, 0];
for i1 in 0.03:0.001:0.04
    for i2 in prange
        for i3 in prange
            for i4 in 0.020:0.0005:0.022
                for i5 in -0.23:0.001:-0.22
                    for i6 in prange
                        for i7 in -0.025:0.01:-0.015
                            for i8 in 0.079:0.0001:0.082
                                if loss([i1,i2,i3,i4,i5,i6,i7,i8])<minL
                                    global minp = [i1,i2,i3,i4,i5,i6,i7,i8];
                                    global minL = loss(minp);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
println("min loss: ",minL," - ",string(minp));
=#
