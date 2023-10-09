#======================================================================
    Testing script for the whole package

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 21/08/2023
======================================================================#
using AirfoilInverseDesign;
using Test;

#println("Testing AirfoilInverseDesign...");

#airfoilutils
(airfoil1,airfoil1header) = generatenaca4airfoil("4412",100);
@test sum(isnan.(airfoil1))==0 && size(airfoil1)==(100,2);
exportairfoildat(airfoil1, airfoil1header, joinpath(@__DIR__,"01_naca_4412_generated.dat"));
(airfoil2,airfoil2header) = importairfoilfromfile(joinpath(@__DIR__,"01_naca_4412_generated.dat"));
@test sum(isnan.(airfoil2))==0 && airfoil2==airfoil1;
@test 0.119<=maximumthickness(airfoil1)<=0.121;

#panel
(cptarget,_,_,_,_) = panel1(airfoil2,0);
@test sum(isnan.(cptarget))==0 && maximum(abs.(cptarget))<2;

#mgm
(airfoil3,status) = mgm(cptarget, airfoil2);
@test sum(isnan.(airfoil3))==0 && startswith(status,"COMPLETED") && airfoil3==airfoil2;
params = [0.27, -0.79, 0.4, 1.8, -0.1, 0.04, -0.275, 3.3, 5.5, +0.175];
cpgenerated = cpgen10h(params,airfoil3[:,1]);
@test sum(isnan.(cpgenerated))==0 && -0.80<=minimum(cpgenerated[:,2])<=-0.78
