#======================================================================
    Panel methods predict the pressure distribution around an airfoil
    under the assumptions of potential flow

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 21/08/2023
======================================================================#
using LinearAlgebra;


"""
PANEL1 calculates the cp distribution around a 2D airfoil and calculates
the main aerodynamic coefficients under the hypothesis of potential flow
(inviscid, incompressible, irrotational, steady)
Flat panels, linearly-varying vortex singularity
    (cp,CL,CM,CD,λ) = panel1(myairfoil,-5:0.5:15)
INPUT:
    airfoil: Nx2 matrix containing the surface points coordinates (TE->upper->LE->lower->TE)
    α: angle of attack in degrees
OUTPUT:
    cp: Nx2 matrix containing the pressure coefficients on the airfoil nodes
    CL: lift coefficient
    CM: pitching moment coefficient respect to the 25% of the chord
    CD: drag coefficient (must be close to zero)
    λ: vector containing the vortex singularity density magnitude on the airfoil nodes
"""
function panel1(airfoil, α)
    #rotate Vinf
    Vinf = zeros(2,length(α));
    for (i,AoA) in enumerate(α)
        Vinf[:,i] = [cos(AoA*pi/180), sin(AoA*pi/180)];
    end

    #locate the control points in the middle of each panel
    P = 0.5*(airfoil[2:end,:]+airfoil[1:end-1,:]);
    nP = size(P,1);     #number of control points

    #calculate the tangent and normal vector of each panel
    t = airfoil[2:end,:]-airfoil[1:end-1,:];
    li = sqrt.(t[:,1].^2+t[:,2].^2);
    t ./= li;
    n = similar(t);
    #n[:,1] = t[:,2];
    #n[:,2] = -t[:,1];
    n[:,1] = -t[:,2];
    n[:,2] = t[:,1];

    #define and populate the aerodynamic influence matrix
    A = zeros(nP+1,nP+1);
    b = zeros(nP+1,length(α));
    for i=1:nP
        for j=1:nP
            #calculating the effect that the j-th panel has on the i-th panel
            #l = norm(airfoil[j+1,:]-airfoil[j,:]);      #current panel length
            l = li[j];
            r = P[i,:]-P[j,:];                          #relative position between control points
            x = dot(r, t[j,:]);
            y = dot(r, n[j,:]);

            tmp1 = 0;
            tmp2 = -pi;
            if i!=j
                tmp1 = log((y^2+(x+0.5*l)^2)/(y^2+(x-0.5*l)^2));
                tmp2 = atan((x+0.5*l)/y)-atan((x-0.5*l)/y);
            end
            cax = -(1/(4*pi*l))*((l-2*x)*tmp2 + y*tmp1);
            cbx = -(1/(4*pi*l))*((l+2*x)*tmp2 - y*tmp1);
            cay = (1/(8*pi*l))*((l-2*x)*tmp1 - 4*y*tmp2 + 4*l);
            cby = (1/(8*pi*l))*((l+2*x)*tmp1 + 4*y*tmp2 - 4*l);

            A[i,j] += dot(cax*t[j,:]+cay*n[j,:], n[i,:]);
            A[i,j+1] += dot(cbx*t[j,:]+cby*n[j,:], n[i,:]);
        end
    end
    b[1:end-1,:] = n*Vinf;                              #equivalent to: dot(Vinf,n[i,:]) for each i
    A[nP+1,1] = 1;                                      #Kutta condition
    A[nP+1,nP+1] = 1;


    #solve linear system (all α at the same time) and calculates the pressure distribution
    λ = A\b;
    cp = zeros(nP+1,length(α)+1);
    cp[:,1] = airfoil[:,1];
    cp[:,2:end] = @. 1-λ^2;

    #calculate the main aerodynamic coefficients
    CL = zeros(length(α));
    CD = zeros(length(α));
    CM = zeros(length(α));
    for i=1:nP
        Fres = li[i]*(n[i,:]*transpose(3 .- λ[i+1,:].^2 .- λ[i,:].^2 .- λ[i+1,:].*λ[i,:])./3);      #resulting force on the i-th panel
        #Fres = (li[i]*n[i,:])*transpose(0.5.*(cp[i+1,:].+cp[i,:]));
        CL += Fres[2,:];
        CD += Fres[1,:];
        CM -= (P[i,1]-0.25)*Fres[2,:] - P[i,2]*Fres[1,:];       #equivalent to: -cross([P[i,1]-0.25, P[i,2], 0], [Fres[1,1], Fres[2,1], 0])
    end

    return (cp,CL,CM,CD,λ);
end
