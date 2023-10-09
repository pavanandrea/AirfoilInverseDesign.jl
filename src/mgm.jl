#======================================================================
    The Modified Garabedian-McFadden method generates the airfoil that
    matches a given pressure distribution (inverse design)

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 22/08/2023
======================================================================#
using LinearAlgebra;
#using Roots;
include("panel.jl");


"""
MGM uses an elastic-membrane approach to gradually change the shape of
an airfoil such that the resulting pressure distribution matches the
target as closely as possible
    (airfoil) = mgm(cptarget, airfoil0)
INPUT:
    cptarget: Nx2 containing the target pressure distribution
    airfoil0: Nx2 matrix containing the surface points coordinates of the starting airfoil
OPTIONAL INPUTS:
    itmax: maximum number of iterations (default: 100)
    toll: desired tolerance on the pressure distribution (default: 1e-4)
    constants: vector of size 3 contining the weights of the equation (default: [14,0,6])
OUTPUT:
    airfoil: Nx2 matrix containing the surface points coordinates of the resulting airfoil
    status: string with the status message (eg: OK, ERROR, etc...)
"""
function mgm(cptarget, airfoil0, itmax=100, toll=1e-2, constants=[14,0,6])
    status = "ONGOING";

    #find the leading edge
    airfoil = copy(airfoil0);
    t = airfoil[2:end,:]-airfoil[1:end-1,:];
    li = sqrt.(t[:,1].^2+t[:,2].^2);
    t ./= li;
    #LEidx = findfirst(t[:,1].>0);
    LEidx = findfirst(abs.(airfoil[:,1]).==minimum(abs.(airfoil[:,1])));

    #transform global x coordinates to local s curvilinear coordinates
    s = 0*cptarget[:,1];
    for i=2:lastindex(s)
        s[i] = s[i-1] + li[i-1];
    end

    #elastic-membrane iterations
    (A,B,C) = constants;
    (cpi,_,_,_,_) = panel1(airfoil,0);
    for i=1:itmax
        #=
        #calculate the Δy steps
        a = @. 2*C/((s[3:end]-s[1:end-2])*(s[2:end-1]-s[1:end-2]));
        b = @. A - B/(s[3:end]-s[2:end-1]) - 2*C/((s[3:end]-s[1:end-2])*(s[3:end]-s[2:end-1])) - 2*C/((s[3:end]-s[1:end-2])*(s[2:end-1]-s[1:end-2]));
        c = @. B/(s[3:end]-s[2:end-1]) + 2*C/((s[3:end]-s[1:end-2])*(s[3:end]-s[2:end-1]));
        d = cptarget[2:end-1,2] - cpi[2:end-1,2];
        =#

        #calculate the Δy steps for the upper surface
        a = @. 2*C/((s[3:LEidx]-s[1:LEidx-2])*(s[2:LEidx-1]-s[1:LEidx-2]));
        b = @. A - B/(s[3:LEidx]-s[2:LEidx-1]) - 2*C/((s[3:LEidx]-s[1:LEidx-2])*(s[3:LEidx]-s[2:LEidx-1])) - 2*C/((s[3:LEidx]-s[1:LEidx-2])*(s[2:LEidx-1]-s[1:LEidx-2]));
        c = @. B/(s[3:LEidx]-s[2:LEidx-1]) + 2*C/((s[3:LEidx]-s[1:LEidx-2])*(s[3:LEidx]-s[2:LEidx-1]));
        d = cptarget[2:LEidx-1,2] - cpi[2:LEidx-1,2];
        M = Tridiagonal(a[2:end],b,c[1:end-1]);
        Δyupper = M\d;
        
        #calculate the Δy steps for the lower surface
        a = @. 2*C/((s[LEidx+2:end]-s[LEidx:end-2])*(s[LEidx+1:end-1]-s[LEidx:end-2]));
        b = @. A - B/(s[LEidx+2:end]-s[LEidx+1:end-1]) - 2*C/((s[LEidx+2:end]-s[LEidx:end-2])*(s[LEidx+2:end]-s[LEidx+1:end-1])) - 2*C/((s[LEidx+2:end]-s[LEidx:end-2])*(s[LEidx+1:end-1]-s[LEidx:end-2]));
        c = @. B/(s[LEidx+2:end]-s[LEidx+1:end-1]) + 2*C/((s[LEidx+2:end]-s[LEidx:end-2])*(s[LEidx+2:end]-s[LEidx+1:end-1]));
        d = cptarget[LEidx+1:end-1,2] - cpi[LEidx+1:end-1,2];
        M = Tridiagonal(a[2:end],b,c[1:end-1]);
        Δylower = M\d;

        #update airfoil coordinates
        airfoil[2:LEidx-1,2] += Δyupper;
        airfoil[LEidx+1:end-1,2] -= Δylower;
        if maximum(abs.(airfoil[:,2]))>1.0
            status = "ERROR while updating airfoil: unable to reach convergence";
            break;
        end

        #check if the pressure distribution is close enough to the target
        (cpi,_,_,_,_) = panel1(airfoil,0);
        cprms = sqrt(sum((cptarget[:,2]-cpi[:,2]).^2)/length(cpi[:,2]));
        if cprms<=toll
            status = "COMPLETED after reaching convergence in "*string(i)*" iterations";
            break;
        end
    end

    if !startswith(status,"ERROR") && !startswith(status,"COMPLETED")
        status = "COMPLETED after reaching the maximum number of iterations";
    end
    return (airfoil,status);
end


"""
CPGEN10H generates a pressure distribution suitable for incompressible low-Re airfoils
    cp = cpgen10h(params,xcp)
INPUT:
    params: Vector of size 10 containing the parameters
OPTIONAL INPUT:
    xcp: Vector of size N containing the x-coordinates of the generated cp distribution
OUTPUT:
    cp: Nx2 matrix containing the generated pressure distribution
"""
function cpgen10h(params, xcp=vcat(reverse([0.003:0.003:0.015,0.02:0.02:0.98,0.995,1]),[0:0.003:0.015,0.02:0.02:0.98,0.995,1]))
    if length(params)!=10
        error("ERROR in cpgen10h: unsupported number of parameters (must be 10, but ",length(params)," are given)");
    end

    #find the leading edge
    LEidx = findfirst(abs.(xcp).==minimum(abs.(xcp)));

    #parameters description:
    #params[1] = x-coordinate of the minimum cp (upper surface)
    #params[2] = minimum cp value (upper surface)
    #params[3] = cp(x) function curvature before the minimum cp (upper surface)
    #params[4] = cp(x) function curvature after the minimum cp (upper surface)
    #params[5] = cp value at x/c=0.9 (upper surface)
    #params[6] = x-coordinate of the minimum cp (lower surface)
    #params[7] = minimum cp value (lower surface)
    #params[8] = cp(x) function curvature before the minimum cp (lower surface)
    #params[9] = cp(x) function curvature after the minimum cp (lower surface)
    #params[10] = cp value at x/c=0.9 (lower surface)

    #define the parametrized cp function
    hickshenne(x,t1,t2,t3) = @. (t2-1)*sin(pi*x^(log(0.5)/log(t1)))^t3;
    function parametrizedcp(x,t)
        if x<=0.0075
            #straight line between 0 and 0.0075
            return @. 1 + x * hickshenne(0.0075,t[1],t[2],t[3])/0.0075;
        elseif 0.0075<x<=t[1]
            #Hicks-Henne curve between 0.0075 and cpmin
            return @. 1 + hickshenne(x,t[1],t[2],t[3]);
        elseif t[1]<x<=t[1]+0.15
            #Hicks-Henne curve between cpmin and params[1]+0.15
            return @. hickshenne(x,t[1],t[2]+1,t[4]);
        elseif t[1]+0.15<x<=0.9
            #cubic spline between params[1]+0.15 and 0.9
            cp4 = hickshenne(t[1]+0.15,t[1],t[2]+1,t[4]);
            dcp4 = (cp4-hickshenne(t[1]+0.15-1e-6,t[1],t[2]+1,t[4]))/1e-6;
            #5 = find_zero(Ai -> t[5]-1-Ai*0.740969719^(-0.4144031946148317*Ai - 0.22596643509454467), t[5]);
            A5 = -1;
            ΔA5 = 1;
            for i=1:1000
                #find the value of the constant 2A using a fixed-point iteration method
                ΔA5 = (t[5]-1)/(0.740969719^(-0.4144031946148317*A5 - 0.22596643509454467)) - A5;
                A5 += ΔA5;
                if abs(ΔA5)<1e-9
                    break;
                end
            end
            C5 = -0.4144031946148317*A5 - 0.22596643509454467;
            dcp5 = (t[5]-(1+A5*sin(pi*(1-0.9+1e-6)^(log(0.5)/log(0.3)))^C5))/1e-6;
            a = (dcp4+dcp5-2*(t[5]-cp4)/(0.9-(t[1]+0.15)))/(((t[1]+0.15)-0.9)^2);
            b = 0.5*(dcp5-dcp4)/(0.9-(t[1]+0.15)) - 1.5*(0.9+(t[1]+0.15))*a;
            c = dcp4 - 3*a*(t[1]+0.15)^2 - 2*b*(t[1]+0.15);
            d = cp4 - a*(t[1]+0.15)^3 - b*(t[1]+0.15)^2 - c*(t[1]+0.15);
            return @. a*x^3 + b*x^2 + c*x + d;
        else
            #Hicks-Henne curve after 0.9
            #A = find_zero(Ai -> t[5]-1-Ai*0.740969719^(-0.4144031946148317*Ai - 0.22596643509454467), t[5]);    #from the equation cp(x=0.9)=params[5]
            A = -1;
            ΔA = 1;
            for i=1:1000
                #find the value of the constant A using a fixed-point iteration method
                ΔA = (t[5]-1)/(0.740969719^(-0.4144031946148317*A - 0.22596643509454467)) - A;
                A += ΔA;
                if abs(ΔA)<1e-9
                    break;
                end
            end
            C = -0.4144031946148317*A - 0.22596643509454467;                                                    #from the trailing edge coefficients regression
            return @. 1 + A*sin(pi*(1-x)^(log(0.5)/log(0.3)))^C;
        end
        #=elseif 0.9<x<=0.995
            A = find_zero(Ai -> t[5]-1-Ai*0.740969719^(-0.4144031946148317*Ai - 0.22596643509454467), t[5]);    #from the equation cp(x=0.9)=params[5]
            C = -0.4144031946148317*A - 0.22596643509454467;                                                    #from the trailing edge coefficients regression
            return @. 1 + A*sin(pi*(1-x)^(log(0.5)/log(0.3)))^C;
        else
            A = find_zero(Ai -> t[5]-1-Ai*0.740969719^(-0.4144031946148317*Ai - 0.22596643509454467), t[5]);    #from the equation cp(x=0.9)=params[5]
            C = -0.4144031946148317*A - 0.22596643509454467;                                                    #from the trailing edge coefficients regression
            cpend = @. 1 + A*sin(pi*(1-x)^(log(0.5)/log(0.3)))^C;
            return @. cpend + (x-0.995)*(1-cpend)/0.005;
        end=#
    end

    #generate pressure distribution
    cp = zeros(length(xcp),2);
    for i=1:lastindex(xcp)
        cp[i,1] = xcp[i];
        if i<=LEidx
            #upper surface
            cp[i,2] = parametrizedcp(xcp[i],params[1:5]);
        else
            #lower surface
            cp[i,2] = parametrizedcp(xcp[i],params[6:10]);
        end
    end
    return cp;
end
