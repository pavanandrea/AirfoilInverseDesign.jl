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
            cp[i,2] = parametrizedcp(xcp[i],atanh.(params[1:5]));
        else
            #lower surface
            cp[i,2] = parametrizedcp(xcp[i],atanh.(params[6:10]));
        end
    end
    return cp;
end


"""
CPGEN08N generates a pressure distribution suitable for incompressible low-Re airfoils
    cp = cpgen08n(params,xcp)
INPUT:
    params: Vector of size 8 containing the parameters
OPTIONAL INPUT:
    xcp: Vector of size N containing the x-coordinates of the generated cp distribution
OUTPUT:
    cp: Nx2 matrix containing the generated pressure distribution
"""
function cpgen08n(params, xcp=[1.0, 0.9698463103929542, 0.883022221559489, 0.75, 0.5868240888334652, 0.4131759111665348, 0.25, 0.116977778440511, 0.030153689607045786, 0.0, 0.024471741852423234, 0.09549150281252627, 0.20610737385376343, 0.3454915028125263, 0.5, 0.6545084971874737, 0.7938926261462365, 0.9045084971874737, 0.9755282581475768, 1.0])
    if length(params)!=8
        error("ERROR in cpgen08n: unsupported number of parameters (must be 8, but ",length(params)," are given)");
    end

    #define decoder weights
    W3 = [0.8015374 -0.5990157 -0.18087512 0.07731245 -0.61930573 0.9293534 -0.19602673 0.13007699; -0.29019997 0.9920272 0.71272105 -0.7493334 2.1605203 1.3799862 -0.6781261 0.3144694; 0.3120384 0.9963041 -0.021765098 -0.30608675 0.22327766 -0.14415248 -0.040516786 0.7197953; 0.48338494 0.114065066 -1.1007736 -1.3659017 0.74279034 -0.8320154 0.41117513 -0.27497062; -0.03915875 -1.6294022 1.2582008 0.69180304 0.57606536 -0.2581074 0.08828592 1.1297721; 0.4388396 -0.9549064 -0.5630029 0.22292994 -0.42703012 -0.20571807 0.13356975 0.8022728; 0.052747495 -0.3368468 0.4924171 -0.069299094 -0.027123263 -0.63098717 0.8597031 -1.2650431; -0.49966058 -0.08987702 -0.24684016 -0.92968196 -0.79355043 -0.5527268 -0.9303465 -0.609024; 0.56279653 -0.47440648 0.32607323 0.36880827 0.31699705 -0.24380437 -0.45596933 -0.35867018; -0.38696793 -0.5112233 -0.88200855 -1.1777732 0.110958196 -0.1616046 0.6527641 0.30737633; -0.4410916 0.19914505 -1.8297246 0.58823967 -0.48786512 -0.080890626 -0.31103334 -0.37956002; 0.51055956 -0.22722387 -0.73076046 0.6376682 0.42858747 0.84136873 1.1915189 -0.16951586];
    b3 = [-0.2156199, -0.49294743, 0.10166748, -0.15939446, 0.35445994, -0.062598385, -0.08431804, 0.50139964, 0.108566605, 0.43357626, 0.46754467, 0.13667274];
    W4 = [0.46298432 0.05487591 0.92299837 0.4495963 -0.004154214 1.8390393 -1.0655475 0.30677167 -2.101477 0.5998454 0.03361331 0.46520144; 1.6305175 0.11217032 0.57364523 0.38772467 0.72628886 -0.56692594 1.1045278 -0.66148275 -2.0011213 -0.817466 0.99658155 -0.8364645; 0.18854636 0.8501607 0.22556071 -1.3277937 -0.4875738 1.3204964 1.445012 -0.5932118 -0.07199024 0.6174115 0.4869855 -0.34616622; -0.3623916 0.4524758 -0.49128267 -0.9988175 -0.33110884 1.0453947 0.4875834 0.24946986 -0.13157608 0.1590036 -0.23223439 0.7577924; 0.0965649 -0.48250452 -1.3817704 0.13922286 0.21137147 -0.51885486 -0.8035133 0.2484077 -0.32427886 -0.31352815 -0.69997966 1.0101993; 1.3571289 -0.43798366 -1.2108278 0.76176196 0.8043557 -1.1001334 -0.56415844 -0.09558436 -0.23842992 -0.44448146 -0.09486207 0.090029284; 1.2223284 0.13774107 -1.2479779 0.5160777 0.49340913 0.3432822 -0.19075334 -0.06873247 0.2462638 -0.56910115 0.12215397 -0.43131852; 0.24817738 0.035547975 -1.4594449 0.18602628 -0.15844753 1.5079646 -0.1716524 -0.32755655 0.6259562 -0.5258306 0.06233403 -0.28183416; -0.5319181 -0.66467136 -0.33306456 0.62747324 0.32585028 0.28222403 -0.021394959 -0.14376217 -0.82706434 -0.7031471 0.19796063 1.0532387; -0.61770123 -0.3153936 0.7783516 -0.03124668 0.63449305 -0.23114257 -0.47163412 2.9473145 0.04140415 -0.5501679 -0.49426335 1.9287983; 0.5467197 0.284589 0.4099826 0.5798313 -0.41161525 -0.88769466 -0.7412688 -0.3737598 0.7581543 -0.2589115 1.4092202 0.015682686; -1.3567795 0.041861802 0.77765805 -0.85450524 -1.0396692 1.148504 -1.1717135 0.8170388 1.9231626 0.15324353 -1.8504131 1.1310815; 0.247526 0.18262105 1.4685688 -0.003622221 0.9597872 -0.30669418 0.039034065 -0.15295282 0.43790862 0.053462002 0.32271734 -0.25999004; 0.5100938 0.04978484 1.0571047 0.6711011 1.6574868 -0.9057639 -0.024199104 -0.31598258 -0.7219456 -0.47063053 0.8988295 -0.52387077; -0.335661 -0.15945472 0.03990918 0.71332747 1.1463879 -0.59195304 -0.879647 -0.08470902 -0.8485535 -0.7574944 0.33112985 -0.14858116; -1.1876434 -0.34438142 -0.7575905 0.24489324 0.022507783 -0.0014716168 -1.6437414 0.095365904 -0.1490304 -0.5308374 -0.58065677 0.16010721; -1.1522465 -0.1436593 -0.38389352 -0.6709555 -0.96265644 0.43724823 -1.2833337 -0.10771872 0.91867363 0.31725866 -0.7919274 -0.14657986; -0.1290102 0.08474075 0.6969439 -1.0732764 -0.9293083 0.1314546 -0.1562085 -0.32019627 0.96871644 0.7339256 -0.19261773 -0.5633551; 1.483871 -0.46445262 1.311482 0.40102404 0.5280848 -1.4111837 0.31014872 0.37155768 -1.4324379 -1.0280437 0.35404935 -0.1889745; 0.48745742 0.04596589 0.9921412 0.38951686 -0.004798859 1.7941666 -1.033959 0.2786738 -2.0161383 0.69613683 0.075280525 0.4447912];
    b4 = [0.036489867, 0.7901641, 0.24798653, -0.2605599, -0.49974605, -0.50623745, -0.31224188, -0.26184863, -0.45990568, -0.8503921, -0.23336938, -0.2969707, -0.4439577, -0.29128766, -0.08473854, 0.09436655, 0.13930129, 0.06559583, 0.24155912, -0.02949089];
    cpi = W4*tanh.(W3*params.+b3).+b4;      #pressure distribution on the default 20 nodes

    #interpolate pressure distribution on the given nodes
    cp = zeros(length(xcp),2);
    cp[:,1] = xcp;
    cp[:,2] = cpi;

    #
    #TODO: cp[:,2] = interp1(...)
    #
    #LEidx = findfirst(abs.(xcp).==minimum(abs.(xcp)));    #index of the leading edge (the point that separates the upper from the lower surface)
    #cpiLEidx = 10;
    #...

    return cp;
end
