#======================================================================
    This file contains various utilities to work with airfoils:
    - import/export dat files;
    - naca 4-digits airfoil generator;

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 21/08/2023
======================================================================#


"""
IMPORTAIRFOILFROMFILE reads an airfoil geometry from a dat file
    (myairfoil,myairfoilname) = importairfoilfromfile("myairfoil.dat");
INPUT:
    filein: absolute or relative path to a textual file
OUTPUT:
    airfoil: Nx2 matrix containing the surface points coordinates
    header: first line of the file, usually containing the name of the airfoil
"""
function importairfoilfromfile(filein)
    if !isfile(filein)
        error("ERROR: file '",filein, "' not found");
        return ([],"");
    end
    lines = readlines(filein);
    lines = lines[lines.!=""];
    header = lines[1];
    airfoil = zeros(length(lines)-1,2);
    for i=1:length(lines)-1
        airfoil[i,:] = parse.(Float64, split(lines[i+1]));
    end
    return (airfoil, header);
end


"""
EXPORTAIRFOILDAT saves an airfoil geometry as a DAT file
    exportairfoildat(myairfoil,"my airfoil","myairfoil.dat");
INPUT:
    airfoil: Nx2 matrix containing the surface points coordinates
    header: first line of the file, usually containing the name of the airfoil
    fileout: absolute or relative path of the target file
OUTPUT:
    none
"""
function exportairfoildat(airfoil, header, fileout)
    fileio = open(fileout, "w");
    println(fileio, header);
    for i=1:size(airfoil,1)
        println(fileio, airfoil[i,1]," ",airfoil[i,2]);
    end
    close(fileio);
end


"""
GENERATENACA4AIRFOIL generates a four-digit NACA airfoil
    (myairfoil,myairfoilheader) = generatenaca4airfoil("4412",100);
INPUT:
    digits: String of 4 numbers that defines the airfoil
    N: number of points used to discretize the airfoil surface
OUTPUT:
    airfoil: Nx2 matrix containing the surface points coordinates (TE->upper->LE->lower->TE)
    header: formatted name of the airfoil
"""
function generatenaca4airfoil(digits::String, N)
    #specification available at: https://en.wikipedia.org/wiki/NACA_airfoil
    m = parse(Float64, string(digits[1]))/100;          #maximum camber
    p = parse(Float64, string(digits[2]))/10;           #maximum camber location
    t = parse(Float64, string(digits[3:4]))/100;        #maximum thickness

    #generate a symmetric airfoil
    nupper = floor(Int, N/2);                           #number of nodes on the upper surface
    nlower = N-nupper+1;                                #number of nodes on the lower surface (add one point because the origin will be discarded later)
    #xupper = @. 0.5*(1-cos(pi*(2*(1:nupper)-1)/(2*nupper-1)));     #Gauss-Chebyshev-Lobatto nodes
    #xlower = @. 0.5*(1-cos(pi*(2*(1:nlower)-1)/(2*nlower-1)));
    #x = 0.5.-0.5.*cos.(LinRange(0,2*pi,N));             #cosine-spaced nodes
    #xupper = x[1:nupper];
    #xlower = x[nupper+1:end];
    xupper = 0.5.-0.5.*cos.(LinRange(0,pi,nupper));     
    xlower = 0.5.-0.5.*cos.(LinRange(0,pi,nlower));
    yt(x) = @. 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x^2+0.2843*x^3-0.1036*x^4);  #equation for a symmetrical 4-digits airfoil
    #yupper = yt.(xupper);
    #ylower = yt.(xlower);

    #superpose the camber line
    function yc(x)                                      #equation for the mean camber line
        if p==0
            #symmetric airfoil
            return 0*x;
        end

        if 0<=x<=p
            return @. (m/(p^2))*(2*p*x-x^2);
        else
            return @. (m/((1-p)^2))*(1-2*p+2*p*x-x^2);
        end
    end
    function θ(x)                                       #equation for the camber line slope angle
        if p==0
            #symmetric airfoil
            return 0*x;
        end

        if 0<=x<=p
            return @. atan((2*m/(p^2))*(p-x));
        else
            return @. atan((2*m/((1-p)^2))*(p-x));
        end
    end
    xU = @. xupper-yt.(xupper)*sin(θ(xupper));
    xU[end] = 1.0;
    yU = @. yc(xupper)+yt.(xupper)*cos(θ(xupper));
    yU[end] = 0.0;
    xL = @. xlower+yt.(xlower)*sin(θ(xlower));
    xL[end] = 1.0;
    yL = @. yc(xlower)-yt.(xlower)*cos(θ(xlower));
    yL[end] = 0.0;

    #append points in a counter-clockwise manner
    airfoil = zeros(N,2);
    header = "NACA "*digits*" airfoil";
    airfoil[1:nupper,1] = reverse(xU);
    airfoil[1:nupper,2] = reverse(yU);
    airfoil[nupper+1:end,1] = xL[2:end];
    airfoil[nupper+1:end,2] = yL[2:end];
    return (airfoil, header);
end


"""
MAXIMUMTHICKNESS calculate the maximum thickness of an airfoil
    t = maximumthickness(airfoil)
INPUT:
    airfoil: Nx2 matrix containing the surface points coordinates
OUTPUT:
    t: Float64 value of the maximum thickness evaluated at the airfoil nodes
"""
function maximumthickness(airfoil)
    LEidx = findfirst(abs.(airfoil[:,1]).==minimum(abs.(airfoil[:,1])));    #index of the airfoil leading edge (the point that separates the upper from the lower surface)
    thickness = zeros(length(airfoil[:,2]));         #thickness of the airfoil on every node, initialized with the y-coordinates of the nodes themselves
    for i=1:LEidx
        #upper surface
        for j=LEidx:lastindex(thickness)-1
            if airfoil[j,1]<=airfoil[i,1]<=airfoil[j+1,1]
                #thickness = yupper - interpolated_ylower
                thickness[i] = airfoil[i,2] - (airfoil[j,2] + (airfoil[i,1]-airfoil[j,1])*(airfoil[j+1,2]-airfoil[j,2])/(airfoil[j+1,1]-airfoil[j,1]));
            end
        end
    end
    for i=LEidx+1:lastindex(thickness)
        #lower surface
        for j=1:LEidx
            if airfoil[j+1,1]<=airfoil[i,1]<=airfoil[j,1]
                #thickness = interpolated_yupper - ylower
                thickness[i] = airfoil[i,2] - (airfoil[j+1,2] + (airfoil[i,1]-airfoil[j+1,1])*(airfoil[j,2]-airfoil[j+1,2])/(airfoil[j,1]-airfoil[j+1,1]));
                break;
            end
        end
    end
    return maximum(thickness);
end
