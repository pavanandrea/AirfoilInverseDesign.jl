#======================================================================
    AirfoilInverseDesign.jl is a simple Julia package that implements
    the Modified Garabedian-McFadden Method (MGM) for the airfoil
    inverse design.
    By providing a target pressure distribution, this package finds
    the airfoil that generates it.
    Also, by parametrizing the pressure distribution, this package
    generates the optimal airfoil for a given application.
    Advantages of the "flow feature parametrization":
    - Lower dimensionality of the searching space;
    - Robust and smooth geometry variations;
    - Avoids locally over-optimized geometries;
    - Easy performance constraints (e.g. on the operating CL value);
    See the documentation for more details.

    Author: Andrea Pavan
    License: MIT
    Date: 21/08/2023
======================================================================#
module AirfoilInverseDesign
export importairfoilfromfile, exportairfoildat, generatenaca4airfoil, maximumthickness;
export mgm, cpgen10h, cpgen08n;
export panel1;
export neldermead;

include("airfoilutils.jl");
include("mgm.jl");
#include("panel.jl");       #already included in "mgm.jl"

end # module AirfoilInverseDesign
