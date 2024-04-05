# FIRM
Flow-Induced Reconfiguration Model based on the Corotational Finite Element Method

## Description
This software is based on the implementation of the corotational framework
for shells from Caselli & Bisegna (2013). Moreover, the finite element assembly
solver is based on an anynonymous contributor's implementation. A relaxed Newton-Raphson
and an arc-length method solver is implemented to allow for the resolution of 
flow-induced reconfiguration problems. Analysis functions are given as starting points
for further analysis. Moreover, the addition of computing loads is presented.
Some meshing files are given and are based off a working installation of GMSH
using python. Multiple verification cases have been tested and compared with
analytical/numerical/experimental results and give good agreement, suggesting
an adequate formulation of a simple reconfiguration solver.

## Prerequisites
You will need a functional MATLAB license and installation. Octave is currently
not supported due to the lack of some required functionalities.

(Optional) In order to mesh your structures in a similar manner as has been
done here, a functional installation of python and GMSH will be required.
However, other techniques can be used to generate the relevant meshes.

(Optional) We implemented functions to generate vtk files that can be visualized
using paraview. While installing paraview is not required, it is recommended for
the visualization.

## Making sure the software is installed correctly
If you have reasons to believe the software is malfunctioning, run the following 
files and judge if the FEM results are close enough to what is expected
(either through plot or console prints).

1. Run verification_cases/1-corotationalImplementationVerification/corotationalImplementationTest.m (Verifies corotational implementation)
2. Run verification_cases/2-clampedBoundaryVerification/clampedBoundaryVerification.m (Verifies finite element assembly)
3. Run verification_cases/3-followerDistributionVerification/followerDistributionVerification.m (Verifies follower loads and nonlinear solver)
4. Run verification_cases/4-loadReadVerification/loadReadVerificationStretch.m and loadReadVerificationBend.m (Verifies displacement control and load reading capabilities)
5. Run verification_cases/5-hingedRoof/hingedRoof.m and readresults.m once the simulation goes over the buckling mechanics (verifies the arc-length method)
6. Run verification_cases/6-plateReconfigurationVerification/plateReconfigurationCurve.m (Verifies reconfiguration loads and analysis, might take a long time)

## How to use

Create a new folder for your files and create a matlab script. In this script,
choose the adequate dimensions for your problem and properties (fluid and solid),
and call the adequate meshing tool that returns the coordinates (nodes) and
triangulation (connectivity table). Using a for loop, solve for each load
incrementally and analyse the results using the analysis functions or your own.
Look at verification cases given here for examples of how to run the given software.

## References:

### Meshing
C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009.

### Arc length method
Kadapa, C. (2021). A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics. Engineering Structures, 234, 111755. https://doi.org/10.1016/j.engstruct.2020.111755

### Pressure field formulation
F. Gosselin, E. de Langre, and B. A. Machado-Almeida, Drag reduction of flexible plates by reconfiguration, Journal of
Fluid Mechanics 10.1017/S0022112009993673 (2010).

### Corotational formulation
Caselli, F. and Bisegna, P. (2013), Polar decomposition based corotational framework for triangular shell elements with distributed loads. Int. J. Numer. Meth. Engng, 95: 499-528. https://doi.org/10.1002/nme.4528

### Exporting to VTK
Yeh, J. (2016), vtkwrite: vtkwrite writes 3D Matlab array into VTK file format. https://www.mathworks.com/matlabcentral/fileexchange/47814-vtkwrite-exports-various-2d-3d-data-to-paraview-in-vtk-file-format

## Contact
Author's email: danick.lamoureux@polymtl.ca

## License
Copyright (C) 2024 Danick Lamoureux

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program. If not, see https://www.gnu.org/licenses/.