%Function calling the python script to generate a mesh of a perforated plate using gmsh
%Based on Jin et al. (2020)'s description
%Y. Jin, J.-T. Kim, S. Cheng, O. Barry, and L. P. Chamorro, On the distinct drag, reconfiguration and wake of perforated
%structures, Journal of Fluid Mechanics 10.1017/jfm.2020.98 (2020)
%Author: Danick Lamoureux
%Creation date: 2023-05-23

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [coords, connec] = mesh_perforated_plate(L, W, a, w, mesh_size, show)
    %Inputs:
    %   L: Length of the plate
    %   W: Width of the plate
    %   a: Side length of the pockets
    %   w: Width between the pockets
    %   mesh_size: Size of the elements
    %   show: To display or not
    %
    %Outputs:
    %   coords: Nodes coordinates matrix
    %   connec: Triangulation
    if show
        [flag, ~] = system("python ../../mesh/mesh_perforated_plate.py " + L + " " + W + " " + a + " " + w + " " +  mesh_size + " True");
    else
        [flag, ~] = system("python ../../mesh/mesh_perforated_plate.py " + L + " " + W + " " + a + " " + w + " " +  mesh_size + " False");
    end
    
    if flag ~= 0
        error("Mesh was not produced");
    end
    coords = table2array(readtable("nodes.csv"));
    connec = table2array(readtable("connec.csv"));
end

