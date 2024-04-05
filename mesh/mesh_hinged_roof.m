%Function calling the python script to generate a mesh of a hinged roof using
%gmsh
%Author: Danick Lamoureux
%Creation date: 2023-12-07

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [coords, connec] = mesh_hinged_roof(L, theta, R, mesh_size, show)
    %Inputs:
    %   L: Length of the plate
    %   theta: Angle of the roof
    %   R: Radius of the roof
    %   mesh_size: Size of the elements
    %   show: To display or not
    %
    %Outputs:
    %   coords: Nodes coordinates matrix
    %   connec: Triangulation
    if show
        [flag, ~] = system("python ../../mesh/mesh_hinged_roof.py " + L + " " + theta + " " + R + " " + mesh_size + " True");
    else
        [flag, ~] = system("python ../../mesh/mesh_hinged_roof.py " + L + " " + theta + " " + R + " " + mesh_size + " False");
    end
    
    if flag ~= 0
        error("Mesh was not produced");
    end
    coords = table2array(readtable("nodes.csv"));
    connec = table2array(readtable("connec.csv"));
end

