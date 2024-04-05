%Function calling the python script to generate a mesh of a fraction of a disk using
%gmsh
%Author: Danick Lamoureux
%Creation date: 2023-05-21

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [coords, connec] = mesh_fraction_disk(R, t, mesh_size, fraction, show)
    %Inputs:
    %   R: Radius of the circle
    %   t: Thickness of the disk
    %   mesh_size: Size of the elements
    %   fraction: Denominator of the fraction of the disk (2 means half a
    %   disk, 3 a third of the disk, etc).
    %   show: To display or not
    %
    %Outputs:
    %   coords: Nodes coordinates matrix
    %   connec: Triangulation
    if show
        [flag, ~] = system("python ../../mesh/mesh_fraction_disk.py " + R + " " + t + " " + mesh_size + " " + fraction + " True");
    else
        [flag, ~] = system("python ../../mesh/mesh_fraction_disk.py " + R + " " + t + " " + mesh_size + " " + fraction + " False");
    end
    
    if flag ~= 0
        error("Mesh was not produced");
    end
    coords = table2array(readtable("nodes.csv"));
    connec = table2array(readtable("connec.csv"));
end

