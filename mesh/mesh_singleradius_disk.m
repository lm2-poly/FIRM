%Function calling the python script to generate a mesh of a disk with a single cut using
%gmsh as defined by Schouveiler & Boudaoud (2006)
%Based on Schouveiler & Boudaoud (2006)'s description
%L. Schouveiler and A. Boudaoud, The rolling up of sheets in a steady flow, J. Fluid Mech. 10.1017/S0022112006000851
%(2006).
%Author: Danick Lamoureux
%Creation date: 2023-05-20

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [coords, connec] = mesh_singleradius_disk(R, tcut, mesh_size, show)
    %Inputs:
    %   R: External radius of the disk
    %   tcut: Width of the cut
    %   mesh_size: Size of the elements
    %   show: To display or not
    %
    %Outputs:
    %   coords: Nodes coordinates matrix
    %   connec: Triangulation
    if show
        [flag, ~] = system("python ../../mesh/mesh_singleradius_disk.py " + R + " " + tcut + " " +mesh_size + " True");
    else
        [flag, ~] = system("python ../../mesh/mesh_singleradius_disk.py " + R + " " + tcut + " " +mesh_size + " False");
    end
    
    if flag ~= 0
        error("Mesh was not produced");
    end
    coords = table2array(readtable("nodes.csv"));
    connec = table2array(readtable("connec.csv"));
end

