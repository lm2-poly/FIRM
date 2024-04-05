%Function calling the python script to generate a mesh of a parallel slits 
% kirigami pattern on a rectangle using gmsh (also called Ribbon kirigami)
%Based on Marzin et al. (2022)'s description
%T. Marzin, K. Le Hay, E. de Langre, and S. Ramananarivo, Flow-induced deformation of kirigami sheets, Physical Review
%Fluids 10.1103/PhysRevFluids.7.023906 (2022)
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

function [coords, connec] = mesh_parallel_slits_kirigami(L, W, Ls, dx, dy, cut_thickness, t, mesh_size, show)
    %Inputs:
    %   L: Length of the plate
    %   W: Width of the plate
    %   Ls: Length of the slits
    %   dx: Width between the slits in the x direction
    %   dy: Width between the slits in the y direction
    %   cut_thickness: Thickness of the slits (should be small)
    %   t: Thickness of the sheet
    %   mesh_size: Size of the elements
    %   show: To display or not
    %
    %Outputs:
    %   coords: Nodes coordinates matrix
    %   connec: Triangulation
    if show
        [flag, ~] = system("python ../../mesh/mesh_parallel_slits_kirigami.py " + L + " " + W + " " + Ls + " " + dx + " " + dy + " " + cut_thickness + " " + t + " " + mesh_size + " True");
    else
        [flag, ~] = system("python ../../mesh/mesh_parallel_slits_kirigami.py " + L + " " + W + " " + Ls + " " + dx + " " + dy + " " + cut_thickness + " " + t + " " + mesh_size + " False");
    end
    
    if flag ~= 0
        error("Mesh was not produced");
    end
    coords = table2array(readtable("nodes.csv"));
    connec = table2array(readtable("connec.csv"));
end

