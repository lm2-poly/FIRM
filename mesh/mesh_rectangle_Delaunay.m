%Mesh rectangle using triangles from basic matlab functions
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

function [nodes, elements] = mesh_rectangle_Delaunay(Lx, Ly, nx, ny)
    %Inputs:
    %   Lx: Length of the rectangle
    %   Ly: Width of the rectangle
    %   nx: Number of nodes in x
    %   ny: Number of nodes in y
    %
    %Outputs:
    %   nodes: Matrix (number of nodes x 3), each line is a node
    %   elements: Triangulation matrix (number of elements x 3), each line
    %   is an element, each column is one of the vertex ID
    
    %Nodes in x
    x = linspace(0, Lx, nx);
    %Nodes in y
    y = linspace(0, Ly, ny);
    
    %Meshgrid
    [X,Y] = meshgrid(x, y);
    elements = delaunay(X,Y);
    n = numel(X);
    nodes = zeros(n, 3);
    for i = 1:n
        nodes(i,:) = [X(i), Y(i), 0.];
    end
end


