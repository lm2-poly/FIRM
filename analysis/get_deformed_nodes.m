%Deform the nodes
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

function def_nodes = get_deformed_nodes(coords, displacements)
    %Inputs:
    %   coords: coordinates of the nodes initially
    %   displacements: Displacement matrix
    %Outputs:
    %   def_nodes: Deformed coordinates of the nodes
    def_nodes = coords + displacements(:,[1,2,3]);
end