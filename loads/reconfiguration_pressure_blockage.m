%Reconfiguration load varying with the displacement (follower load)
%This formulation is based on a momentum conservation equation developed by
%Gosselin et al. (2010) - The load is aligned with the vertical axis z
%Here, we adapt the load to a variable blockage formulation through mass conservation
%based on T. Marzin, K. Le Hay, E. de Langre, and S. Ramananarivo, Flow-induced deformation of kirigami sheets, Physical Review
%Fluids 10.1103/PhysRevFluids.7.023906 (2022)
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

function f_load = reconfiguration_pressure_blockage(coord, connec, u, P, fn, S)
    %Inputs:
    %   coord: Undeformed nodes coordinates
    %   connec: Triangulation
    %   u: Displacement vector
    %   P: Pressure (nominal)
    %   fn: Load function
	%	S: Canal area
    %Outputs:
    %   f_load: Load matrix
    
    
    f_load = zeros(length(coord), 6);
    z = [0,0,1];
    [x, ~] = get_displacements(u);
    part_x = 0.5*x;
    def_nodes = get_deformed_nodes(coord, x);
    partdef_nodes = get_deformed_nodes(coord, part_x);
    [~, n_partdef] = get_node_angles(partdef_nodes, connec);
    [node_angles, n_def] = get_node_angles(def_nodes, connec);
    [face_angles, face_areas] = get_face_angles(def_nodes, connec);
    
    frontal_area = sum(face_areas.*cos(face_angles));
    factor = S/(S-frontal_area);
    P = P*factor^2;
    
    for i = 1:length(coord)
        %Making sure the node normal still faces the same direction as
        %before, otherwise the load is neglected - aims at removing
        %reversal deformations
        if dot(n_def(i,:), n_partdef(i,:)) > 0
            f_load(i,3) = -P*fn(node_angles(i));
        end
    end
end