%Reconfiguration load varying with the displacement (follower load) with
%shading (for circular specimen from Gosselin et al. (2010)) - Based on a
%perimeter-reducing formulation
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

function f_load = shaded_reconfiguration_pressure(coord, connec, u, P, fn)
    %Inputs:
    %   coord: Undeformed nodes coordinates
    %   connec: Triangulation
    %   u: Displacement vector
    %   P: Pressure (nominal)
    %   fn: Load function
    %Outputs:
    %   f_load: Load matrix
    
    f_load = zeros(length(coord), 6);
    
    [x, ~] = get_displacements(u);
    part_x = 0.5*x;
    def_nodes = get_deformed_nodes(coord, x);
    partdef_nodes = get_deformed_nodes(coord, part_x);
    [~, n_partdef] = get_node_angles(partdef_nodes, connec);
    [node_angles, n_def] = get_node_angles(def_nodes, connec);
    [face_angles, ~] = get_face_angles(def_nodes, connec);
    for i = 1:length(coord)
        current_def_R = sqrt(def_nodes(i,1)^2+def_nodes(i,2)^2);
        current_undef_R = sqrt(coord(i,1)^2+coord(i,2)^2);

        shade = current_def_R/current_undef_R;
        
        if dot(n_def(i,:), n_partdef(i,:)) > 0
            f_load(i,3) = -P*shade*fn(node_angles(i));
        end
    end
end