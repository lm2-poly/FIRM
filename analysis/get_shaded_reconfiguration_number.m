%Computes the reconfiguration number with the shade from the circular
%specimen from Gosselin et al. (2010)
%The shading is based on a perimeter-reducing approach
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

function R = get_shaded_reconfiguration_number(coord, connec, u, fn, Ri)
    %Inputs:
    %   coord: Nodes' initial coordinates
    %   connec: Triangulation
    %   u: Displacement vector
    %   fn: Load function
    %   Ri: Inner radius
    %Outputs:
    %   R: Reconfiguration number
    
    totalarea = -pi*Ri^2;
    cosarea = -pi*Ri^2;
    
    [x, ~] = get_displacements(u);
    def_nodes = get_deformed_nodes(coord, x);
    [angles, areas] = get_face_angles(def_nodes, connec);
    
    for i = 1:length(connec)
        xi_undef = (1/3)*(coord(connec(i,1),:) + coord(connec(i,2),:) + coord(connec(i,3),:));
        xi_def = (1/3)*(def_nodes(connec(i,1),:) + def_nodes(connec(i,2),:) + def_nodes(connec(i,2),:));

        current_def_R = sqrt(xi_def(1)^2+xi_def(2)^2);
        current_undef_R = sqrt(xi_undef(1)^2+xi_undef(2)^2);

        shade = current_def_R/current_undef_R;
        
        cosarea = cosarea + areas(i)*fn(angles(i))*abs(cos(angles(i)))*shade;
        totalarea = totalarea + areas(i);
    end
    
    R = cosarea/totalarea;
end