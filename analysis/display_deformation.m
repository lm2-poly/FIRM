%Display the deformed nodes in a surface
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

function display_deformation(connec, def_nodes, color, equalaxis)
    %Inputs:
    %   connec: Triangulation
    %   def_nodes: Deformed nodes
    %   color: The array to color the plot
    %   equalaxis: To put equal axes or not (and other axes options)
    
    %Plotting the final deformation
    trisurf(connec, def_nodes(:, 1), def_nodes(:, 2), def_nodes(:, 3), color, "EdgeAlpha", 0, "EdgeColor", [220,220,220]/256);% 'EdgeColor', 'none');
    shading interp;
    xlabel('x [mm]');
    ylabel('y [mm]');
    zlabel('z [mm]');
    if equalaxis
        axis('equal');
        axis off
        set(gcf,'color','white')
    end
end