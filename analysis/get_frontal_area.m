%Get the area of the sheet perpendicular to the flow (used for tunnel
%testing when blockage is considered)
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

function A = get_frontal_area(coord, connec)
    %Inputs:    
    %   coord: Coordinates of the nodes (deformed or not)
    %   connec: Triangulation
    %Outputs:
    %   A: Total area
    
    A = 0;
    for i = 1:length(connec)
        G1 = coord(connec(i,1),:);
        G2 = coord(connec(i,2),:);
        G3 = coord(connec(i,3),:);
        
        v1 = G2 - G1;
        v2 = G3 - G1;
        n = cross(v1, v2);
        A = A+ 0.5*norm(n);
    end
end