%Clamps the center of a structure in (0,0) with a cut-off radius
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

function [K, r] = centerclamp(coord, connec, K, r, Ri)
    %Inputs:
    %   coord: Node coordinnates
    %   connec: Triangulation
    %   Ri: Cut-off radius
    %   K: Stiffness matrix
    %   r: Load vector
    %
    %Outputs:
    %   K: Modified stiffness matrix
    %   r: Modified load vector
    
    %Applying boundary conditions
    clamped = false;
    for i = 1:length(coord)
        if (abs(coord(i,1))^2 + abs(coord(i,2))^2) <= Ri^2
            clamped = true;
            for k = 1:6
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end
        end
        if ismember(i, connec) == false
            for k = 1:6
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end
        end
    end
    if ~clamped
        fprintf("No fixed nodes");
        exit()
    end
end