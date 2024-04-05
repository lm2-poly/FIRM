%Applies particular boundary conditions for symmetrical multistable beams
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

function [K, r] = pivotsymmetry(coord, connec, K, r, L)
        %Inputs:
    %   coord: Node coordinnates
    %   K: Stiffness matrix
    %   r: Load vector
    %   L: Length at which to clamp
    %
    %Outputs:
    %   K: Modified stiffness matrix
    %   r: Modified load vector
    
    %Applying boundary conditions
    for i = 1:length(coord)
        used = false;
        for j = 1:length(connec)
            if i == connec(j,1) || i == connec(j,2) || i == connec(j,3)
                used = true;
                break;
            end
        end
        if used
            if coord(i,1) == 0
                for k = 1:3
                    m = 6*(i-1)+k;
                    K(m,:) = 0;
                    K(m,m) = 1.;
                    r(m) = 0.;
                end
            end
            if coord(i,1) == L
                for k = [1,5,6]
                    m = 6*(i-1)+k;
                    K(m,:) = 0;
                    K(m,m) = 1.;
                    r(m) = 0.;
                end
            end
        else
            for k = 1:6
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end
        end
    end
end