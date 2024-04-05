%Get displacement vector from the displacement matrix (can be used for
%forces and residual too as it is analoguous)
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

function [u] = matrix2vec(x)
    %Inputs:
    %   x: Displacement matrix (or forces matrix)
    %Outputs:
    %   u: Displacement vector (or force vector)

    %Displacements
    u = zeros(length(x)*6, 1);
    for i = 1:length(x)
        for k = 1:6
            m = 6*(i-1) + k;
            u(m) = x(i,k);
        end
    end
end