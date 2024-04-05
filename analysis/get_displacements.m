%Get displacement matrix from the displacement vector (can be used for
%forces and residual too as it is analoguous) - does not return rotation
%values
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

function [x, disp_norms] = get_displacements(u)
    %Inputs:
    %   u: Displacement vector (or forces vector)
    %Outputs:
    %   x: Displacement matrix (or force matrix)
    %   disp_norms: Norm of the displacement of each node (or norm of the
    %   forces); does not consider angular displacements (or forces)

    %Displacements
    x = zeros(length(u)/6, 6);
    %Magnitude of displacements
    disp_norms = zeros(length(x), 1);
    for i = 1:length(x)
        for k = 1:6
            m = 6*(i-1) + k;
            x(i,k) = u(m);
        end
        disp_norms(i) = norm(x(i,[1,2,3]));
    end
end