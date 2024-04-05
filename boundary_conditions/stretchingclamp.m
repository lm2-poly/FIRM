%Applies the boundary conditions required in a stretching test
% Could currently be replaced by bothsideclamp
%Author: Danick Lamoureux
%Creation date: 2023-05-23

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [K, r] = stretchingclamp(K, r, nodes, nonmembernodes)
    %Inputs:
    %   K: Stiffness matrix
    %   r: Load vector
    %   nodes: nodes to clamp
    %   nonmembernodes: Nodes that are not part of an element - clamped
    %
    %Outputs:
    %   K: Modified stiffness matrix
    %   r: Modified load vector
    indices = [];
    for i = 1:length(nodes)
        for k = 1:6
            m = 6*(nodes(i)-1)+k;
            K(m,:) = 0;
            K(m,m) = 1.;
            r(m) = 0.;
        end
    end
    for i = 1:length(nonmembernodes)
        for k = 1:6
            m = 6*(nonmembernodes(i)-1)+k;
            K(m,:) = 0;
            K(m,m) = 1.;
            r(m) = 0.;
        end
    end
end