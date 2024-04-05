%Getting the best fit plane equations through the nodes
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

function [node_planes] = get_node_planes(def_nodes, connec)
    %Inputs:
    %   def_nodes: Deformed nodes coordinates
    %   connec: Triangulation
    %Outputs:
    %   node_planes: Nodal planes coefficients

    node_planes = zeros(length(def_nodes),3);
    for i = 1:length(def_nodes)
        connected_nodes = [];
        for j = 1:length(connec)
            if any(connec(j,:) == i)
                for k = 1:3
                    connected_nodes(end+1) = connec(j,k);
                end
            end
        end
        connected_nodes = unique(connected_nodes);
        %This method is taken from C. Audet and W. Hare, Derivative-Free and Blackbox Optimization (Springer, 2017).
        if length(connected_nodes) ~= 0
            Y = [];
            f = [];
            for j = 1:length(connected_nodes)
                Y(:,end+1) = def_nodes(connected_nodes(j),[1,2])';
                f(end+1,1) = def_nodes(connected_nodes(j), 3);
            end
            %Interpolation
            M = [ones(length(connected_nodes),1), Y'];
            a = (M'*M)\(M'*f);
            node_planes(i,:) = [-a(2), -a(3), 1, -a(1)];
        else
            node_planes(i,:) = [0,0,1,0];
    end
end