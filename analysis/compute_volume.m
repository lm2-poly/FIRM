%Compute the volume a surface makes in regards to a reference point (used for the buckling shell example)
%Author: Danick Lamoureux
%Creation date: 2024-02-08

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function V = compute_volume(connec, def_nodes, ref)
    %Inputs:
    %   connec: Triangulation
    %   def_nodes: Deformed nodes matrix
    %   ref: Point of reference
    %Output:
    %   V: Volume
    V = 0;
    for i = 1:length(connec(:,1))
        G1 = ref;
        G2 = def_nodes(connec(i,1),:);
        G3 = def_nodes(connec(i,2),:);
        G4 = def_nodes(connec(i,3),:);
        
        v1 = G2 - G1;
        v2 = G3 - G1;
        v3 = G4 - G1;
        
        axb = cross(v2, v3);
        c = v1;
        
        dV = abs(dot(c, axb))/6;
        V = V + dV;
    end
end