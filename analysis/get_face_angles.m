%Getting the angles made by the faces
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

function [face_angles, face_areas] = get_face_angles(def_nodes, connec)
    %Inputs:    
    %   def_nodes: Deformed nodes coordinates
    %   connec: Triangulation
    %Outputs:
    %   face_angles: Angles of the triangulation faces
    %   face_areas: Areas of the triangulation faces

    face_angles = zeros(length(connec), 1);
    face_areas = zeros(length(connec), 1);
    for i = 1:length(connec)
        G1 = def_nodes(connec(i,1),:);
        G2 = def_nodes(connec(i,2),:);
        G3 = def_nodes(connec(i,3),:);
        
        v1 = G2 - G1;
        v2 = G3 - G1;
        n = cross(v1, v2);
        face_areas(i) = 0.5*norm(n);
        z = [0, 0, 1];
        
        cosa = dot(n, z)/(norm(n)*norm(z));
        face_angles(i) = acos(cosa);
    end
end