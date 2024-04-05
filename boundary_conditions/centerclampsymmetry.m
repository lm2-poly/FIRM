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

function [K, r] = centerclampsymmetry(coord, connec, K, r, Ri, fraction)
    %Inputs:
    %   coord: Node coordinnates
    %   connec: Triangulation
    %   Ri: Cut-off radius
    %   K: Stiffness matrix
    %   r: Load vector
    %   fraction: Denominator of the disk fraction
    %
    %Outputs:
    %   K: Modified stiffness matrix
    %   r: Modified load vector
    
    %Angle of the symmetry end
    theta = 2*pi/fraction;
    
    %Applying boundary conditions
    clamped = false;
    for i = 1:length(coord)
        %Inner clamp
        if abs(coord(i,1))^2 + abs(coord(i,2))^2 <= Ri^2
            clamped = true;
            for k = 1:6
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end
        
        %Symmetry on the x = 0        
        elseif round(atan2(coord(i,2), coord(i,1)), 4) == 0
            for k = [2, 4, 6]
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end        
        %Symmetry on the other side
        elseif round(atan2(coord(i,2), coord(i,1)), 4) == round(theta, 4)
            %Translation perpendicular to the plane = 0
            %-Tx*sin(theta) + Ty*cos(theta) = 0
            n = 6*(i-1)+1;
            m = 6*(i-1)+2;
            if abs(cos(theta)) > abs(sin(theta))
                K(m,:) = 0;
                K(m,m) = cos(theta);
                K(m,n) = -sin(theta);
                r(m) = 0.;
            else
                K(n,:) = 0;
                K(n,m) = cos(theta);
                K(n,n) = -sin(theta);
                r(n) = 0.;
            end
            %Rotations in the plane = 0
            %Rz = 0
            for k = [6]
                m = 6*(i-1)+k;
                K(m,:) = 0;
                K(m,m) = 1.;
                r(m) = 0.;
            end
            %Rx*cos(theta) + Ry*sin(theta) = 0
            m = 6*(i-1)+4;
            n = 6*(i-1)+5;
            if abs(cos(theta)) > abs(sin(theta))
                K(m,:) = 0;
                K(m,m) = cos(theta);
                K(m,n) = sin(theta);
                r(m) = 0.;
            else
                K(n,:) = 0;
                K(n,m) = cos(theta);
                K(n,n) = sin(theta);
                r(n) = 0.;
            end
        end
            
        
        %Non member nodes
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