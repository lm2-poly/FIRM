%This function returns the block diagonal tensor collecting 6 copies of R
%This code is taken as is from an anonymous contributor's code with no modification. All
%credit regarding this function goes to this anonymous collaborator.

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [ R_diese ] = diese( R )
% Returns the block diagonal tensor collecting 6 copies of R

% INPUT:
%  R: rotation tensor (3,3)
%  nn: number of nodes (1,1)

% OUTPUT:
%  R_diese: xxx (18,18)

R_diese = zeros(18, 18);

v = 1:3;
for ii = 1:6
    R_diese(v,v) = R;
    v = v + 3;
end

end