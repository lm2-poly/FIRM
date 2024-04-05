%Function to save all the data to a .mat file
%Author: Danick Lamoureux
%Date: 2024-02-28

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

%Data to save:
%   Coord
%   Connec
%   Displacements
%   Angles
%   Forces
%   Residual

function savedata(filename, coord, connec, displacements, angles, forces, residual)
    save(filename, 'coord', 'connec', 'displacements', 'angles', 'forces', 'residual');
end