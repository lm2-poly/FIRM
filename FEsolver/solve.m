%Solve one iteration of the finite element problem
%Author: Danick Lamoureux
%Creation date: 2023-06-03

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [du, r, f] = solve(coord, connec, dofTable, u, d_load, f_load, bc, E, nu, h)
    % INPUT:
    %  coord: nodes of the current element
    %  connec: Triangulation
    %  dofTable: To convert from vector to matrices
    %  u: Current value of dofs in global frame
    %  d_load: Dead load matrix
    %  f_load: Follower load matrix
    %  bc: Boundary condition function
    %  E: Young's modulus
    %  nu: Poisson's ratio
    %  h: Plate's thickness

    % OUTPUT:
    %  du: nodal value variation (in global frame)
    %  r: Residual vector
    %  f: Load vector

    %Assembling matrix and vector
    [K, r, f] = assemble(coord, connec, dofTable, u, d_load(u), f_load(u), E, nu, h);
    [K, r] = bc(K, r);

    %Solving
    du = -K\r;
end
