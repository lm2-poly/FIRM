%Unit test to verify the  implementation for reconfiguration
%Solves the deflection of symmetrical multistable beam
%Author: Danick Lamoureux
%Creation date: 2024-03-15

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

close all;
clear all;
clc;

% Geometry
Lx = 50/2; % Length in direction 1 (before applying symmetry, if any)
Ly = 10; % ... in direction 2 (before applying symmetry, if any)
h = 0.1; % Thickness

dh = 20;

% Material properties
E = 3E3; % Young's modulus
nu = 0.; % Poisson's coefficient

%Load definition
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 1.225*(10^(-3))^3;

%Meshes
%Reuse meshes
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
% Mesh
mesh_size = 1;
[coord, connec] = mesh_rectangle(Lx, Ly, mesh_size, true);
L = Lx;
for i = 1:length(coord)
    s = coord(i,1)/L;
    coord(i,3) = dh*(3*s-s^3)/2;
end

%Initial guess
u = zeros(length(coord)*6, 1);

CYmax = 1.;
lambda0 = 1E-2;

%Actual velocity and pressure
U = sqrt(2*D*CYmax/(CD*rho*(2*Lx)^3));
P = 0.5*rho*CD*U^2;

%Load
f_load = @(u) reconfiguration_pressure(coord, connec, u, P, @(a) cos(a)^2);
d_load = zeros(length(coord), 6);
d_load = @(u) d_load;

%Number of degrees of freedom
nn = length(coord);
ndof = 6*nn;

% Indexation table for d.o.f
dofTable = zeros(6,nn);
dofTable(:) = 1:1:ndof;

%Solving using Newton
arclength(coord, connec, dofTable, u, lambda0, @(u) d_load(u), ...
                @(u, ite) f_load(u), @(K, r) pivotsymmetry(coord, connec, K, r, Lx),...
                E, nu, h, 1, 1E4, 100, 1E-9, 1, true, true);