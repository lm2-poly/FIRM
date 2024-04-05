%Unit test to verify the  implementation for reconfiguration of a more
%complex structure: a disk
%Solves the deflection of a center clamped disk under reconfiguration
%Compares with Schouveiler's model (Schouveiler and Eloy., 2020)
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

close all;
clear all;
clc;

% Geometry
Re = 70E-3; % Length in direction 1 (before applying symmetry, if any)
Ri = 6E-3; % ... in direction 2 (before applying symmetry, if any)
h = 75E-6; % Thickness


% Material properties
%B = 180E-6;
nu = 0.3; % Poisson's coefficient
E = 3E9; %Young's modulus
B = E*h^3/(12*(1-nu^2));

n = 2;

%Load definition
CD = 1.0;
rho = 1.225;

%Meshes
% Mesh
[coord, connec] = mesh_disk(Re, h, 6E-3, true);
%[coord, connec] = mesh_fraction_disk(Re, h, 4E-3, n, true);
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));

%Initial guess
u = zeros(length(coord)*6, 1);

CYmax = 1;
lambda0 = 1E-3;
U = sqrt(log(Re/Ri)*B*CYmax/(CD*rho*Re^3));
P = 0.5*rho*CD*U^2;

%Load
fn = @(a) cos(a)^2;
f_load = @(u) -reconfiguration_pressure(coord, connec, u, P, fn);
d_load = @(u) zeros(length(coord), 6);

 %Number of degrees of freedom
nn = length(coord);
ndof = 6*nn;

% Indexation table for d.o.f
dofTable = zeros(6,nn);
dofTable(:) = 1:1:ndof;

%Solving using Arc length solver
arclength(coord, connec, dofTable, u, lambda0, @(u) d_load(u), ...
                @(u, ite) f_load(u), @(K, r) centerclamp(coord, connec, K, r, Ri),...
                E, nu, h, 1, 1E9, 100, 1E-12, 1, true, true);