%Unit test to verify the arc length implementation
%Solves the buckling of a hemispherical shell
%This test is based on Gorissen et al. (2020)
%B. Gorissen, D. Melancon, N. Vasios, M. Torbati, and K. Bertoldi, Inflatable soft jumper inspired by shell snapping, Science
%Robotics 10.1126/scirobotics.abb1967 (2020)
%Author: Danick Lamoureux
%Creation date: 2023-05-20

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
R = 100E-3; % Radius of the plate
h = R/8.5; %Thickness of the plate
theta = deg2rad(60); %Aperture of the sphere

Rh = R*sin(theta);
H = R - Rh/tan(theta);

% Material properties
E = 3102.75; % Young's modulus
nu = 0.3; % Poisson's coefficient

% Mesh
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
[coord, connec] = mesh_disk(Rh, h, 10E-3, true); %We form the sphere from a disk

for i = 1:length(coord)
    coord(i,3) = real(sqrt(R^2 - coord(i,1)^2 - coord(i,2)^2)) - R*cos(theta);
end

for i = 1:length(coord)
    if (coord(i,1)^2+coord(i,1)^2) < 0.01*Rh^2
        mid = i;
        break;
    end
end

%Initial guess (is not really important)
x0 = zeros(length(coord)*6, 1);
lambda0 = 0.1;

%Maximum load
p = 10;

f_load = zeros(length(coord), 6);
d_load = zeros(length(coord), 6);
f_load(:,3) = -p;

%Number of degrees of freedom
nn = length(coord);
ndof = 6*nn;

% Indexation table for d.o.f
dofTable = zeros(6,nn);
dofTable(:) = 1:1:ndof;

%Solve using arc length
arclength(coord, connec, dofTable, x0, lambda0, @(u) d_load, ...
                @(u, ite) f_load, @(K, r) roundroller(coord, connec, K, r, Rh),...
                E, nu, h, 1, 1500, 50, 1E-12, 1, true, true);
