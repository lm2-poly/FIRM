%Unit test to verify the arc-length method
%Solves the deflection of a hinged roof under point load
%This test is based on Caselli and Bisegna (2013) and Kadapa (2021)
%F. Caselli and P. Bisegna, Polar decomposition based corotational framework for triangular shell elements with distributed
%loads: COROTATIONAL FRAMEWORK WITH DISTRIBUTED LOADS, International Journal for Numerical Methods
%in Engineering 10.1002/nme.4528 (2013)
%C. Kadapa, A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for
%nonlinear structural mechanics, Engineering Structures 10.1016/j.engstruct.2020.111755 (2021)
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
L = 254;
R = 2540; % Radius of the plate
theta = 0.1;
h = 6.35; % Thickness

% Material properties
E = 3102.75; % Young's modulus
nu = 0.3; % Poisson's coefficient

% Mesh
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
mesh_size = 8;
[coord, connec] = mesh_hinged_roof(L, theta, R, mesh_size, true);

%Initial guess (is not really important)
x0 = zeros(length(coord)*6, 1);
lambda0 = 0.01;

%Maximum load
F = -1000;

d_load = zeros(length(coord), 6);
realA = 0;
for i = 1:length(coord)
    if coord(i,1) == 0 && coord(i,2) == L
        for j = 1:length(connec)
            if ismember(i, connec(j,:))
                inode = i;
                x1 = coord(i,1); y1 = coord(i,2);
                is = connec(j,:);
                is = is(is~=i);
                x2 = coord(is(1), 1); y2 = coord(is(1),2);
                x3 = coord(is(2), 1); y3 = coord(is(2),2);
                J = [x2 - x1, y2 - y1;
                    x3 - x1, y3 - y1];
                d_load(i,3) = 2*F/abs(det(J));
            end
        end
    end
end

%No follower load
f_load = zeros(length(coord), 6);

%Number of degrees of freedom
nn = length(coord);
ndof = 6*nn;

% Indexation table for d.o.f
dofTable = zeros(6,nn);
dofTable(:) = 1:1:ndof;

%Solve using arc length
arclength(coord, connec, dofTable, x0, lambda0, @(u) d_load, ...
                @(u, ite) f_load, @(K, r) symmetrichinge(coord, connec, L, R, theta, K, r),...
                E, nu, h, 1, 1500, 50, 1E-11, 1, true, true);