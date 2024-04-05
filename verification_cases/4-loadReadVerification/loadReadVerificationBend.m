%Unit test to verify the finite element solver and corotational
%implementation and make sure we read adequate loads under bending
%Solves the deflection of a one-sided clamped plate under tip load
%This test is initially based on Caselli and Bisegna (2013), but had some
%modifications to allow for displacement control
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
Lx = 1; %Length of the plate
Ly = 10; % Width of the plate
h = 10E-1; % Thickness

% Material properties
E = 1.2E3; % Young's modulus
nu = 0.; % Poisson's coefficient

% Mesh
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
[coord, connec] = mesh_rectangle(Lx, Ly, 0.5, false);

%What node is the tip of the plate
for i = 1:length(coord)
    if coord(i,2) == Ly
        tip = i;
    end
end

%Initial guess
u = zeros(length(coord)*6, 1);

%Calculate the required displacements
P = linspace(0, 1, 100);
I = Lx*h^3/3;

%Displacements
du = linspace(0, 0.7, 1000)*Ly;

%Data to follow during incrementation
F = zeros(length(du),1);

figure(3);
hold on;

for k = 1:length(du)
    for j = 1:length(coord)
        if coord(j,2) == Ly
            m = 6*(j-1)+3;
            u(m) = du(k);
        end
    end

    %Load (dead)
    d_load = zeros(length(coord), 6);
    
    %No follower load
    f_load = zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solve using Newton under displacement control
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load, @(u) f_load,...
        @(K, r) clamp_dispcontrol(coord, connec, Ly, K, r), E, nu, h, 1E-12, 100, true, true, 1, true);
    
    %Saving the data
    [x, ~] = get_displacements(u);
    [f, ~] = get_displacements(f);
    [r, ~] = get_displacements(r);
    
    savedata("data/solution_"+num2str(k)+".mat", coord, connec, x, [], f, r);
    vtkwrite("data/vtk/solution_"+num2str(k)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f(:,1), f(:,2), f(:,3), 'vectors', 'moment', f(:,4), f(:,5), f(:,6),...
        'vectors', 'residualLin', r(:,1), r(:,2), r(:,3), 'vectors', 'residualRot', r(:,4), r(:,5), r(:,6));
    
    %FEM solution    
    F(k) = 0;
    for i = 1:length(coord)
        if coord(i,2) == Ly
            F(k) = F(k) + f(i,3);
        end
    end
    fprintf("Evaluation = %i of %i \t dL = %e mm\t F = %e N\n", k, length(du), du(k), F(k));

    
end
f1 = figure(3);
hold on;
plot(du, F);
legend("FEM", "location", "best");
xlabel('$w$ [mm]');
ylabel('$F$ [N]');
f1.Units = "centimeters";
f1.Position = [1 1 9, 5];
figureSaver(f1, "bendRead.pdf");
