%Unit test to verify the finite element solver and corotational implementation
%Solves the deflection of a one-sided clamped plate under uniform follower load
%This test is based on Caselli and Bisegna (2013)
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

%Comparison models
modelw = "modelw.csv";
modelw = table2array(readtable(modelw));

modelu = "modelu.csv";
modelu = table2array(readtable(modelu));

% Geometry
Lx = 10; %Length of the plate
Ly = 1; % Width of the plate
h = 10E-1; % Thickness

% Material properties
E = 1.2E6; % Young's modulus
nu = 0.; % Poisson's coefficient

% Mesh
%Reuse old meshes
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));

%Remesh your own
[coord, connec] = mesh_rectangle(Lx, Ly, 0.1, false);

%What node is the tip of the plate
for i = 1:length(coord)
    if coord(i,1) == Lx
        tip = i;
    end
end

%Initial guess (is not really important)
x0 = zeros(length(coord)*6, 1);

%Defining the applied loads
I = Ly*h^3/12;
q0 = E*I/Lx^3;
F = linspace(0, 40*q0, 100);

%Data to follow during incrementation
w = [];
u = [];

for m = 1:length(F)
    %Current load applied
    f = F(m);

    %Load (follower)
    f_load = zeros(length(coord), 6);
    for i = 1:length(coord)
        f_load(i, 3) = f;
    end
    
    %No dead load
    d_load = zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solve using Newton
    [x0, f0, r0] = newton(coord, connec, dofTable, x0, @(u) d_load, @(u) f_load,...
        @(K, r) onesideclamp(coord, K, r), E, nu, h, 1E-9, 100, true, true, 0.4, false);

    %FEM solution
    [x, disp_norms] = get_displacements(x0); 
    [f0, ~] = get_displacements(f0); 
    [r0, ~] = get_displacements(r0); 
    def_nodes = get_deformed_nodes(coord, x);
    
    %Saving the data
    savedata("data/solution_"+num2str(m-1)+"_P_"+num2str(f)+".mat", coord, connec, x, [], f0, r0);
    vtkwrite("data/vtk/solution_"+num2str(m-1)+"_P_"+num2str(f)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f0(:,1), f0(:,2), f0(:,3), 'vectors', 'moment', f0(:,4), f0(:,5), f0(:,6),...
        'vectors', 'residualLin', r0(:,1), r0(:,2), r0(:,3), 'vectors', 'residualRot', r0(:,4), r0(:,5), r0(:,6));
    
    %Tip deflections
    w(end+1) = x(tip, 3);
    u(end+1) = -x(tip,1);
    fprintf("Load = %e\t Deflection = %e\n", f, w(end));
end

%Plotting the tip deflections (to compare with the article)
f1 = figure();
hold on;
plot(w/Lx, F*Lx^3/(E*I), 'blue');
scatter(modelw(:,1), modelw(:,2), 'blue');
plot(u/Lx,F*Lx^3/(E*I), 'red');
scatter(modelu(:,1), modelu(:,2), 'red');
xlabel('$\frac{\Delta}{L}$', 'interpreter', 'latex');
ylabel('$\frac{FL^3}{EI}$', 'interpreter', 'latex');
legend('$w_{tip, FEM}$', '$w_{tip, ref}$', '$-u_{tip, FEM}$', '$-u_{tip, ref}$', 'interpreter', 'latex', 'location','best');
xlim([0,1]);
ylim([0,40]);
f1.Units = "centimeters";
f1.Position = [1 1 9, 5];
figureSaver(f1, "followerDistributionCurve.pdf");
