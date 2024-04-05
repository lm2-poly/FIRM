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


%Load definition
CD = 1.0;
rho = 1.225;

%Meshes
% Mesh
%[coord, connec] = mesh_fraction_disk(Re, h, 8E-3, 2, true);
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

%Initial guess
u = zeros(length(coord)*6, 1);

CY = [1E-2:1E-2:1E-1, 1E-1:5E-2:1E0, 1E0:5E-1:1E1, 1E1:2E0:1E2, 1E2:1E1:1E3, 1E3:1E2:1E4];
CY(end+1) = 0.1; CY(end+1) = 1; CY(end+1) = 10; CY(end+1) = 100; CY(end+1) = 1000;
CY = unique(CY);
R = zeros(size(CY));
for k = 1:length(CY)
    Cy = CY(k);
    U = sqrt(log(Re/Ri)*B*Cy/(CD*rho*Re^3));
    P = 0.5*rho*CD*U^2;

    %Load
    fn = @(a) cos(a)^2;
    f_load = @(u) reconfiguration_pressure(coord, connec, u, P, fn);
    d_load = @(u) zeros(length(coord), 6);

     %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solving using Newton
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load(u), @(u) f_load(u), ...
        @(K, r) centerclamp(coord, connec, K, r, Ri), E, nu, h, 10^(-8), 100, true, true, 0.5, false);
    
    %Saving the data
    [x, ~] = get_displacements(u);
    [f, ~] = get_displacements(f);
    [r, ~] = get_displacements(r);
    
    def_nodes = get_deformed_nodes(coord, x);
    
    [nodes_angles, ~] = get_node_angles(def_nodes, connec);
    
    savedata("Newton/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".mat", coord, connec, x, nodes_angles, f, r);
    vtkwrite("Newton/vtk/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f(:,1), f(:,2), f(:,3), 'vectors', 'moment', f(:,4), f(:,5), f(:,6),...
        'vectors', 'residualLin', r(:,1), r(:,2), r(:,3), 'vectors', 'residualRot', r(:,4), r(:,5), r(:,6),...
        'scalars', 'angles', nodes_angles);
    
    %FEM solution
    R(k) = get_reconfiguration_number(coord, connec, u, fn, Ri);
    fprintf("Evaluation = %i of %i \t CY = %e m/s \t Reconfiguration number = %e\n", k, length(CY), CY(k), R(k));
    
    figure(3);
    hold on;
    scatter(Cy, R(k), 'k');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('$C_D C_Y$');
    ylabel('$\mathcal{R}$');
end