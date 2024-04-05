%Unit test to verify the  implementation for reconfiguration
%Solves the deflection of disk cut along multiple radii
%Compares with Gosselin's model (Gosselin et al., 2009)
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

filename = "modelDisk.csv";
model = table2array(readtable(filename));

% Geometry
Re = 3.7;
Ri = 0.9;
N = 36;
tcut = 0.001;
h = 1E-2; % Thickness

% Material properties
B = 404E-3;
nu = 0.0;
E = B*12*(1-nu^2)/h^3; % Young's modulus

%Load definition
beta = Re/Ri;
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 1.225*(10^(-3))^3;

%Meshes
% Mesh
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
[coord, connec] = mesh_slitted_disk(Re, Ri, N, tcut, 0.15, true);

%Initial guess
u = zeros(length(coord)*6, 1);

%Flow velocity (Cauchy numbers)
CY = logspace(-1, 3, 100);
CY(end+1) = 1;
CY(end+1) = 10;
CY(end+1) = 100;
CY = unique(CY);

%Reconfiguration number table
R = zeros(size(CY));
figure(3);
hold on;
plot(model(:,1), model(:,2));
for k = 1:length(CY)
    Cy = CY(k);
    
    %Actual velocity and pressure
    U = sqrt(2*D*Cy/(CD*rho*(Re-Ri)^3*beta));
    P = 0.5*rho*CD*U^2;

    %Load
    %f_load = @(u) reconfiguration_pressure(coord, connec, u, P, @(a) cos(a)^2);
    f_load = @(u) shaded_reconfiguration_pressure(coord, connec, u, P, @(a) cos(a)^2);
    d_load = @(u) zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solving using Newton
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load(u), @(u) f_load(u), ...
        @(K, r) centerclamp(coord, connec, K, r, Ri), E, nu, h, 1E-12, 100, false, true, 0.75, false);
    
    %Saving the data
    [x, ~] = get_displacements(u);
    [f, ~] = get_displacements(f);
    [r, ~] = get_displacements(r);
    
    def_nodes = get_deformed_nodes(coord, x);
    
    [nodes_angles, ~] = get_node_angles(def_nodes, connec);
    
    savedata("shaded/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".mat", coord, connec, x, nodes_angles, f, r);
    vtkwrite("shaded/vtk/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f(:,1), f(:,2), f(:,3), 'vectors', 'moment', f(:,4), f(:,5), f(:,6),...
        'vectors', 'residualLin', r(:,1), r(:,2), r(:,3), 'vectors', 'residualRot', r(:,4), r(:,5), r(:,6),...
        'scalars', 'angles', nodes_angles);
    
    %FEM solution
    R(k) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    fprintf("Evaluation = %i of %i \t Cauchy = %e\t Reconfiguration number = %e\n", k, length(CY), Cy, R(k));
    %Plotting the reconfiguration number curve
    figure(3);
    hold on;
    scatter(CY, R, 'k');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('$C_D C_Y$');
    ylabel('$\mathcal{R}$');
    legend("Gosselin's slitted disk model", "Gosselin's plate model", "FEM", 'location', 'best');

end
