%Unit test to verify the  implementation for reconfiguration
%Solves the deflection of a disk cut along one radius
%Compares with Schouveiler's model and data (Schouveiler & Boudaoud, 2006)
%L. Schouveiler and A. Boudaoud, The rolling up of sheets in a steady flow, J. Fluid Mech. 10.1017/S0022112006000851
%(2006)
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

filename = "SchouveilerModel.csv";
model = table2array(readtable(filename));

% Geometry
Re = 100;
Ri = 1;
tcut = 0.1;
h = 0.1; % Thickness

% Material properties
B = 1.92;
nu = 0.3;
E = B*12*(1-nu^2)/h^3; % Young's modulus

%Load definition
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 997*(10^(-3))^3;

%Meshes
% Mesh
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
[coord, connec] = mesh_singleradius_disk(Re, tcut, 8, true);

bestangle = pi;
for i = 1:length(coord)
    if coord(i,1)^2 + coord(i,2)^2 == Re^2
        if abs(atan2(coord(i,2), coord(i,1)) - pi) < bestangle
            masternode = i;
            bestangle = abs(atan2(coord(i,2), coord(i,1)) - pi);
        end
    end
end

%Initial guess
u = zeros(length(coord)*6, 1);

%Flow velocity (Cauchy numbers)
CY = logspace(-1, 2, 100);
CY(end+1) = 1;
CY(end+1) = 10;
CY(end+1) = 100;
CY = unique(CY);

%Reconfiguration number table
R = zeros(size(CY));
figure(3);
hold on;
plot(model(:,1), model(:,3));
for k = 1:length(CY)
    Cy = CY(k);
    
    %Actual velocity and pressure
    U = sqrt(log(Re/Ri)*B*Cy/(CD*rho*Re^3));
    P = 0.5*rho*CD*U^2;

    %Load
    f_load = @(u) -reconfiguration_pressure(coord, connec, u, P, @(a) cos(a)^2);
    d_load = @(u) zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solving using Newton
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load(u), @(u) f_load(u), ...
        @(K, r, u) centerclamp(coord, connec, K, r, Ri), ...
        E, nu, h, 10^(-10), 1000, true, true, 0.75, false);
    %Saving the data
    [x, ~] = get_displacements(u);
    [f, ~] = get_displacements(f);
    [r, ~] = get_displacements(r);
    
    def_nodes = get_deformed_nodes(coord, x);
    
    [nodes_angles, ~] = get_node_angles(def_nodes, connec);
    
    savedata("data/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".mat", coord, connec, x, nodes_angles, f, r);
    vtkwrite("data/vtk/solution_"+num2str(k)+"_CY_"+num2str(Cy)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
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
    legend("Schouveiler's momentum conservation model", "FEM", 'location', 'best');

end
