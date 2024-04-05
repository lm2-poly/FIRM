%Unit test to verify the  implementation for reconfiguration and kirigami
%Solves the deflection of kirigami sheet under reconfiguration
%Compares with Marzin's data (Marzin et al., 2022)
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

filename = "MarzinEtAl.csv";
model = table2array(readtable(filename));

% Geometry
cut_thickness = 1E-5;

t = 100E-6;
Ls = 17.2E-3;
dx = 4.4E-3/2;
dy = 3.8E-3;

Nx = 5;
Ny = 31;

L = (Ny+1)*dy;
W = Nx*(Ls+2*dx);
Ssheet = L*W;
S = 0.15*0.15;

% Material properties
E = 4E9; % Young's modulus
nu = 0.3; % Poisson's coefficient

K = 20*E*t^3*dy*Nx/((Ls-2*dx)^3*2*Ny);

%Load definition
D = E*t^3/(12*(1-nu^2));
CD = 2.0;
rho = 997;

%Meshes
%mesh_size = 2E-3;
%coord = table2array(readtable("nodes.csv"));
%connec = table2array(readtable("connec.csv"));
mesh_size = 10E-3;
[coord, connec] = mesh_parallel_slits_kirigami(L, W, Ls, dx, dy, cut_thickness, t, mesh_size, true);

for i = 1:length(coord)
    if coord(i,2) == L
        tip = i;
        break;
    end
end

%Initial guess
u = zeros(length(coord)*6, 1);

%Flow velocity (Cauchy numbers)
CY = logspace(-4, 3, 300);
relax = logspace(log10(0.05), log10(0.1), length(CY));

%Reconfiguration number table
R = zeros(size(CY));
figure(3);
hold on;


scatter(model(:,1), model(:,2));
for k = 1:length(CY)
    Cy = CY(k);
    
    %Actual velocity and pressure
    U = sqrt(K*Cy/(rho*W));
    P = 0.5*rho*CD*U^2;
    %Load
    %f_load = @(u) reconfiguration_pressure(coord, connec, u, P, @(a) cos(a)^2);
    f_load = @(u) reconfiguration_pressure_blockage(coord, connec, u, P, @(a) cos(a)^2, S);
    d_load = @(u) zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solving using Newton
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load(u), @(u) f_load(u), ...
        @(K, r) bothsidesclamp(coord, connec, L, K, r), E, nu, t, 1E-8, 250, true, true, 0.15, false);
    
    %Saving the data
    [x, ~] = get_displacements(u);
    [f, ~] = get_displacements(f);
    [r, ~] = get_displacements(r);
    
    def_nodes = get_deformed_nodes(coord, x);
    
    [nodes_angles, ~] = get_node_angles(def_nodes, connec);
    
    savedata("variableblocked/solution_"+num2str(k)+"_U_"+num2str(U)+".mat", coord, connec, x, nodes_angles, f, r);
    vtkwrite("variableblocked/vtk/solution_"+num2str(k)+"_U_"+num2str(U)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f(:,1), f(:,2), f(:,3), 'vectors', 'moment', f(:,4), f(:,5), f(:,6),...
        'vectors', 'residualLin', r(:,1), r(:,2), r(:,3), 'vectors', 'residualRot', r(:,4), r(:,5), r(:,6),...
        'scalars', 'angles', nodes_angles);

    %FEM solution
    [x, ~] = get_displacements(u);
    y(k) = max(abs(x(:,3)))/L;
    fprintf("Evaluation = %i of %i \t CY = %e cm/s \t Ymax = %e\n", k, length(CY), Cy, y(k));

    figure(3);
    hold on;
    scatter(Cy, y(k), 'k');
    xlabel('$C_Y$');
    set(gca, 'XScale', 'log')
    ylabel('$\mathcal{Y}$');

end
