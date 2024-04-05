%Unit test to verify the current implementation can perform kirigami
%traction tests
%T. Marzin, Flow-gamis : Interaction of folds and cuts with a flow, Ph.D. thesis (2023), thse de doctorat dirige par Langre,
%Emmanuel de Ingnierie, mcanique et nergtique Institut polytechnique de Paris 2023.
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
cut_thickness = 1E-4;

t = 100E-6;
Ls = 17.2E-3;
dx = 4.4E-3/2;
dy = 3.8E-3;

Nx = 5;
Ny = 31;

L = (Ny+1)*dy;
W = Nx*(Ls+2*dx);

% Material properties
E = 4E9; % Young's modulus
nu = 0.3; % Poisson's coefficient

%Meshes
%coord = table2array(readtable("MarzinEtAl/nodes.csv"));
%connec = table2array(readtable("MarzinEtAl/connec.csv"));
mesh_size = 10E-3;
[coord, connec] = mesh_parallel_slits_kirigami(L, W, Ls, dx, dy, cut_thickness, t, mesh_size, true);

bc_nodes = [];
for i = 1:length(coord)
    if coord(i,2) == 0 || coord(i,2) == L
        bc_nodes(end+1) = i;
    end
end

nonmembernodes = [];
for i = 1:length(coord)
    if ismember(i, connec) == false
        nonmembernodes(end+1) = i;
    end
end

for i = 1:length(coord)
    if coord(i,2) == L
        tip = i;
        break;
    end
end

%Initial guess
u = zeros(length(coord)*6, 1);

%Loads and displacements
du = logspace(-4, log10(L), 300);%(0, L, 300);
relax = logspace(0, -3, length(du));
F = zeros(size(du));
f = zeros(length(coord), 6);
figure(3);
hold on;

for k = 1:length(du)
    for j = 1:length(coord)
        if coord(j,2) == max(coord(:,2))
            m = 6*(j-1)+2;
            u(m) = du(k);
        end
    end

    %Load
    f_load = @(u) zeros(length(coord), 6);
    d_load = @(u) zeros(length(coord), 6);

    %Number of degrees of freedom
    nn = length(coord);
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Solving using Newton
    [u, f, r] = newton(coord, connec, dofTable, u, @(u) d_load(u), @(u) f_load(u), ...
        @(K, r) stretchingclamp(K, r, bc_nodes, nonmembernodes), E, nu, t, 1E-8, 250, ...
        true, true, 0.005, false);
    
    %FEM solution
    [x, disp_norms] = get_displacements(u); 
    [f, ~] = get_displacements(f); 
    [r, ~] = get_displacements(r); 
    def_nodes = get_deformed_nodes(coord, x);
    
    %Saving the data
    savedata("MarzinEtAl/solution_"+num2str(k)+".mat", coord, connec, x, [], f, r);
    vtkwrite("MarzinEtAl/vtk/solution_"+num2str(k)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
        'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
        'vectors', 'force', f(:,1), f(:,2), f(:,3), 'vectors', 'moment', f(:,4), f(:,5), f(:,6),...
        'vectors', 'residualLin', r(:,1), r(:,2), r(:,3), 'vectors', 'residualRot', r(:,4), r(:,5), r(:,6));
    
    dL(k) = max(abs(x(:,2)));
    
    F(k) = 0;
    for i = 1:length(coord)
        if coord(i,2) == L
            F(k) = F(k) + f(i,2);
        end
    end
    fprintf("Evaluation = %i of %i \t dL = %e mm\t F = %e N\n", k, length(du), du(k), F(k));
    figure(3);
    hold on;
    scatter(dL(k), F(k), 'k');
    xlabel('$\Delta L$ [mm]');
    ylabel('$F$ [N]');

end
