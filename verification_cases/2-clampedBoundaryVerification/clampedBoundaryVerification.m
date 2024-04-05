%Unit test to verify the finite element assembly and solver implementation
%without considering large displacements with the corotational formulation
%Solves the deflection of a clamped plate under uniform load in comparison
%with the Kirchoff plate solution.
%This test is based on the solution presented in J. Bleyer, Numerical Tours of Computational Mechanics with FEniCS (2018).
%https://comet-fenics.readthedocs.io/en/latest/demo/reissner_mindlin_quads_reduced_integration/reissner_mindlin_quads.py.html
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
Lx = 1; % Length in direction x
Ly = 1; % Width of the plate in y
h = 1E-3; % Thickness

% Material properties
E = 1000; % Young's modulus
nu = 0.3; % Poisson's coefficient

%Meshes sizes
N = [5,10,20,40,80];

%Data to follow during mesh refinements
diff = [];
nnodes = [];

%Refinement analysis
for m = 1:length(N)
    %Number of nodes for this refinement
    n = N(m);
    
    % Mesh
    nx = n; % Number of nodes in x
    ny = n; % ... y

    %Meshing using matlab's Delaunay
    [coord, connec] = mesh_rectangle_Delaunay(Lx, Ly, nx, ny);
    
    %Number of nodes
    nn = size(coord, 1);
    nnodes(end+1) = nn;

    %Load defined by a function of the thickness
    f = h^3;
    
    %The load is a dead load
    d_load = zeros(length(coord), 6);
    for i = 1:length(coord)
        d_load(i, 3) = f;
    end
    
    %Follower loads are zero
    f_load = zeros(length(coord), 6);

    %Initial guess of the geometry (does not matter for linear analysis)
    x0 = zeros(length(coord)*6, 1);
    
    %Number of degrees of freedom
    ndof = 6*nn;

    % Indexation table for d.o.f
    dofTable = zeros(6,nn);
    dofTable(:) = 1:1:ndof;
    
    %Assemble the stiffness matrix and load vector
    [K, r] = assemble(coord, connec, dofTable, x0, d_load, f_load, E, nu, h);
    
    %Applying the clamp boundary conditions
    [K, r] = allclamp(coord, Lx, Ly, K, r);
    
    %Solving for the displacement field
    u = -K\r;

    %Reference solution
    %Bending stiffness
    D = E*h^3/(12*(1-nu^2));
    %Deflection
    wk = 1.265319087e-3*f/D;
    fprintf("Kirchhoff deflection: %e\n", wk)
    
    %FEM solution
    %Deformed nodes
    def_nodes = zeros(size(coord));
    %Displacement magnitude
    disp_norms = zeros(length(def_nodes), 1);
    %Displacements
    x = zeros(length(def_nodes), 3);
    for i = 1:length(coord)
        for k = 1:3
            m = 6*(i-1) + k;
            def_nodes(i,k) = coord(i,k) + u(m);
            x(i,k) = u(m);
        end
        disp_norms(i) = norm(x(i,:));
    end
    %FEM result
    wf = max((x(:,3)));
    fprintf("Current deflection: %e\n", wf);
    
    %Absolute dimensionless difference between the two solutions
    diff(end+1) = abs((wf-wk)/wk);
    fprintf("Difference = %e\n", diff(end));
    
    save("DataNNodes="+nn+".mat", "coord", "connec", "u");
end

%Plotting the convergence plot
f1 = figure(1);
hold on;
plot(nnodes, diff, "-ok");
xlabel('Number of nodes');
ylabel('$e = \frac{w_{max, FEM} - w_{max, theory}}{w_{max, theory}}$', 'interpreter', 'latex');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
grid on;
f1.Units = "centimeters";
f1.Position = [1 1 9.81, 8.01];
figureSaver(f1, "Convergence.pdf");

%Order of convergence for the error
p = polyfit(log(nnodes), log(diff), 1);
fprintf("Order of convergence = %e. You should expect a convergence of order -1.\n", p(1));
title("Order of convergence = " + p(1));

%Plotting the 3D deformation
f2 = figure(2);
trisurf(connec, def_nodes(:, 1), def_nodes(:, 2), def_nodes(:, 3), disp_norms, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
axis('tight', 'equal');
f2.Units = "centimeters";
f2.Position = [1 1 8.81, 7.01];
figureSaver(f2, "clampedDeformation.png");
