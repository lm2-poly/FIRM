%Nonlinear solver: Newton's implementation
%Convergence according to the displacement and/or the reconfiguration could
%be added, but it was found that the residual was adequate
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

function [u, f, r] = newton(coord, connec, dofTable, u, d_load, f_load, bc, E, nu, h, eps, N, display_conv, display, theta, variable_relax)
    
    %Inputs:
    %   coord: Nodes coordinates matrix
    %   connec: Triangulation
    %   dofTable: Degree of freedom table
    %   u: Initial estimate of the displacement
    %   d_load: Dead loads matrix (functions of u)
    %   f_load: Follower load matrix (functions of u)
    %   bc: Boundary condition function
    %   E: Young's modulus
    %   nu: Poisson's ratio
    %   h: Plate's thickness
    %   eps: Tolerance for the residual vector
    %   N: Number of allowed iterations
    %   display_conv: Display the convergence plot
    %   display: To display or not the deformed nodes
    %   theta: Relaxation parameter
    %   variable_relax: To activate variable relaxation
    %
    %Outpus:
    %   u: Final displacement field
    %   f: Final internal load field
    %   r: Final residual
    
    %Clear display if required
    if display_conv
        figure(1);
        clf(1);
    end
    if display
        figure(2);
        clf(2);
    end
    
    %Displacement variation container
    var_r = realmax;
    
    %Number of iterations
    ite = 0;
    
    ites = [];
    residuals = [];
    
    theta0 = theta;    
    
    [du, r, f] = solve(coord, connec, dofTable, u, d_load, @(u) f_load(u), bc, E, nu, h);
    
    while (var_r > eps && ite < N) || ite < 2        
        if variable_relax
            %Change relaxation if required
            theta = theta0;
            [du1, r1, f1] = solve(coord, connec, dofTable, u+theta*du, d_load, @(u) f_load(u), bc, E, nu, h);
            var_r1 = norm(r1)/length(r1);
            i = 0;
            while var_r1 > var_r && i < 6
                theta = theta/2;
                [du1, r1, f1] = solve(coord, connec, dofTable, u+theta*du, d_load, @(u) f_load(u), bc, E, nu, h);
                var_r1 = norm(r1)/length(r1);
                i = i + 1;
                fprintf("Relaxation %i\t Residual = %f\n", i, var_r1);
            end
            var_r1 = norm(r1)/length(r1);
        end
        
        %Adding the variation to the previously obtained displacement
        u = u + theta*du;
        
        [du, r, f] = solve(coord, connec, dofTable, u, d_load, @(u) f_load(u), bc, E, nu, h);
        
        %Variation
        var_r = norm(r)/length(r);
        
        %Incrementing
        ite = ite + 1;

        fprintf("Iteration = %i\t Relaxation = %f\t Norm of the residual %e\n", ite, theta, var_r);
        
        %Display
        if display_conv
            ites(end+1) = ite;
            residuals(end+1) = var_r;
            figure(1);
            clf(1);
            hold on;
            plot(ites, residuals, '-k');
            xlabel('Iterations');
            grid on;
            set(gca, 'YScale', 'log')
        end
        if display
            figure(2);
            [x, disp_norms] = get_displacements(u);
            def_nodes = get_deformed_nodes(coord, x);
            display_deformation(connec, def_nodes, disp_norms, true);
            view(-45, 30);
        end
    end
end