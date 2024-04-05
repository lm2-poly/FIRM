%Nonlinear solver: Arc length's implementation
%See C. Kadapa, « A simple extrapolated predictor for overcoming the 
% starting and tracking issues in the arc-length method for nonlinear 
% structural mechanics », Engineering Structures, vol. 234, p. 111755, 
% may 2021, doi: 10.1016/j.engstruct.2020.111755.
%Also based on the supplementaty material given from the previous article.
%However, this is a separate implementation to allow for follower loads.
%Author: Danick Lamoureux
%Creation date: 2023-10-22

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function arclength(coord, connec, dofTable, dispinit, ...
                        loadIncr, d_load, f_load, bc, E, nu, h, psi, nIncr,...
                        maxiter, eps, theta, display_conv, display)
    %Inputs
    %coord: Node coordinates
    %connec: Connectivity table
    %dofTable: Degrees of freedom table
    %dispinit: Initial displacement guess
    %loadIncr: Initial load increment (lambda0)
    %d_load: Dead loads
    %f_load: Follower loads
    %bc: Boundar condition function
    %E: Young's modulus
    %nu: Poisson's ratio
    %h: Thickness
    %psi: Arc-length shape parameter
    %nIncr: Number of allowed increments
    %maxiter: Maximum number of iterations per increment
    %eps: Residual error allowed
    %theta: Relaxation factor (typically 1)
    %display_conv: True or false, to display the residual or not
    %display: To display the deformed shape or not
    %Output: None, saves the data to be read in post-processing
    
    %Clear display if required
    if display_conv
        figure(2);
        clf(2);
    end
    if display
        figure(3);
        clf(3);
    end

    %Displacement containers
    u = dispinit;
    uPrev = u;
    uPrev2 = u;
    uPrev3 = u;
    uPrev4 = u;
    
    %Arc length history
    ds = loadIncr;
    dsPrev = ds;
    dsMax = ds*2;
    dsMin = ds/2;
    
    %Load factor history
    lambda = loadIncr;
    lambdaPrev = 0;
    lambdaPrev2 = 0;
    lambdaPrev3 = 0;
    lambdaPrev4 = 0;
    
    %Convergence history
    converged = false;
    convergedPrev = false;
    
    loadStepConverged = 0;
    
    for n = 1:nIncr
        if display_conv
            figure(1);
            clf(1);
        end
        
        %Predictor step
        if n > 1
            alpha = ds/dsPrev;
            lambda = (1+alpha)*lambdaPrev - alpha*lambdaPrev2;
            u = (1+alpha)*uPrev - alpha*uPrev2;
        end
        Du = u - uPrev;
        Dlambda = lambda - lambdaPrev;
        
        fprintf("Increment = %i of %i\t Estimated load factor = %f\n", n, nIncr, lambda);
        
        convergedPrev = converged;
        converged = false;
        
        ites = [0];
        residuals = [realmax];
		incrcounter = 0;
		
		%Corrector step
        for k = 1:maxiter
            [KT, r, Fext] = assemble(coord, connec, dofTable, u, lambda*d_load(u), lambda*f_load(u, n), E, nu, h);
            [~, Fext] = bc(KT, Fext);
            Fext = Fext/lambda;
            [KT, r] = bc(KT, r);
            
            FtF = Fext'*Fext;
            
            %Calculating A, a and b
            if n == 1
                A = 0.0;
                a = 0.0*Du;
                b = 1.0;
            else
                A = Du'*Du + psi*Dlambda*Dlambda*FtF - ds*ds;
                a = 2*Du;
                b = 2*psi*Dlambda*FtF;
            end
            
            %Residual norm
            norm_r = norm(r, 2);
            residual = sqrt(norm_r^2 + A^2)/(length(r)+1);
            fprintf("Iteration = %i of %i\t Residual = %e\n", k, maxiter, residual);
            if display_conv
                ites(end+1) = k;
                residuals(end+1) = residual;
                figure(1);
                hold on;
                plot(ites(2:end), residuals(2:end), "-k");
                xlabel('Iterations');
                ylabel("Residual");
                grid on;
                set(gca, 'YScale', 'log')
            end
            if residual < eps && k > 3
                converged = true;
                break;
            elseif residual > 1E5*residuals(end-1) || residual > 100
                break;
            end
            
            if residual > residuals(end-1)
                incrcounter = incrcounter + 1;
            end
            
            %Calculating du
            du1 = KT\Fext;
            du2 = KT\r;
            
            dlambda = (a'*du2 - A)/(b + a'*du1);
            du = -du2 + dlambda*du1;
            
            Du = Du + theta*du;
            Dlambda = Dlambda + theta*dlambda;
            u = u + theta*du;
            lambda = lambda + theta*dlambda;
        end
        
        if converged == false && ds == dsMin
            error("Could not converge. Reduce the minimum arc length or change initial load factor");
        end
        
        %Solution update
        if converged
            if n == 1
                ds = sqrt(Du'*Du + psi*Dlambda*Dlambda*FtF);
                dsMax = 10*ds;
                dsMin = ds/4096;
            end
            lambdaPrev4 = lambdaPrev3;
            lambdaPrev3 = lambdaPrev2;
            lambdaPrev2 = lambdaPrev;
            lambdaPrev = lambda;
            
            dsPrev = ds;
            
            uPrev4 = uPrev3;
            uPrev3 = uPrev2;
            uPrev2 = uPrev;
            uPrev = u;
            
            if convergedPrev
                if ds == dsMax
                    dsMax = 1.1*dsMax;
                end
                ds = min(max(1.5*ds, dsMin), dsMax);
            end
            
            if display_conv
                figure(2);
                hold on;
                scatter(n, lambda, 'ok');
                xlabel('Increments');
                ylabel("Load factor");
                grid on;
            end
            if display
                figure(3);
                [x, disp_norms] = get_displacements(u);
                def_nodes = get_deformed_nodes(coord, x);
                display_deformation(connec, def_nodes, disp_norms, true);
            end
            
            loadStepConverged = loadStepConverged + 1;
            save("arclength/"+num2str(n)+".mat", "coord", "connec", "u", "lambda", "Fext", "E", "nu", "h", "r");
            
            
            [x, ~] = get_displacements(u); 
            [f0, ~] = get_displacements(Fext); 
            [r0, ~] = get_displacements(r); 
            [nodes_angles, ~] = get_node_angles(def_nodes, connec);
            vtkwrite("arclength/vtk/solution_"+num2str(n)+".vtk", 'polydata', 'triangle', coord(:,1), coord(:,2), coord(:,3), connec, ...
                    'vectors', 'displacement', x(:,1), x(:,2), x(:,3), 'vectors', 'rotations', x(:,4), x(:,5), x(:,6), ...
                    'vectors', 'force', f0(:,1), f0(:,2), f0(:,3), 'vectors', 'moment', f0(:,4), f0(:,5), f0(:,6),...
                    'vectors', 'residualLin', r0(:,1), r0(:,2), r0(:,3), 'vectors', 'residualRot', r0(:,4), r0(:,5), r0(:,6),...
                    'scalars', 'angles', nodes_angles);
        else            
            if convergedPrev
                ds = max(ds/2, dsMin);
            else
                    ds = max(ds/4, dsMin);
            end
        end
        
    end
end