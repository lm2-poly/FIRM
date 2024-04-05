%Assemble the stiffness matrix and load vector
%This is heavily based on an anonymous contributor's code for the finite element assembly
%Editor: Danick Lamoureux
%Modification date: 2023-06-14

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [ K, r, f ] = assemble( coord, connec, dofTable, u, d_load, f_load, E, nu, h )
    % INPUT:
    %  coord: Node matrix
    %  connec: Triangulation
    %  dofTable: Degree of freedom table
    %  u: Initial deformation estimate
    %  d_load: Dead loads matrix
    %  f_load: Follower load matrix
    %  E: Young's modulus
    %  nu: Poisson's ratio
    %  h: Plate's thickness

    % OUTPUT:
    %  K: tangent stiffness matrix; (ndof,ndof)
    %  r: residual vector; (ndof,1)
    %  f: Load vector; (ndof, 1) - MODIFIED

    ne = size(connec, 1);
    ndof = length(u);

    % Initialize
    lgn = zeros(18*18*ne,1);
    col = zeros(18*18*ne,1);
    val = zeros(18*18*ne,1);
    r = zeros(ndof, 1);
    f = zeros(ndof, 1); %-ADDED

    % ...
    v = 1:18;
    w = 1:3;
    kk = 1:1:(18*18);
    for eid = 1:ne

        % Unpack
        current_connec = connec(eid,:);
        current_coord = coord(current_connec,:);
        current_dof = dofTable(:,current_connec);
        current_u = u(current_dof(:));
        current_d_load = d_load(current_connec,:);
        current_f_load = f_load(current_connec,:);
        % Contribution of current element
        
        [re, Ke, fe] = corot_filter(current_coord', current_u, current_d_load', current_f_load', E, nu, h);

        % Add contributions to the residual vector and load vector
        r(current_dof(:)) = r(current_dof(:)) + re;
        f(current_dof(:)) = f(current_dof(:)) + fe; %-ADDED
        % ---
        ii = repmat(current_dof(:), 1, 18);
        jj = ii';
        lgn(kk) = ii(:);
        col(kk) = jj(:);
        val(kk) = Ke(:);

        % ...
        v = v + 18;
        w = w + 3;
        kk = kk + (18*18);

    end

    % Form sparse tangent stiffness matrix
    K = sparse(lgn, col, val, ndof, ndof, 18*18*ne);

end