%This applies the corotational filter to the element assembly
%This is based on an anonymous contributor's code as well as Caselli and Bisegna
%(2013)'s corotational formulation.
%Currently only supports linear isotropic materials and we assume 3 gauss
%points for the integration
%Editor: Danick Lamoureux
%Modification date: 2023-05-20

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function [ q, K,q_e ] = corot_filter( coord, a, d_load, f_load, E, nu, h )
    % INPUT:
    %  coord: nodes of the current element
    %  a: value of dofs (in global frame); (18,1)
    %  d_load: Dead load matrix
    %  f_load: Follower load matrix
    %  E: Young's modulus
    %  nu: Poisson's ratio
    %  h: Plate's thickness

    % OUTPUT:
    %  q: nodal residual vector (in global frame); (18,1)
    %  K: consistent tangent stiffness matrix (in global frame); (18,18)
    %  q_e: Applied external forces; (18,1) - MODIFIED

    %This is essentially the same content as what is presented in
    %corotationalImplementationTest.m with the add-ons for computing the
    %loads
    
    %Material
    % type of constitutive law
    type='linear_elastic_isotropic_plane_stress';
    % constitutive parameters
    param=struct('E', E, 'nu',nu);
    % constitutive properties
    constitutive=struct('type',type, 'param',param);
    mat_el=struct('n_g',3, 'constitutive',constitutive, 'th',h);

    % Geometry
    [ ...
        ... % reference triangle: centroid G and element reference frame (e,h,n) (in global frame)
        ~, ref2glob, ...
        ... % reference triangle: nodal coordinates (in element reference frame), area, side lengths
        coor, area, l23, l31, l12, ...
        ... % deformed triangle: centroid G_def and attatched frame (e_def,h_def,n_def) (in global frame)
        ~, ~, ...
        ... % deformed triangle: nodal coordinates (in frame attatched to deformed triangle), area, side lengths
        ~, area_def, l23_def, l31_def, l12_def, ...
        ... % deformed triangle: attatched frame (in element reference frame)
        e_def, h_def, n_def, ...
        ... % nodal parameters (in element reference frame)
        a, ...
        ... % dead and follower loads (in element reference frame)
        d_load, f_load]=geometry( ...
        ... % nodal coordinates and nodal parameters (in global frame)
        coord, a, ...
        ... % dead and follower loads (the former in global frame; the latter in element reference frame)
        d_load, f_load);

    % Filter in
    [ ... % filtered nodal parameters (to be passed to the core element)
        bar_a, ...
        ... % tensor relating delta_tld_a to delta_a
        M, ...
        ... % tensors needed to filter out the rigid rototranslation t,hat_R
        ~, hat_R, T, hat_G, hat_P, ...
        ... % tensors needed to filter out the rigid rotation chk_R
        chk_R, chk_G, chk_P, ...
        ... % tensor relating delta_bar_a to delta_chk_a
        B, ...
        ... % derivative of deformed area with respect to bar_a, divided by deformed area
        b, ...
        ... % number of load conditions and loads
        n_load, load]=filter_in( ...
        ... % nodal coordinates of reference triangle
        coor, ...
        ... % nodal parameters
        a, ...
        ... % frame attached to deformed triangle
        e_def, h_def, n_def, ...
        ... % area and side lengths of reference triangle
        area, l23, l31, l12, ...
        ... % area and side lengths of deformed triangle
        area_def, l23_def, l31_def, l12_def, ...
        ... % dead and follower loads
        d_load, f_load);

    [ ... % material tangent stiffness matrix
        bar_K_bar_q_i, ...
        ... % nodal forces equivalent to internal forces
        bar_q_i, ...
        ... % opposite of nodal forces equivalent to applied loads
        bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, ...
        ... % resultants of applied loads
        f_load, m_G_load, M_G_load]= DKT_OPT_shell( ...
        ... % nodal coordinates and nodal parameters
        coor, bar_a, ...
        ... % number of load conditions and applied loads
        n_load, load, ...
        ... % material parameters
        mat_el);

    % Filter out
    [ ... % nodal residual vector and consistent tangent stiffness tensor
        q, K] = filter_out( ...
        ... % nodal parameters and filtered nodal parameters
        a, bar_a, ...
        ... % tensor relating delta_tld_a to delta_a
        M, ...
        ... % tensors needed to filter out the rigid rototranslation t,hat_R
        hat_R, T, hat_G, hat_P, ...
        ... % tensors needed to filter out the rigid rotation chk_R
        chk_R, chk_G, chk_P, ...
        ... % tensor relating delta_bar_a to delta_chk_a
        B, ...
        ... % derivative of deformed area with respect to bar_a, divided by deformed area
        b, ...
        ... % nodal forces equivalent to internal forces and material tangent stiffness tensor
        bar_q_i, bar_K_bar_q_i, ...
        ... % opposite of nodal forces equivalent to applied loads
        bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, ...
        ... % resultants of applied loads
        f_load, m_G_load, M_G_load);

        %Internal forces (See equation 99 from Caselli & Bisegna, 2013)
        q_i = M'*diese(hat_R)*hat_P'*diese(chk_R)*chk_P'*B'*bar_q_i;
        %External efforts if the residual is null (ideal case) - MODIFIED
        q_l = q_i;

        % Transformation from element reference frame to global frame
        ref2glob_diese = diese(ref2glob);
        q = ref2glob_diese*q;
        q_e = ref2glob_diese*q_l; %-MODIFIED
        K = ref2glob_diese*K*ref2glob_diese';

end