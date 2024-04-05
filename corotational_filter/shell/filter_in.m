% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % filtered nodal parameters (to be passed to the core element)
    bar_a, ...
    ... % tensor relating delta_tld_a to delta_a
    M, ...
    ... % tensors needed to filter out the rigid rototranslation t,hat_R
    t, hat_R, T, hat_G, hat_P, ...
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
    d_load, f_load)

% This routine implements the closed-form formulas derived in Section 2.
% It is called by the assembler before the core-element routine,
% and has the purpose of removing rigid body motions.
% In particular, it transforms the element nodal parameters a
% into their filtered counterpart bar_a,
% and the dead [respectively, follower] distributed loads into their pulled-back counterparts
% R'l_d [respectively, l_f*area_def/area, to be passed to the core element.
%
% input:
%
% coor(3x3)                 nodal coordinates of reference triangle
%                           coor(:,i) contains the coordinates of node V_i
%                           the third row of coor (z-coordinates) in zero
%
% a(18x1)                   nodal parameters
%                           a=(u_1;theta_1;u_2;theta_2;u_3;theta_3)
%
% e_def(3x1), h_def(3x1), n_def(3x1) frame attached to deformed triangle
%
% area, l23, l31, l12       reference area and side lengths
%
% area_def, l23_def, l31_def, l12_def	deformed area and side lengths
%
% d_load(6x3), f_load(6x3)  dead and follower loads, linearly varying over the triangle
%                           d_load(:,i) and f_load(:,i) are values at node V_i
%
% output:
%
% bar_a(18x1)               filtered nodal parameters
%                           a=(bar_u_1;bar_theta_1;bar_u_2;bar_theta_2;bar_u_3;bar_theta_3)
%
% M(18x18)                  tensor relating delta_tld_a to delta_a
%                           delta_tld_a = M*delta_a
%
% t(3x1), hat_R(3x3)        rigid rototranslation (t,hat_R)
%
% T(3x18), hat_G(3x18), hat_P(18x18) tensors needed to filter out the rigid rototranslation t,hat_R
%                           delta_hat_a = hat_P*shp_hat_R'*delta_tld_a
%
% chk_R(3x3)                rigid rotation (chk_R)
%
% chk_G(3x18), chk_P(18x18) tensors needed to filter out the rigid rotation chk_R
%                           delta_chk_a = chk_P*shp_chk_R'*delta_hat_a
%
% B(18x18)                  tensor relating delta_bar_a to delta_chk_a
%                           delta_bar_a = B*delta_chk_a
%
% b(1,18)                   derivative of deformed area with respect to bar_a, divided by deformed area
%
% n_load, load(6,3,n_load)  number of load conditions (n_load=2) and loads
%                           load(:,i,l) are intensities at node V_i for load condition l
%                           l=1: dead loads; l=2: follower loads
%
% Note: all the above quantities are expressed in the element reference frame (G,e,h,n)
%
% This code is part of a Matlab toolkit distributed as supplementary material of the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528
%
% Authors' e-mail addresses:
% caselli@ing.uniroma2.it (Federica Caselli)
% bisegna@uniroma2.it (Paolo Bisegna)
%
% (C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% number of element nodes
n_el_n = 3;

% indices of translational and rotational dofs
[ind_transl, ind_rot]=indices(n_el_n);
% nodal displacements
% u_i(:,i) is the displacement of node i
u_i=a(ind_transl);
% nodal rotation vectors
% theta_i(:,i) is the rotation vector of node i
theta_i=a(ind_rot);

% trasformation from nodal rotation vectors to nodal rotation tensors
% and from variations of nodal rotation vectors to nodal spins
[ ... % nodal rotation tensors
    R_i, ...
    ... % tensor relating delta_tld_a to delta_a
    M]=rot2spin( ...
    ... % nodal rotation vectors
    theta_i, ...
    ... % number of nodes
    n_el_n);

% filtering out the rigid rototranslation t,hat_R
[ ... % tensors needed to filter out the rigid rototranslation t,hat_R
    t, hat_R, T, hat_G, hat_P, ...
    ... % nodal displacements and rotation tensors after filtering out the rigid rototranslation t,hat_R
    hat_u_i, hat_R_i]=hat( ...
    ... % nodal coordinates of reference triangle
    coor, ...
    ... % nodal displacements and rotation tensors
    u_i, R_i, ...
    ... % area and side lengths of deformed triangle
    area_def, l23_def, l31_def, l12_def, ...
    ... % frame attached to deformed triangle
    e_def, h_def, n_def);

% filtering out the rigid rotation chk_R
[ ... % tensors needed to filter out the rigid rotation chk_R
    chk_R, chk_G, chk_P, ...
    ... % core nodal displacements and rotation tensors
    bar_u_i, bar_R_i]=chk( ...
    ... % nodal coordinates of reference triangle
    coor, ...
    ... % nodal displacements and rotation tensors after filtering out the rigid rototranslation t,hat_R
    hat_u_i, hat_R_i, ...
    ... % area and side lengths of reference triangle
    area, l23, l31, l12, ...
    ... % area and side lengths of deformed triangle
    area_def, l23_def, l31_def, l12_def);

% transformation from filtered nodal rotation tensors to filtered nodal rotation vectors
% and from filtered nodal spins to variations of filtered nodal rotation vectors
[ ... % core nodal rotation vectors
    bar_theta_i, ...
    ... % tensor relating delta_bar_a to delta_chk_a
    B]=spin2rot( ...
    ... % core nodal rotation tensors
    bar_R_i, ...
    ... % number of nodes
    n_el_n);

% filtered nodal parameters
bar_a=[bar_u_i; bar_theta_i];
% arrange into a 18x1 vector
bar_a=bar_a(:);

% derivative of deformed area with respect to bar_a, divided by deformed area
[b]=b_fun(coor, bar_u_i, area_def);

% overall rotation tensor
R=hat_R*chk_R;

% two load conditions: dead loads and follower loads
n_load=2;
% loads enter equations (91), (92), (93), (105), (113)
% dead loads are acted upon by Rt
load(:,:,1)=blkdiag(R',R')*d_load;
% follower loads are rescaled for area deformation
load(:,:,2)=f_load*area_def/area;

end
