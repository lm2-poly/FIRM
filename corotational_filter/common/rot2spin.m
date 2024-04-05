% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % nodal rotation tensors
    R_i, ...
    ... % tensor relating delta_tld_a to delta_a
    M]=rot2spin( ...
    ... % nodal rotation vectors
    theta_i, ...
    ... % number of nodes
    n)

% This function transforms nodal rotation vectors into nodal rotation tensors
% and variations of nodal rotation vectors into nodal spins
% the nodal displacements are unaffected
%
% input:
%
% theta_i(3xn)              nodal rotation vectors
%                           theta_i(:,i) is relevant to node i
%
% n                         number of nodes
%
% output:
%
% R_i(3x3xn)                nodal rotation tensors
%                           R_i(:,:,i) is relevant to node i
%
% M(6nx6n)                  tensor relating delta_tld_a to delta_a;
%                           it acts as identity on translation dofs;
%                           and relates delta_omega_i (nodal spin) to delta_theta_i on rotational dofs
%                           i.e., delta_omega_i = M_i*delta_theta_i, with d_omega_i=ax(delta_R_i*R_i^t)
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

% indices of translational and rotational dofs
[ind_transl, ind_rot]=indices(n);

% initialize nodal rotation tensor
R_i=zeros(3,3,n);

% M=diag(I,M_1,...,I,M_n): initialize
M=zeros(6*n);

for i=1:n
    % translation dofs of node i
    ind_transl_i=ind_transl(:,i);
    % rotation dofs of node i
    ind_rot_i=ind_rot(:,i);
    
    % M acts as identity on translation dofs
    M(ind_transl_i,ind_transl_i)=eye(3);
    
    % adapted from:
    % B. Nour-Omid and C.C. Rankin
    % Finite rotation analysis and consistent linearization using projectors
    % Computer Methods in Applied Mechanics and Engineering 93 (1991) 353-384
    % doi: http://dx.doi.org/10.1016/0045-7825(91)90248-5
    % and
    % Costin Pacoste,
    % Co-rotational flat facet triangular elements for shell instability analyses
    % Comput. Methods Appl. Mech. Engrg. 156 (1998) 75-110
    % doi: 10.1016/S0045-7825(98)80004-2
    % and
    % C.A. Felippa,  B. Haugen
    % A unified formulation of small-strain corotational finite elements: I. Theory
    % Comput. Methods Appl. Mech. Engrg. 194 (2005) 2285-2335
    % doi: 10.1016/j.cma.2004.07.035
    
    % extract rotation vector
    theta=theta_i(:,i);
    % rotation tensor, given by expm(spin(theta));
    % Nour-Omid, Rankin, eq. (A.1); Felippa, Haugen, eq. (89)
    % equation (1)
    R_i(:,:,i)=expm_spin(theta);

    % tensor M, such that: delta_tld_a = M*delta_a, equation (35)
    % Nour-Omid, Rankin, eq. (A.11); Felippa, Haugen, eq. (100)
    % cf. also Pacoste, eq. (18)
    norm_theta=norm(theta);
    spin_theta=spin(theta);
    spin_theta_sq=spin_theta*spin_theta;
    c2=evaluate('(1-cos(x))/x^2',norm_theta);
    c3=evaluate('(x-sin(x))/x^3',norm_theta);
    % equation (33)
    M(ind_rot_i,ind_rot_i)=eye(3)+c2*spin_theta+c3*spin_theta_sq;
    
end

end
