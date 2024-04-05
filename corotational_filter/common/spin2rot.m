% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % core nodal rotation vectors
    bar_theta_i, ...
    ... % tensor relating delta_bar_a to delta_chk_a
    B]=spin2rot( ...
    ... % core nodal rotation tensors
    bar_R_i, ...
    ... % number of nodes
    n)

% This function transforms filtered nodal rotation tensors into filtered nodal rotation vectors 
% and filtered nodal spins into variations of filtered nodal rotation vectors
% the core nodal displacements are unaffected
%
% input:
%
% bar_R_i(3x3xn)            core nodal rotation tensors
%                           bar_R_i(:,:,i) is relevant to node i
%
% n                         number of nodes
%
%
% output
%
% bar_theta_i(3xn)          core nodal rotation vectors
%                           (they are parameters passed to the core element)
%                           bar_tehta_i(:,i) is relevant to node i
%
% B(6nx6n)                  tensor relating delta_bar_a to delta_chk_a;
%                           it acts as identity on translation dofs;
%                           and relates delta_bar_theta_i to delta_chk_omega_i (core nodal spin) on rotational dofs
%                           i.e., delta_bar_theta_i = B_i*delta_chk_omega_i, with delta_chk_omega_i=ax(delta_bar_R_i*bar_R_i^t)
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

% initialize filtered nodal rotation vectors
bar_theta_i=zeros(3,n);

% B=diag(I,B_1,...,I,B_n): initialize
B=zeros(6*n);

for i=1:n
    % translation dofs of node i
    ind_transl_i=ind_transl(:,i);
    % rotation dofs of node i
    ind_rot_i=ind_rot(:,i);
    
    % B acts as identity on translation dofs
    B(ind_transl_i,ind_transl_i)=eye(3);
    
    % adapted from:
    % B. Nour-Omid and C.C. Rankin
    % Finite rotation analysis and consistent linearization using projectors
    % Computer Methods in Applied Mechanics and Engineering 93 (1991) 353-384
    % doi: http://dx.doi.org/10.1016/0045-7825(91)90248-5
    % and
    % C.A. Felippa,  B. Haugen
    % A unified formulation of small-strain corotational finite elements: I. Theory
    % Comput. Methods Appl. Mech. Engrg. 194 (2005) 2285-2335
    % doi: 10.1016/j.cma.2004.07.035
    
    % extract core rotation vector
    % bar_theta_i is ax(log(bar_R_i))
    % Nour-Omid, Rankin, eq. (A.3); Felippa, Haugen, eq. (98)
    % equation (10)
    bar_theta_i(:,i)=ax_log(bar_R_i(:,:,i));
    % tensor B_i, such that: delta_bar_theta_i = B_i*delta_chk_omega_i, with delta_chk_omega_i=ax(delta_bar_R_i*bar_R_i^t)
    theta=bar_theta_i(:,i);
    norm_theta=norm(theta);
    spin_theta=spin(theta);
    % Nour-Omid, Rankin, eq. (A.14); Felippa, Haugen, eq. (101)
    % equation (72)
    % eta(x)=(1-x/2*cot(x/2))/x^2;
    eta=evaluate('eta',norm_theta);
    B(ind_rot_i,ind_rot_i)= eye(3)-spin_theta/2+eta*spin_theta*spin_theta;
    
end

end
