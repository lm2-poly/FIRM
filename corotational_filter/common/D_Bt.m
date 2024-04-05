% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [K_Bt]=D_Bt(v, theta_i, n)

% This function computes the tensor:
% d_(B^t*v)_d_a, with constant v
%
% input:
%
% v(6nx1)                   constant vector
%
% theta_i(3xn)              core nodal rotation vectors (passed to the core element)
%                           theta_i(:,i) is relevant to node i
%
% n                         number of nodes
%
% output:
%
% K_Bt(6nx6n)               tensor d_(B^t*v)_d_a, with constant v
%
% adapted from:
% B. Nour-Omid and C.C. Rankin
% Finite rotation analysis and consistent linearization using projectors
% Computer Methods in Applied Mechanics and Engineering 93 (1991) 353-384
% http://dx.doi.org/10.1016/0045-7825(91)90248-5
% and
% C.A. Felippa,  B. Haugen
% A unified formulation of small-strain corotational finite elements: I. Theory
% Comput. Methods Appl. Mech. Engrg. 194 (2005) 2285-2335
% doi: 10.1016/j.cma.2004.07.035
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

% indices of rotational dofs
[~, ind_rot]=indices(n);

% K_Bt=diag(0,K_1,...,0,K_n): initialize
K_Bt=zeros(6*n);

for i=1:n
    % rotation dofs of node i
    ind_rot_i=ind_rot(:,i);
    % relevant components of constant vector v
    m=v(ind_rot_i);
    spin_m=spin(m);
    
    % extract code rotation vector
    theta=theta_i(:,i);
    norm_theta=norm(theta);
    spin_theta=spin(theta);
    spin_theta_sq=spin_theta*spin_theta;
    % Nour-Omid, Rankin, eq. (A.14); Felippa, Haugen, eq. (101)
    % equation (72)
    % eta(x)=(1-x/2*cot(x/2))/x^2;
    eta=evaluate('eta',norm_theta);
    B_i=eye(3)-spin_theta/2+eta*spin_theta_sq;
    
    % Nour-Omid, Rankin, eq. (A.17);  Felippa, Haugen, eq. (102)
    % equation (B.5)
    % mu(x)=1/theta*d_eta/d_x=(x^2+4*cos(x)+x*sin(x)-4)/(4*x^4*(sin(x/2))^2)
    mu=evaluate('mu',norm_theta);
    % useful products
    m_tensor_theta=m*theta';
    theta_dot_m=theta'*m;
    % equations (B.3) and (B.4)
    K_Bt(ind_rot_i,ind_rot_i)= ( ...
        eta*(theta_dot_m*eye(3)+m_tensor_theta'-2*m_tensor_theta) ...
        +mu*spin_theta_sq*m_tensor_theta-spin_m/2 ...
        )*B_i;
    
end
end
