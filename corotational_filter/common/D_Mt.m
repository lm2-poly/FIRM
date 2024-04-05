% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [K_Mt]=D_Mt(v, theta_i, n)

% This function computes the tensor:
% d_(M^t*v)_d_a, with constant v
%
% input:
%
% v(6nx1)                   constant vector
%
% theta_i(3xn)              nodal rotation vectors
%                           theta_i(:,i) is relevant to node i
%
% n                         number of nodes
%
% output:
%
% K_Mt(6nx6n)               tensor d_(M^t*v)_d_a, with constant v
%
% adapted from:
% Costin Pacoste,
% Co-rotational flat facet triangular elements for shell instability analyses
% Comput. Methods Appl. Mech. Engrg. 156 (1998) 75-110
% doi: dx.doi.org/10.1016/S0045-7825(98)80004-2
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

% K_Mt=diag(0,K_1,...,0,K_n): initialize
K_Mt=zeros(6*n);

for i=1:n
    % rotation dofs of node i
    ind_rot_i=ind_rot(:,i);
    % relevant components of constant vector v
    m=v(ind_rot_i);
    spin_m=spin(m);
    
    % extract rotation vector
    theta=theta_i(:,i);
    norm_theta=norm(theta);
    % Pacoste, eq. (74)
    % singular terms around norm_theta=0 have been rearranged in order to
    % obtain functions with removable singularities, which in turn are computed
    % via Taylor expansions around zero
    % equation (B.2)
    c1=evaluate('(1-cos(x))/x^2',norm_theta);
    c2=evaluate('(x-sin(x))/x^3',norm_theta);
    c3=evaluate('(x*sin(x)+2*cos(x)-2)/x^4',norm_theta);
    c4=evaluate('(3*sin(x)-x*cos(x)-2*x)/x^5',norm_theta);
    % useful products
    m_tensor_theta=m*theta';
    theta_tensor_theta=theta*theta';
    theta_dot_m=theta'*m;
    % equations (B.1) and (B.2)
    K_Mt(ind_rot_i,ind_rot_i)= ...
        +c1*(spin_m-m_tensor_theta) ...
        +c2*((m_tensor_theta+m_tensor_theta')+theta_dot_m*eye(3)) ...
        +c3*(spin_m*theta_tensor_theta) ...
        +c4*(theta_dot_m*theta_tensor_theta);
    
end

end
